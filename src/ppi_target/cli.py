"""Six-stage CLI with immutable runs, validated reuse and source-isolated evidence."""
from pathlib import Path
import argparse,hashlib,importlib,importlib.metadata,json,shutil,sys,traceback,uuid
from .provenance import load_config,new_run,write_json,sha256,now
from .annotation import load,propagate

NAMES={'01':'01_interactions','02':'02_network','03':'03_network_analysis','04':'04_functional_context','05':'05_biological_evidence','06':'06_target_profiles'}
MODULES={'01':'biogrid','02':'network','03':'analysis','04':'enrichment','05':'evidence','06':'profiles'}
DEPS={'01':[],'02':['01'],'03':['02'],'04':['02','03'],'05':['03'],'06':['03','04','05']}

def fingerprint(stage,cfg,identity,upstream,root):
    sections={'01':['biogrid'],'02':[],'03':['network'],'04':['go'],'05':['bgee','blast','wormbase'],'06':['plotting']}[stage]
    selected={k:cfg[k] for k in sections}
    # Paths need not match across machines, but actual local snapshot bytes must.
    paths=[]
    if stage=='01':paths=[cfg['biogrid'][x]['path'] for x in ('tab3','gene_info')]
    if stage=='04':paths=[cfg['go'][x]['path'] for x in ('gaf','obo')]
    if stage=='05':paths=[cfg['bgee']['source']['path'],cfg['blast']['target_gzip']]
    if stage=='05' and cfg['wormbase'].get('snapshot'):
        spec=cfg['wormbase']['snapshot']
        paths += [str(Path(spec['directory'])/name) for name in sorted(spec['files'])]
    sources={p:sha256(root/p) if (root/p).exists() else 'missing' for p in paths}
    packages={p:importlib.metadata.version(p) for p in ['networkx','numpy','scipy','pandas','python-louvain']}
    content={'stage':stage,'method_version':'2.0.0','identity':identity,'parameters':selected,'sources':sources,
             'upstream':{k:upstream[k] for k in DEPS[stage]},'packages':packages}
    return hashlib.sha256(json.dumps(content,sort_keys=True).encode()).hexdigest()


def verified_stage(root,previous,stage,key):
    info=previous.get('stages',{}).get(stage,{})
    if info.get('status') not in (('SUCCESS','PARTIAL','FAILED') if stage=='05' else ('SUCCESS','PARTIAL')) or info.get('cache_key')!=key:return False
    prefix=f"results/runs/{previous['run_id']}/{NAMES[stage]}/"
    outputs=[x for x in previous.get('outputs',[]) if x['path'].startswith(prefix)]
    if not outputs:return False
    return all((root/x['path']).is_file() and sha256(root/x['path'])==x['sha256'] for x in outputs)


def stage_arg(value):
    value=str(value).zfill(2)
    if value not in NAMES:raise argparse.ArgumentTypeError('Stage must be 1–6')
    return value


def main(default_stage=None):
    p=argparse.ArgumentParser(description='Physical PPI target profiles; optional OG and independent evidence failures.')
    p.add_argument('--config',default='config/config.yaml');p.add_argument('--stage',type=stage_arg,default=default_stage)
    p.add_argument('--from-stage',type=stage_arg,default='01');p.add_argument('--to-stage',type=stage_arg,default='06')
    p.add_argument('--source',choices=['all','wormbase','bgee','blast'],default='all')
    p.add_argument('--reuse',action='store_true',help='Default: reuse compatible, checksum-validated results')
    p.add_argument('--force',action='store_true',help='Recompute selected local stages; does not refresh external sources')
    p.add_argument('--refresh',action='store_true',help='New external snapshots for selected sources; never overwrite prior successes')
    p.add_argument('--resume',help='Use this run as the compatible source instead of current.json')
    a=p.parse_args();cp,root,cfg=load_config(a.config)
    pins=root/'config/snapshot_pins.local.json'
    if pins.exists():
        for k,v in json.loads(pins.read_text()).items():cfg[k]=v
    rd,manifest=new_run(root,cp,cfg);print('RUN_ID='+manifest['run_id'],flush=True)
    previous={};pointer=root/'results/current.json'
    if a.resume:previous=json.loads((root/'results/runs'/a.resume/'run_manifest.json').read_text())
    elif pointer.exists():previous=json.loads((root/'results/runs'/json.loads(pointer.read_text())['run_id']/'run_manifest.json').read_text())
    manifest['execution_commands']=[{'at':now(),'argv':sys.argv}];manifest['package_versions']={d.metadata['Name']:d.version for d in importlib.metadata.distributions()}
    requested={a.stage} if a.stage else {s for s in NAMES if a.from_stage<=s<=a.to_stage}
    required=set(requested)
    for s in reversed(NAMES):
        if s in required:required.update(DEPS[s])
    history=[previous] if a.resume else []
    if not a.resume:
        for mf in sorted((root/'results/runs').glob('*/run_manifest.json'),reverse=True):
            if mf.parent!=rd:
                try:history.append(json.loads(mf.read_text()))
                except (ValueError,OSError):pass
    cfg['_source']=a.source;cfg['_reuse_sources']=[]
    keys={}
    try:
        rows,canonical,ogqc=load(root,cfg,rd/'annotation')
        canonical.to_csv(rd/'annotation/candidate_annotations.tsv',sep='\t',index=False)
        manifest['orthogroup_annotation']=ogqc;manifest['input_files']=[{'path':str((root/cfg['input']['path']).resolve()),'sha256':sha256(root/cfg['input']['path'])}]
        if a.refresh:
            from .refresh import refresh
            cfg=refresh(root,cfg,requested,a.source,rd)
            manifest['resolved_configuration']=cfg
        for stage in NAMES:
            if stage not in required:
                manifest['stages'][stage]={'status':'NOT_RUN','reason':'Not selected'};continue
            keys[stage]=fingerprint(stage,cfg,ogqc['identity_checksum'],keys,root)
            # Target profiles/rendering include annotations and source status, unlike topology.
            if stage=='06':keys[stage]=hashlib.sha256((keys[stage]+canonical.to_csv(index=False)+json.dumps(manifest['stages'].get('05',{}).get('summary',{}),sort_keys=True)).encode()).hexdigest()
            compatible=[m for m in history if verified_stage(root,m,stage,keys[stage])]
            if stage=='05':
                compatible.sort(key=lambda m:sum(x.get('status')=='SUCCESS' for x in m['stages']['05'].get('summary',{}).get('sources',[])),reverse=True)
            previous=next(iter(compatible),{})
            reusable=bool(previous) and not (a.force and stage in requested) and not (a.refresh and stage in ('01','04','05'))
            dest=rd/NAMES[stage]
            if reusable and stage!='05':
                shutil.copytree(root/'results/runs'/previous['run_id']/NAMES[stage],dest)
                manifest['stages'][stage]={**previous['stages'][stage],'reused_from':previous['run_id']}
                print(f'[Stage {stage}] REUSED validated artifacts',flush=True)
            else:
                if stage not in requested and not reusable:raise ValueError('Missing compatible upstream Stage '+stage+'; run full pipeline or include this stage')
                if stage=='05' and previous:
                    # Preserve successful source snapshots; rerun only requested failed/deferred sources.
                    shutil.copytree(root/'results/runs'/previous['run_id']/NAMES[stage],dest)
                    old=json.loads((dest/'evidence_summary.json').read_text())
                    cfg['_reuse_sources']=[x['source'] for x in old['sources'] if x['status']=='SUCCESS' and x['source'] not in cfg.get('_refresh_sources',[]) and not (a.force and stage in requested and (a.source=='all' or a.source==x['source']))]
                    for name in ('wormbase','bgee','blast'):
                        if name not in cfg['_reuse_sources'] and (dest/name).exists():shutil.rmtree(dest/name)
                code=rd/'software'/NAMES[stage];code.mkdir(parents=True)
                for f in (root/'src/ppi_target').glob('*.py'):shutil.copy2(f,code/f.name)
                module=importlib.import_module('ppi_target.'+MODULES[stage])
                manifest['stages'][stage]={'status':'RUNNING','started_at':now()}
                write_json(rd/'run_manifest.json',manifest)
                module.execute(root,cfg,rd,manifest)
                if a.refresh:
                    keys[stage]=fingerprint(stage,cfg,ogqc['identity_checksum'],keys,root)
                    write_json(pins,{k:cfg[k] for k in ('biogrid','go','bgee','blast')})
                    manifest['resolved_configuration']={k:v for k,v in cfg.items() if not k.startswith('_')}
                manifest['stages'][stage].update(cache_key=keys[stage],finished_at=now())
            if stage=='04':
                renames={'enrichment_all.tsv':'go_enrichment.tsv','community_summary.tsv':'community_go_summary.tsv','background.tsv':'gene_annotation_coverage.tsv'}
                for old,new in renames.items():
                    if (dest/old).exists():(dest/old).rename(dest/new)
                if not (dest/'source_manifest.json').exists():
                    write_json(dest/'source_manifest.json',{'sources':[json.loads(f.read_text()) for f in (dest/'raw').glob('*.source_manifest.json')]})
            propagate(rd,canonical,ogqc)
            write_json(rd/'run_manifest.json',manifest)
        state=manifest['stages'].get('05',{}).get('status')
        manifest['final_status']='PARTIAL' if state in ('PARTIAL','FAILED') else 'SUCCESS'
    except Exception as exc:
        manifest['final_status']='FAILED';manifest['errors'].append({'at':now(),'error_type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()})
        if 'stage' in locals():manifest['stages'][stage]={'status':'FAILED','reason':str(exc)}
        raise
    finally:
        manifest['finished_at']=now()
        manifest['outputs']=[{'path':str(f.relative_to(root)),'sha256':sha256(f)} for f in sorted(rd.rglob('*')) if f.is_file() and f.name!='run_manifest.json' and '__pycache__' not in f.parts]
        write_json(rd/'run_manifest.json',manifest)
        if manifest['final_status']!='FAILED':
            write_json(pointer,{'run_id':manifest['run_id']})
            for stage,name in NAMES.items():
                if (rd/name).exists():
                    link=root/'results'/name
                    if link.is_symlink():link.unlink()
                    if not link.exists():link.symlink_to(Path('runs')/manifest['run_id']/name,target_is_directory=True)
        print('RUN STATUS:',manifest['final_status'],flush=True)
    return 0
