"""Three independent evidence modules; failures never become biological absence."""
import json,re
import pandas as pd
from . import bgee,human_homology
from .sources import fetch,independent_sources
from .provenance import write_json,now


def wormbase(root,config,run,hubs,out):
    if config['wormbase'].get('snapshot'):
        from .wormbase import from_snapshot
        return from_snapshot(root,config['wormbase']['snapshot'],hubs,out)
    if config['wormbase'].get('deferred',True) and config.get('_source')!='wormbase':
        meta={'source':'WormBase','status':'NOT_RUN','reason':'Deferred by user after HTTP 403; no current WormBase evidence generated','retrieved_at':None,'identifiers':hubs.WBGene.tolist()}
        write_json(out/'source_manifest.json',meta)
        return {'status':'NOT_RUN','records':0,'reason':meta['reason']}
    rows=[];summaries=[];errors=[];manifests=[]
    for gene in hubs.WBGene:
        try:
            response,meta=fetch(f'https://rest.wormbase.org/rest/widget/gene/{gene}/phenotype',out/'raw'/gene,[gene]);manifests.append(meta)
            payload=response.json();fields=payload.get('fields')
            if not isinstance(fields,dict) or not any(x in fields for x in ('phenotype','phenotype_not_observed')):raise ValueError('Unrecognized phenotype schema; not an empty biological result')
            local=[]
            for field,observed in [('phenotype',True),('phenotype_not_observed',False)]:
                items=fields.get(field,{}).get('data') or []
                if not isinstance(items,list):raise ValueError('Unexpected phenotype list structure')
                for item in items:
                    ph=item.get('phenotype') or {};ev=item.get('evidence') or {}
                    row={'WBGene':gene,'phenotype_id':ph.get('id',''),'phenotype_label':ph.get('label',''),'observed':observed,
                         'evidence_json':json.dumps(ev),'supporting_references':';'.join(sorted(set(re.findall(r'WBPaper\d+',json.dumps(ev)))))}
                    local.append(row)
            rows+=local
            summaries.append({'WBGene':gene,'wormbase_record_status':'evidence_available' if local else 'no_evidence_returned',
                              'wormbase_observed_records':sum(r['observed'] for r in local),'wormbase_not_observed_records':sum(not r['observed'] for r in local),
                              'wormbase_essentiality':'NA; phenotype category IDs require validation','wormbase_retrieval':str(out.relative_to(root))})
        except Exception as e:
            errors.append({'WBGene':gene,'error':str(e),'error_type':type(e).__name__})
            # Uniform access failure: do not send dozens of known blocked requests.
            if getattr(e,'http_status',None) in (401,403):
                if not summaries:raise
                break
    if not summaries:raise ValueError('All WormBase queries failed; no scientific result generated')
    pd.DataFrame(rows,columns=['WBGene','phenotype_id','phenotype_label','observed','evidence_json','supporting_references']).to_csv(out/'phenotypes.tsv',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(out/'candidate_summary.tsv',sep='\t',index=False)
    state='PARTIAL' if errors or len(summaries)!=len(hubs) else 'SUCCESS'
    write_json(out/'source_manifest.json',{'source':'WormBase','status':state,'retrievals':manifests,'errors':errors})
    if errors:write_json(out/'error.json',{'status':state,'errors':errors})
    return {'status':state,'records':len(rows),'reason':'Phenotype classification unvalidated; no essentiality assertion','retrieved_at':now(),'database_version':'unresolved'}


def execute(root,config,run,manifest):
    out=run/'05_biological_evidence';out.mkdir(exist_ok=True)
    hubs=pd.read_csv(run/'03_network_analysis/hubs.tsv',sep='\t',keep_default_na=False)
    source=config.get('_source','all');names=['wormbase','bgee','blast'] if source=='all' else [source]
    def bgee_call(folder):
        if 'bgee' in config.get('_refresh_sources',[]):
            from .refresh import refresh_evidence
            refresh_evidence(root,config,run,'bgee')
        frame,stats,meta=bgee.execute(root,config,run,hubs)
        return {'status':'SUCCESS','records':stats['detail_rows'],'database_version':meta.get('database_version','unresolved'),'retrieved_at':meta['retrieved_at']}
    names=[name for name in names if name not in config.get('_reuse_sources',[])]
    def blast_call(folder):
        if 'blast' in config.get('_refresh_sources',[]):
            from .refresh import refresh_evidence
            refresh_evidence(root,config,run,'blast')
        return human_homology.execute(root,config,run,hubs,folder)
    callbacks={'wormbase':lambda folder:wormbase(root,config,run,hubs,folder),
               'bgee':bgee_call,'blast':blast_call}
    outcomes=independent_sources(names,callbacks,out,hubs.WBGene)
    # Preserve compatible previously successful sources copied by the runner.
    prior=out/'evidence_status.tsv'
    if prior.exists():
        old=pd.read_csv(prior,sep='\t',keep_default_na=False).to_dict('records')
        outcomes=[r for r in old if r['source'] not in names]+outcomes
    for name in ('wormbase','bgee','blast'):
        if name not in {r['source'] for r in outcomes}:outcomes.append({'source':name,'status':'NOT_RUN','records':0,'reason':'Source not selected'})
    pd.DataFrame(outcomes).to_csv(prior,sep='\t',index=False)
    statuses=[r['status'] for r in outcomes]
    state='SUCCESS' if all(x=='SUCCESS' for x in statuses) else ('PARTIAL' if any(x in ('SUCCESS','PARTIAL') for x in statuses) else 'FAILED')
    summary={'status':state,'sources':outcomes,'candidate_count':len(hubs)}
    write_json(out/'evidence_summary.json',summary)
    # Always export the exact candidate list for deferred/manual evidence collection.
    manual=hubs[['WBGene','UniProt','Orthogroup','Gene_symbol','community','degree','hub_frequency']].copy()
    manual['wormbase_url']=manual.WBGene.map(lambda x:'https://wormbase.org/species/c_elegans/gene/'+x)
    manual.to_csv(out/'wormbase_manual_genes.tsv',sep='\t',index=False)
    template=manual[['WBGene','Gene_symbol']].copy()
    for col in ('source','database_release','retrieved_at','phenotype_id','phenotype_label','observed_or_not_observed','evidence_type_RNAi_mutant','allele_id','reference','source_url','notes'):template[col]=''
    template.to_csv(out/'wormbase_manual_evidence_template.tsv',sep='\t',index=False)
    manifest['stages']['05']={'status':state,'summary':summary}
    return summary
