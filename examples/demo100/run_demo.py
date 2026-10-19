#!/usr/bin/env python3
"""Run the 100-input demonstration in an isolated local workspace."""
import argparse,json,shutil,subprocess,sys
from pathlib import Path

repo=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(repo/'src'))
from ppi_target.provenance import load_config,write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace',default=str(repo/'local/demo100-reproduction'))
    parser.add_argument('--resources-config',default=str(repo/'config/config.example.yaml'),help='Use a configured pinned-resource file for exact snapshot reuse; otherwise pass --refresh')
    parser.add_argument('--refresh',action='store_true',help='Acquire external sources; results may change with database releases')
    args=parser.parse_args()
    workspace=Path(args.workspace).resolve()
    if workspace==repo:parser.error('Use an isolated workspace, not the repository root')
    cp,base,cfg=load_config(args.resources_config)
    for section,keys in [('biogrid',['tab3','gene_info']),('go',['gaf','obo']),('bgee',['source'])]:
        for key in keys:
            for field in ('path','metadata','retrieval_metadata'):
                if field in cfg[section][key]:cfg[section][key][field]=str((base/cfg[section][key][field]).resolve())
    for key in ('target_gzip','target_headers'):cfg['blast'][key]=str((base/cfg['blast'][key]).resolve())
    here=Path(__file__).resolve().parent
    parameters=json.loads((here/'parameters.json').read_text())
    for section,values in parameters.items():cfg[section].update(values)
    cfg['project_root']='..';cfg['input']={'path':'data/raw/candidates.tsv','header':True}
    cfg['annotation']={'optional':True,'orthogroup_file':'data/raw/orthogroup_wbgene_mapping.tsv'}
    cfg['wormbase']={'deferred':False,'snapshot':json.loads((here/'wormbase_snapshot/manifest.json').read_text())}
    cfg['wormbase']['snapshot']['directory']=str(here/'wormbase_snapshot')
    cfg['demo']={'title':'100-input demonstration','publication_note':'This demonstration contains only a sample of a dataset that will be published in full at a later date.'}
    (workspace/'data/raw').mkdir(parents=True,exist_ok=True)
    shutil.copyfile(here/'input.tsv',workspace/'data/raw/candidates.tsv')
    if not (workspace/'src').exists():(workspace/'src').symlink_to(repo/'src',target_is_directory=True)
    config=workspace/'config/config.yaml';write_json(config,cfg)
    command=[sys.executable,str(repo/'run_pipeline.py'),'--config',str(config)]
    if args.refresh:command.append('--refresh')
    return subprocess.call(command,cwd=repo)

if __name__=='__main__':raise SystemExit(main())
