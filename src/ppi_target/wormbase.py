"""Explicit archived phenotype import; uncovered genes remain missing evidence."""
import json,re,shutil
from pathlib import Path
import pandas as pd
from .provenance import sha256,write_json


def from_snapshot(root,spec,hubs,out):
    directory=root/spec['directory'];raw=out/'raw';raw.mkdir(parents=True,exist_ok=True)
    rows=[];summaries=[];sources=[];missing=[]
    for gene in hubs.WBGene:
        name=gene+'.json';expected=spec['files'].get(name);file=directory/name
        if expected is None or not file.exists():
            missing.append(gene)
            summaries.append({'WBGene':gene,'wormbase_record_status':'not_in_snapshot','wormbase_essentiality':'unknown; missing archived evidence'})
            continue
        if sha256(file)!=expected:raise ValueError('WormBase snapshot checksum mismatch: '+name)
        payload=json.loads(file.read_text())
        if payload.get('name')!=gene:raise ValueError('WormBase snapshot gene identity mismatch: '+name)
        fields=payload.get('fields',{})
        if not any(k in fields for k in ('phenotype','phenotype_not_observed')):raise ValueError('Unknown WormBase phenotype schema: '+name)
        shutil.copyfile(file,raw/name);sources.append({'file':name,'sha256':expected})
        local=[]
        for field,observed in [('phenotype',True),('phenotype_not_observed',False)]:
            items=fields.get(field,{}).get('data') or []
            if not isinstance(items,list):raise ValueError('Invalid phenotype list: '+name)
            for item in items:
                ph=item.get('phenotype') or {};ev=item.get('evidence') or {}
                local.append({'WBGene':gene,'phenotype_id':ph.get('id',''),'phenotype_label':ph.get('label',''),
                    'observed':observed,'evidence_json':json.dumps(ev),'supporting_references':';'.join(sorted(set(re.findall(r'WBPaper\d+',json.dumps(ev)))))})
        rows.extend(local)
        summaries.append({'WBGene':gene,'wormbase_record_status':'archived_evidence_available' if local else 'archived_empty_response',
            'wormbase_observed_records':sum(x['observed'] for x in local),'wormbase_not_observed_records':sum(not x['observed'] for x in local),
            'wormbase_essentiality':'not classified; inspect phenotype and experimental context',
            'wormbase_retrieval':spec['retrieved_at']})
    pd.DataFrame(rows,columns=['WBGene','phenotype_id','phenotype_label','observed','evidence_json','supporting_references']).to_csv(out/'phenotypes.tsv',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(out/'candidate_summary.tsv',sep='\t',index=False,na_rep='NA')
    status='SUCCESS' if not missing else 'PARTIAL'
    reason=f"Archived snapshot: {len(sources)}/{len(hubs)} hubs covered; absent snapshot files are missing evidence"
    write_json(out/'source_manifest.json',{'source':'WormBase','status':status,'retrieval_mode':'explicit_archived_snapshot',
        'retrieved_at':spec['retrieved_at'],'database_version':spec['database_version'],'files':sources,'missing_genes':missing,'reason':reason})
    return {'status':status,'records':len(rows),'reason':reason,'retrieved_at':spec['retrieved_at'],'database_version':spec['database_version']}
