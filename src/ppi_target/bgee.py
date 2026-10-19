"""Per-hub expression evidence from the current project's official Bgee snapshot."""
import csv,gzip
from collections import defaultdict
import pandas as pd
from .enrichment import resource
from .provenance import write_json


def execute(root,config,run,hubs):
    out=run/'05_biological_evidence/bgee';out.mkdir(parents=True,exist_ok=True);raw=out/'raw';raw.mkdir()
    source,meta=resource(root,config['bgee']['source'],raw,'expression.tsv.gz')
    details=[];ids=set(hubs.WBGene)
    with gzip.open(source,'rt') as f:
        reader=csv.DictReader(f,delimiter='\t')
        for row in reader:
            if row['Gene ID'] in ids:details.append(row)
        fields=reader.fieldnames
    df=pd.DataFrame(details,columns=fields);df.to_csv(out/'expression.tsv',sep='\t',index=False)
    rows=[]
    for gene in sorted(ids):
        g=df[df['Gene ID']==gene];present=g[g.Expression=='present'];absent=g[g.Expression=='absent']
        rows.append({'WBGene':gene,'bgee_record_status':'evidence_available' if len(g) else 'no_matching_record',
                     'bgee_conditions':len(g),'bgee_present_conditions':len(present),'bgee_absent_conditions':len(absent),
                     'bgee_present_anatomies':present['Anatomical entity ID'].nunique(),'bgee_present_stages':present['Developmental stage ID'].nunique(),
                     'bgee_gold_present_conditions':int((present['Call quality']=='gold quality').sum()),
                     'bgee_present_anatomy_labels':'; '.join(sorted(present['Anatomical entity name'].unique())),
                     'bgee_present_stage_labels':'; '.join(sorted(present['Developmental stage name'].unique())),
                     'bgee_snapshot':str(out.relative_to(root))})
    summary=pd.DataFrame(rows);summary.to_csv(out/'candidate_summary.tsv',sep='\t',index=False)
    stats={'hub_count':len(ids),'hubs_with_records':int((summary.bgee_conditions>0).sum()),'detail_rows':len(df),'source_version':meta.get('database_version','unresolved'),
           'interpretation':'Condition counts include ontology-propagated calls; not independent experiments or high-expression classifications'}
    write_json(out/'summary.json',stats)
    write_json(out/'source_manifest.json',meta)
    return summary,stats,meta
