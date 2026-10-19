"""Available-evidence profiles; no composite importance score or druggability claim."""
import json
import pandas as pd
from .provenance import write_json


def integrate(hubs,statuses,summaries):
    result=hubs.copy()
    for source in ('wormbase','bgee','blast'):
        entry=next((r for r in statuses if r['source']==source),{'status':'NOT_RUN','reason':'Not selected'})
        result[source+'_status']=entry['status'];result[source+'_status_reason']=entry.get('reason','')
        frame=summaries.get(source)
        if entry['status'] in ('SUCCESS','PARTIAL') and frame is not None:
            frame=frame.drop(columns=['UniProt','Orthogroup',source+'_status'],errors='ignore')
            if frame.WBGene.duplicated().any():raise ValueError('Duplicate candidate summary for '+source)
            result=result.merge(frame,on='WBGene',how='left',validate='one_to_one')
        # Explicit missing values, never fake negative evidence.
        for col in {'wormbase':['wormbase_observed_records','wormbase_not_observed_records','wormbase_essentiality'],
                    'bgee':['bgee_conditions','bgee_present_conditions'],
                    'blast':['blast_hits','best_human_UniProt','best_human_gene','best_evalue','best_query_coverage','best_subject_coverage']}[source]:
            if col not in result:result[col]=pd.NA
    return result


def execute(root,config,run,manifest):
    out=run/'06_target_profiles';out.mkdir()
    hubs=pd.read_csv(run/'03_network_analysis/hubs.tsv',sep='\t',keep_default_na=False)
    status=pd.read_csv(run/'05_biological_evidence/evidence_status.tsv',sep='\t',keep_default_na=False).to_dict('records')
    summaries={}
    for source in ('wormbase','bgee','blast'):
        file=run/'05_biological_evidence'/source/'candidate_summary.tsv'
        if file.exists():summaries[source]=pd.read_csv(file,sep='\t')
    master=integrate(hubs,status,summaries)
    go=pd.read_csv(run/'04_functional_context/community_go_summary.tsv',sep='\t',keep_default_na=False)
    master=master.merge(go[['community','significant_BP_terms','top3_GO_IDs','top3_terms']],on='community',how='left',validate='many_to_one')
    master.to_csv(out/'candidate_target_profiles.tsv',sep='\t',index=False,na_rep='NA')
    completeness=master[['WBGene','UniProt','Orthogroup','wormbase_status','bgee_status','blast_status']].copy()
    completeness.to_csv(out/'evidence_completeness.tsv',sep='\t',index=False)
    net=json.loads((run/'03_network_analysis/summary.json').read_text());gos=json.loads((run/'04_functional_context/summary.json').read_text())
    bg=summaries.get('bgee',pd.DataFrame());bl=summaries.get('blast',pd.DataFrame())
    summary={'candidate_count':len(master),'network':net,'GO':gos,'evidence':status,
             'Bgee_hubs_with_records':int((bg.bgee_conditions>0).sum()) if 'bgee_conditions' in bg else None,
             'BLAST_hubs_with_hits':int((bl.blast_hits>0).sum()) if 'blast_hits' in bl else None,
             'BLAST_hubs_without_hits':int((bl.blast_hits==0).sum()) if 'blast_hits' in bl else None,
             'orthogroup':manifest.get('orthogroup_annotation',{}),'interpretation':'Candidate target profiles, not validated targets or druggability predictions'}
    write_json(out/'analysis_summary.json',summary)
    from .plotting import generate
    from .report import generate as report
    generate(root,config,run,out)
    report(root,config,run,out,summary,master)
    manifest['stages']['06']={'status':'SUCCESS','summary':{'profiles':len(master),'evidence_incomplete':any(r['status']!='SUCCESS' for r in status)}}
    return summary
