"""Optional OG annotation; never part of network identity or topology cache keys."""
import csv,hashlib,json,re
from collections import defaultdict,Counter
from pathlib import Path
import pandas as pd
from .identifiers import WB,UP
from .provenance import sha256,write_json


def load(root,config,out=None):
    path=(root/config['input']['path']).resolve()
    with path.open(newline='') as f:
        data=list(csv.reader(f,delimiter='\t'))
    if config['input'].get('header',True):
        header=data.pop(0) if data else []
        if not {'WBGene','UniProt'}<=set(header):raise ValueError('Input requires WBGene and UniProt columns')
    else:header=['WBGene','UniProt']
    rows=[]
    for i,r in enumerate(data,2 if config['input'].get('header',True) else 1):
        if len(r)!=len(header):raise ValueError(f'Input row {i}: incorrect column count')
        x=dict(zip(header,r));g=x['WBGene'];u=x['UniProt'];og=x.get('Orthogroup','')
        if not WB.fullmatch(g) or not UP.fullmatch(u):raise ValueError(f'Input row {i}: malformed required identifier')
        if og and not re.fullmatch(r'OG\d+',og):raise ValueError(f'Input row {i}: malformed OG')
        rows.append({'source_row':i,'WBGene':g,'UniProt':u,'Orthogroup':og})
    if not rows:raise ValueError('Input is empty')
    genes={x['WBGene'] for x in rows};ogs=defaultdict(set)
    for r in rows:
        if r['Orthogroup']:ogs[r['WBGene']].add(r['Orthogroup'])
    ogpath=root/config.get('annotation',{}).get('orthogroup_file','data/raw/orthogroup_wbgene_mapping.tsv')
    issues=[];loaded=0;seen=set();fatal=[]
    if ogpath.exists():
        with ogpath.open(newline='') as f:
            reader=csv.DictReader(f,delimiter='\t')
            if set(reader.fieldnames or [])!={'Orthogroup','WBGene'}:fatal.append('OG file requires exactly Orthogroup and WBGene columns')
            else:
                for i,r in enumerate(reader,2):
                    g=r['WBGene'];o=r['Orthogroup'];reason=None
                    if not WB.fullmatch(g or '') or not re.fullmatch(r'OG\d+',o or ''):reason='malformed_or_missing_identifier';fatal.append(f'OG row {i}: {reason}')
                    elif (g,o) in seen:reason='duplicate_row'
                    elif g not in genes:reason='WBGene_not_in_canonical_mapping'
                    else:ogs[g].add(o);loaded+=1
                    seen.add((g,o))
                    if reason:issues.append({'row':i,'WBGene':g,'Orthogroup':o,'reason':reason})
    for g,values in ogs.items():
        if len(values)>1:issues.append({'WBGene':g,'Orthogroup':';'.join(sorted(values)),'reason':'multiple_OG_assignments_preserved'})
    for r in rows:r['Orthogroup']=';'.join(sorted(ogs[r['WBGene']]))
    canonical=pd.DataFrame(rows)[['WBGene','UniProt','Orthogroup']].drop_duplicates()
    available=ogpath.exists() or any(ogs.values())
    qc={'orthogroup_annotation':'available' if available else 'not_provided','file_detected':ogpath.exists(),
        'file':str(ogpath) if ogpath.exists() else None,'sha256':sha256(ogpath) if ogpath.exists() else None,
        'mappings_loaded':loaded,'unique_OGs':len(set().union(*ogs.values())) if ogs else 0,
        'WBGene_with_OG':sum(bool(ogs[g]) for g in genes),'WBGene_without_OG':sum(not ogs[g] for g in genes),
        'issues':issues,'errors':fatal,'OG_aware_analyses':'YES' if available and any(ogs.values()) else 'NO',
        'OG_analysis_status':'SUCCESS' if any(ogs.values()) else 'NOT_RUN','reason':'' if any(ogs.values()) else 'orthogroup annotation not provided'}
    identity=sorted({(r['WBGene'],r['UniProt']) for r in rows})
    qc['identity_checksum']=hashlib.sha256(json.dumps(identity).encode()).hexdigest()
    if out:
        out.mkdir(parents=True,exist_ok=True);write_json(out/'orthogroup_qc.json',qc)
        pd.DataFrame(issues,columns=['row','WBGene','Orthogroup','reason']).to_csv(out/'orthogroup_issues.tsv',sep='\t',index=False)
    if fatal:raise ValueError('; '.join(fatal))
    return rows,canonical,qc


def propagate(run,canonical,qc):
    proteins=canonical.groupby('WBGene').UniProt.apply(lambda x:';'.join(sorted(set(x)))).to_dict()
    node_file=run/'02_network/nodes.tsv'
    if node_file.exists():
        node_table=pd.read_csv(node_file,sep='\t',keep_default_na=False)
        proteins={**dict(zip(node_table.WBGene,node_table.UniProt)),**proteins}
    og=canonical.groupby('WBGene').Orthogroup.apply(lambda x:';'.join(sorted({o for v in x for o in str(v).split(';') if o}))).to_dict()
    for stage in ['01_interactions','02_network','03_network_analysis','04_functional_context','05_biological_evidence','06_target_profiles']:
        for f in (run/stage).rglob('*.tsv') if (run/stage).exists() else []:
            if 'raw' in f.parts or 'database' in f.parts:continue
            frame=pd.read_csv(f,sep='\t',keep_default_na=False)
            if 'WBGene' in frame.columns:
                frame['Orthogroup']=frame.WBGene.map(og).fillna('')
                if 'UniProt' not in frame:frame['UniProt']=frame.WBGene.map(proteins).fillna('')
                frame.to_csv(f,sep='\t',index=False)
    # Annotation-only summaries: one gene may appear in multiple OG rows; no topology collapse.
    metrics=run/'03_network_analysis/node_metrics.tsv'
    if metrics.exists():
        df=pd.read_csv(metrics,sep='\t',keep_default_na=False)
        detail=[]
        for r in df.to_dict('records'):
            for o in str(r.get('Orthogroup','')).split(';'):
                if o:detail.append({'Orthogroup':o,'WBGene':r['WBGene'],'community':r['community'],'degree':r['degree'],'robust_hub':r['robust_hub']})
        target=metrics.parent/'og_summary.tsv'
        if detail:
            d=pd.DataFrame(detail)
            d.groupby('Orthogroup').agg(WBGene_count=('WBGene','nunique'),communities=('community','nunique'),robust_hubs=('robust_hub','sum')).reset_index().to_csv(target,sep='\t',index=False)
            d.to_csv(metrics.parent/'og_gene_membership.tsv',sep='\t',index=False)
        elif target.exists():target.unlink()
        write_json(metrics.parent/'og_status.json',qc)
