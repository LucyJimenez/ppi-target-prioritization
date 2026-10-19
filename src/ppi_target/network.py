"""One physical topology/LCC, primary analysis and paired support sensitivity."""
from collections import Counter, defaultdict
import csv
import json
import math
import networkx as nx
import numpy as np
import pandas as pd
import community as cl
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score
from .provenance import write_json


def table(rows,path,columns=None):
    df=pd.DataFrame(rows,columns=columns)
    df.to_csv(path,sep='\t',index=False)
    return df


def execute(root,config,run,manifest):
    out=run/'02_network';out.mkdir()
    with (run/'01_interactions/interactions.tsv').open() as source:
        records=list(csv.DictReader(source,delimiter='\t'))
    inp=pd.read_csv(run/'01_interactions/input_mapping.tsv',sep='\t',dtype=str)
    inputs=set(inp.WBGene);meta=defaultdict(lambda:{'accessions':set(),'symbols':set()});pairs=defaultdict(list);selfs=[]
    for row in inp.itertuples():meta[row.WBGene]['accessions'].add(row.UniProt)
    for r in records:
        a,b=r['WBGene_A'],r['WBGene_B']
        assert r['Experimental System Type']=='physical' and r['Organism ID Interactor A']==r['Organism ID Interactor B']=='6239'
        assert a and b and (a in inputs or b in inputs)
        for side,n in [('A',a),('B',b)]:
            symbol=r.get(f'Official Symbol Interactor {side}','')
            if symbol and symbol!='-':meta[n]['symbols'].add(symbol)
            for col in ('SWISS-PROT','TREMBL'):
                for acc in r.get(f'{col} Accessions Interactor {side}','').split('|'):
                    if acc and acc!='-':meta[n]['accessions'].add(acc)
        if a==b:selfs.append(r)
        else:pairs[tuple(sorted((a,b)))].append(r)
    table(selfs,out/'self_evidence.tsv',columns=list(records[0]) if records else ['WBGene_A','WBGene_B'])
    g=nx.Graph();g.add_nodes_from(sorted(inputs|set(meta)));edgerows=[]
    for (a,b),rs in sorted(pairs.items()):
        ids=sorted({r['#BioGRID Interaction ID'] for r in rs});pubs=sorted({r['Publication Source'] for r in rs if r['Publication Source'] not in ('','-')});systems=sorted({r['Experimental System'] for r in rs})
        if not pubs:raise ValueError(f'No publication for {a},{b}; weighted analysis blocked')
        d={'evidence_count':len(ids),'publication_count':len(pubs),'experimental_system_count':len(systems),
           'experimental_systems':';'.join(systems),'publications':';'.join(pubs),'interaction_ids':';'.join(ids),
           'support':len(pubs),'distance':1/len(pubs),'raw_rows':len(rs),'duplicate_extra_rows':len(rs)-len(ids)}
        edgerows.append({'WBGene_A':a,'WBGene_B':b,**d});g.add_edge(a,b,**d)
    components=sorted(nx.connected_components(g),key=lambda s:(-len(s),min(s)))
    if not components or len(components[0])<2:raise ValueError('No nontrivial LCC: network analysis not applicable')
    lcc=g.subgraph(sorted(components[0])).copy();nodes=sorted(lcc)
    assert nx.is_connected(lcc) and not list(nx.selfloop_edges(g))
    assert sum(dict(g.degree()).values())==2*g.number_of_edges()
    node_rows=[{'WBGene':n,'UniProt':';'.join(sorted(meta[n]['accessions'])),'Gene_symbol':';'.join(sorted(meta[n]['symbols'])),
               'Orthogroup':'','original_input':n in inputs,'in_LCC':n in lcc,
               'component':next(i+1 for i,s in enumerate(components) if n in s)} for n in sorted(g)]
    nt=table(node_rows,out/'nodes.tsv');table(edgerows,out/'edges.tsv')
    nt[nt.in_LCC].to_csv(out/'lcc_nodes.tsv',sep='\t',index=False)
    table([r for r in edgerows if r['WBGene_A'] in lcc and r['WBGene_B'] in lcc],out/'lcc_edges.tsv')
    nx.write_graphml(g,out/'network.graphml');nx.write_graphml(lcc,out/'lcc.graphml')
    table([{'step':'resolved_evidence','records':len(records),'nodes':len(g),'edges':''},
           {'step':'self_records_removed','records':len(records)-len(selfs),'nodes':len(g),'edges':''},
           {'step':'unique_pairs','records':len(records)-len(selfs),'nodes':len(g),'edges':g.number_of_edges()},
           {'step':'LCC','records':'','nodes':len(lcc),'edges':lcc.number_of_edges()}],out/'cleaning_counts.tsv')
    table([{'component':i+1,'nodes':len(c),'edges':g.subgraph(c).number_of_edges()} for i,c in enumerate(components)],out/'connected_components.tsv')
    table([{'WBGene_A':min(r['WBGene_A'],r['WBGene_B']),'WBGene_B':max(r['WBGene_A'],r['WBGene_B']),**r} for r in records if r['WBGene_A']!=r['WBGene_B']],out/'edge_evidence.tsv')
    summary={'full_nodes':len(g),'full_edges':g.number_of_edges(),'lcc_nodes':len(lcc),'lcc_edges':lcc.number_of_edges(),'components':len(components),'isolates':len(list(nx.isolates(g)))}
    write_json(out/'network_summary.json',summary)
    manifest['stages']['02']={'status':'SUCCESS','summary':summary}
    return summary
