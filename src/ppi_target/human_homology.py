"""New BLASTP with exact FASTA, all available named isoforms and full HSP evidence."""
from collections import defaultdict
from pathlib import Path
import gzip,hashlib,json,re,shutil,subprocess
from io import StringIO
import pandas as pd
from Bio import SeqIO
from .identifiers import read_gene_info
from .sources import fetch
from .provenance import sha256,write_json,now


def union_length(intervals):
    intervals=sorted((min(a,b),max(a,b)) for a,b in intervals)
    total=0;end=0
    for a,b in intervals:
        total+=max(0,b-max(a-1,end));end=max(end,b)
    return total


def execute(root,cfg,run,hubs,out):
    raw=out/'raw';raw.mkdir();queries=out/'queries';queries.mkdir();db=out/'database';db.mkdir()
    mapping=read_gene_info(root/cfg['biogrid']['gene_info']['path'],6239)
    reverse=defaultdict(set)
    for gid,wb in mapping.items():
        if len(wb)==1:reverse[next(iter(wb))].add(gid)
    sequence_map=[];fasta={};failures=[];retrievals=[];status={};entries={};release=None
    for h in hubs.to_dict('records'):
        gene=h['WBGene'];gids=reverse.get(gene,set());count_before=len(sequence_map)
        try:
            if len(gids)!=1:raise ValueError('GeneID mapping is not unique')
            gid=next(iter(gids));params={'query':f'(organism_id:6239) AND (xref:geneid-{gid})','format':'json','size':500}
            response,meta=fetch('https://rest.uniprot.org/uniprotkb/search',raw/(gene+'_search'),[gene,gid],params)
            retrievals.append(meta)
            payload=response.json();results=payload.get('results')
            if results is None:raise ValueError('UniProt response has no results array')
            if response.links.get('next'):raise ValueError('UniProt search requires additional pagination; do not truncate')
            current=meta['database_version']
            if release and release!=current:raise ValueError('UniProt release changed during queries')
            release=current
            for entry in results:
                acc=entry['primaryAccession'];taxon=entry['organism']['taxonId'];refs=entry.get('uniProtKBCrossReferences',[])
                generef={p['value'] for r in refs if r['database']=='WormBase' for p in r.get('properties',[]) if p['key']=='GeneId'}
                gidrefs={r['id'] for r in refs if r['database']=='GeneID'}
                if taxon!=6239 or (gene not in generef and gid not in gidrefs):raise ValueError('Sequence reference does not establish hub identity')
                entries[acc]=entry
                sequences={acc:entry['sequence']['value']}
                isoforms={i for c in entry.get('comments',[]) if c.get('commentType')=='ALTERNATIVE PRODUCTS' for iso in c.get('isoforms',[]) for i in iso.get('isoformIds',[])}
                for iso in sorted(isoforms):
                    r,im=fetch(f'https://rest.uniprot.org/uniprotkb/{iso}.fasta',raw/(gene+'_'+iso),[gene,iso]);retrievals.append(im)
                    if im['database_version']!=release:raise ValueError('Isoform release mismatch')
                    seqs=list(SeqIO.parse(StringIO(r.text),'fasta'))
                    if len(seqs)!=1:raise ValueError('Invalid isoform FASTA')
                    sequences[iso]=str(seqs[0].seq)
                for sid,seq in sequences.items():
                    if sid in fasta and fasta[sid]!=seq:raise ValueError('Conflicting sequence for accession')
                    fasta[sid]=seq
                    sequence_map.append({'WBGene':gene,'UniProt':acc,'Orthogroup':h.get('Orthogroup',''),'query_id':sid,
                         'sequence_version':entry['entryAudit']['sequenceVersion'],'entry_version':entry['entryAudit']['entryVersion'],
                         'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'query_length':len(seq),'source_version':release})
            status[gene]='queries_available' if len(sequence_map)>count_before else 'no_current_sequence_returned'
        except Exception as exc:
            failures.append({'WBGene':gene,'error_type':type(exc).__name__,'message':str(exc)})
            status[gene]='query_retrieval_failed_or_incomplete'
        print('Sequence retrieval',gene,status[gene],flush=True)
    write_json(queries/'entries.json',entries);write_json(queries/'retrieval_errors.json',failures)
    qmap=pd.DataFrame(sequence_map)
    if qmap.empty:raise ValueError('No validated query sequences available; no BLAST scientific result generated')
    qmap=qmap.drop_duplicates();qmap.to_csv(queries/'sequence_mapping.tsv',sep='\t',index=False)
    queryfile=queries/'hub_queries.fasta'
    queryfile.write_text(''.join(f'>{q}\n{fasta[q]}\n' for q in sorted(fasta)))
    target_source=root/cfg['blast']['target_gzip'];target=db/'human_target.fasta'
    with gzip.open(target_source,'rb') as fin,target.open('wb') as fout:shutil.copyfileobj(fin,fout)
    shutil.copyfile(root/cfg['blast']['target_headers'],db/'target.headers')
    headers=(db/'target.headers').read_text()
    m=re.search(r'x-uniprot-release:\s*([^\r\n]+)',headers,re.I);targetrelease=m.group(1) if m else None
    if targetrelease!=release:raise ValueError(f'Query release {release} differs from target release {targetrelease}')
    names={};targetcount=0
    for seq in SeqIO.parse(target,'fasta'):
        if 'OX=9606' not in seq.description:raise ValueError('Nonhuman or taxon-unverified target sequence')
        if seq.id in names:raise ValueError('Duplicate target FASTA identifier')
        names[seq.id]=seq.description;targetcount+=1
    if not targetcount:raise ValueError('Empty human target FASTA')
    for exe in ('blastp','makeblastdb'):
        executable=shutil.which(exe)
        if not executable:raise FileNotFoundError(exe+' is not installed')
    versions={exe:subprocess.check_output([exe,'-version'],text=True) for exe in ('blastp','makeblastdb')}
    for exe in versions:(db/(exe+'_help.txt')).write_text(subprocess.check_output([exe,'-help'],text=True))
    build=['makeblastdb','-in',str(target),'-dbtype','prot','-parse_seqids','-out',str(db/'human')]
    options=cfg['blast'];fields='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp stitle'
    command=['blastp','-task','blastp','-query',str(queryfile),'-db',str(db/'human'),'-evalue',str(options['evalue']),'-matrix',options['matrix'],
        '-gapopen',str(options['gapopen']),'-gapextend',str(options['gapextend']),'-word_size',str(options['word_size']),'-comp_based_stats',str(options['comp_based_stats']),
        '-seg',options['seg'],'-soft_masking',options['soft_masking'],'-num_threads',str(options['threads']),'-max_target_seqs',str(targetcount),'-outfmt','6 '+fields,'-out',str(raw/'blast.tsv')]
    bm={'source':'UniProt + local BLASTP','database_version':release,'target_source':options['target_source'],'retrieved_at':now(),'timezone':'UTC',
        'target_retrieved_at':re.search(r'^date:\s*(.*)$',headers,re.M|re.I).group(1).strip() if re.search(r'^date:',headers,re.M|re.I) else None,
        'target_sequence_count':targetcount,'query_sequence_count':len(fasta),'query_fasta':str(queryfile),'query_sha256':sha256(queryfile),
        'target_fasta':str(target),'target_sha256':sha256(target),'query_retrievals':retrievals,'software':versions,'commands':[build,command],
        'isoform_policy':options['isoforms'],'parameters':options,'warnings':failures,'status':'running'}
    write_json(out/'source_manifest.json',bm)
    with (db/'build.log').open('w') as log:subprocess.run(build,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=600)
    print('BLASTP:',len(fasta),'queries against',targetcount,'human sequences',flush=True)
    with (raw/'blast.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=3600)
    columns=fields.split();rawfile=raw/'blast.tsv'
    hsps=pd.read_csv(rawfile,sep='\t',names=columns) if rawfile.stat().st_size else pd.DataFrame(columns=columns)
    hits=[]
    for (query,subject),d in hsps.groupby(['qseqid','sseqid']):
        best=d.sort_values(['evalue','bitscore'],ascending=[True,False]).iloc[0]
        desc=str(best.stitle);gid=re.search(r'GN=([^ ]+)',desc);acc=subject.split('|')[1] if '|' in subject else subject
        hits.append({'query_id':query,'human_hit_id':subject,'human_gene':gid.group(1) if gid else '', 'human_UniProt':acc,'human_protein':desc,
                     'evalue':float(best.evalue),'bit_score':float(best.bitscore),'percent_identity':float(best.pident),'alignment_length':int(best.length),
                     'query_coverage':100*union_length(zip(d.qstart,d.qend))/int(best.qlen),'subject_coverage':100*union_length(zip(d.sstart,d.send))/int(best.slen),
                     'query_length':int(best.qlen),'subject_length':int(best.slen),'HSP_count':len(d),'identity_scope':'best_HSP','interpretation':'sequence_similarity_candidate'})
    h=pd.DataFrame(hits,columns=['query_id','human_hit_id','human_gene','human_UniProt','human_protein','evalue','bit_score','percent_identity','alignment_length','query_coverage','subject_coverage','query_length','subject_length','HSP_count','identity_scope','interpretation'])
    full=qmap.merge(h,on='query_id',how='inner',validate='many_to_many',suffixes=('','_alignment'))
    full=full.sort_values(['WBGene','evalue','bit_score','query_coverage','subject_coverage','human_hit_id'],ascending=[True,True,False,False,False,True])
    full.to_csv(out/'human_similarity.tsv',sep='\t',index=False)
    summaries=[]
    for gene in hubs.WBGene:
        d=full[full.WBGene==gene];q=qmap[qmap.WBGene==gene]
        row={'WBGene':gene,'blast_query_status':status[gene],'blast_query_sequences':q.query_id.nunique(),'blast_hits':len(d),
             'blast_human_proteins':d.human_hit_id.nunique(),'blast_interpretation':'sequence similarity; not orthology',
             'blast_retrieval':str(out.relative_to(root))}
        if len(d):
            best=d.iloc[0];row.update({f'best_{col}':best[col] for col in ['human_UniProt','human_gene','evalue','bit_score','percent_identity','alignment_length','query_coverage','subject_coverage']})
        summaries.append(row)
    pd.DataFrame(summaries).to_csv(out/'candidate_summary.tsv',sep='\t',index=False)
    state='PARTIAL' if failures or any(v!='queries_available' for v in status.values()) else 'SUCCESS'
    bm.update(status=state,finished_at=now(),hit_rows=len(full),raw_sha256=sha256(rawfile));write_json(out/'source_manifest.json',bm)
    return {'status':state,'records':len(full),'reason':f'{len(failures)} query retrieval failures; missing sequences never mean no similarity' if state=='PARTIAL' else '',
            'database_version':release,'retrieved_at':bm['retrieved_at']}
