"""Explicit, versioned source refresh; successful snapshots are never overwritten."""
from copy import deepcopy
from pathlib import Path
import gzip,json,re,zipfile
from .sources import fetch
from .provenance import sha256,write_json


def download(root,url,directory,name,kind='binary'):
    response,meta=fetch(url,directory/'request')
    data=response.content;dest=directory/name
    if kind=='gzip' and not data.startswith(b'\x1f\x8b'):
        if data.startswith(b'>'):data=gzip.compress(data,mtime=0)
        else:raise ValueError('Expected gzip biological resource, got other content')
    dest.write_bytes(data)
    if kind=='gzip':
        with gzip.open(dest,'rb') as f:f.read(100)
    meta.update(raw_response=str(dest.relative_to(root)),sha256=sha256(dest))
    manifest=directory/'source_manifest.json';write_json(manifest,meta)
    return {'path':str(dest.relative_to(root)),'sha256':sha256(dest),'metadata':str(manifest.relative_to(root)),
            'retrieval_metadata':str(manifest.relative_to(root)),**{k:v for k,v in meta.items() if k not in ('raw_response','sha256')}},response


def refresh(root,cfg,stages,source,run):
    cfg=deepcopy(cfg);base=run/'refreshed_sources';base.mkdir(exist_ok=True)
    if '01' in stages:
        directory=base/'biogrid';directory.mkdir()
        spec,response=download(root,'https://downloads.thebiogrid.org/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip',directory,'organism.zip')
        with zipfile.ZipFile(root/spec['path']) as archive:
            if archive.testzip():raise ValueError('BioGRID archive CRC failure')
            files=[n for n in archive.namelist() if re.fullmatch(r'BIOGRID-ORGANISM-Caenorhabditis_elegans-[\d.]+\.tab3\.txt',Path(n).name)]
            if len(files)!=1:raise ValueError('Cannot identify exactly one C. elegans TAB3 member')
            target=directory/Path(files[0]).name;target.write_bytes(archive.read(files[0]))
        version=re.search(r'elegans-([\d.]+)\.tab3',target.name).group(1)
        spec.update(path=str(target.relative_to(root)),sha256=sha256(target),database_version=version)
        cfg['biogrid']['tab3']=spec
        d=base/'ncbi';d.mkdir()
        cfg['biogrid']['gene_info'],_=download(root,'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Invertebrates/Caenorhabditis_elegans.gene_info.gz',d,'gene_info.gz','gzip')
    if '04' in stages:
        for name,url,fname,kind in [('gaf','https://current.geneontology.org/annotations/gaf/CAEEL-mod.gaf.gz','annotation.gaf.gz','gzip'),('obo','https://current.geneontology.org/ontology/go-basic.obo','go-basic.obo','binary')]:
            d=base/name;d.mkdir();cfg['go'][name],_=download(root,url,d,fname,kind)
    cfg['_refresh_sources']=([source] if source!='all' else ['wormbase','bgee','blast']) if '05' in stages else []
    return cfg


def refresh_evidence(root,cfg,run,source):
    base=run/'refreshed_sources'/source;base.mkdir(parents=True,exist_ok=True)
    if source=='bgee':
        cfg['bgee']['source'],_=download(root,'https://www.bgee.org/ftp/pub/current/download/calls/expr_calls/Caenorhabditis_elegans_expr_advanced_all_conditions.tsv.gz',base,'expression.tsv.gz','gzip')
        # Bgee numerical release is not inferred from URL or mtime.
        cfg['bgee']['source']['database_version']='current; numeric release requires confirmation from release metadata'
    if source=='blast':
        spec,response=download(root,cfg['blast']['target_source'],base,'human_target.fasta.gz','gzip')
        headers=base/'target.headers';headers.write_text('\n'.join(f'{k}: {v}' for k,v in response.headers.items() if k.lower()!='set-cookie'))
        cfg['blast']['target_gzip']=spec['path'];cfg['blast']['target_headers']=str(headers.relative_to(root))
    return cfg
