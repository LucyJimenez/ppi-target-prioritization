"""Shared I/O, identifier validation, downloads and provenance. No analysis entry point."""
__version__ = '0.4.0'
'WBGene identity and explicit protein associations, never accession nodes.'
import csv
import gzip
import re
from collections import Counter, defaultdict
WB = re.compile('WBGene\\d{8}\\Z')
UP = re.compile('(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\\d+)?\\Z')

def read_mapping(path, header=False):
    with open(path, newline='') as f:
        rows = list(csv.reader(f, delimiter='\t'))
    if header:
        if not rows or rows.pop(0) != ['WBGene', 'UniProt']:
            raise ValueError('Expected exact header WBGene<TAB>UniProt')
    if not rows:
        raise ValueError('Empty input mapping')
    valid, errors = ([], [])
    for i, row in enumerate(rows, start=2 if header else 1):
        if len(row) != 2 or not WB.fullmatch(row[0]) or (not UP.fullmatch(row[1])):
            errors.append({'source_row': i, 'raw': row, 'reason': 'invalid_columns_or_identifier_syntax'})
        else:
            valid.append({'source_row': i, 'WBGene': row[0], 'UniProt': row[1]})
    return (valid, errors)

def mapping_qc(rows):
    pairs = Counter(((r['WBGene'], r['UniProt']) for r in rows))
    genes, proteins = (defaultdict(set), defaultdict(set))
    for g, p in pairs:
        genes[g].add(p)
        proteins[p].add(g)
    return {'valid_rows': len(rows), 'unique_WBGene': len(genes), 'unique_UniProt': len(proteins), 'unique_pairs': len(pairs), 'duplicate_pair_extra_rows': sum((n - 1 for n in pairs.values())), 'one_gene_multiple_accessions': {k: sorted(v) for k, v in genes.items() if len(v) > 1}, 'one_accession_multiple_genes': {k: sorted(v) for k, v in proteins.items() if len(v) > 1}, 'orthogroups': None, 'orthogroup_status': 'not_provided'}

def read_gene_info(path, taxon):
    """A GeneID is resolvable only if all matching current xrefs agree."""
    mapping = defaultdict(set)
    with gzip.open(path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not {'#tax_id', 'GeneID', 'dbXrefs'} <= set(reader.fieldnames or []):
            raise ValueError('NCBI gene_info schema missing required columns')
        for row in reader:
            if row['#tax_id'] == str(taxon):
                mapping[row['GeneID']].update(re.findall('WBGene\\d{8}(?!\\d)', row['dbXrefs']))
    return dict(mapping)

def resolve(gene_id, mapping):
    values = mapping.get(gene_id, set())
    if len(values) == 1:
        return (next(iter(values)), 'resolved')
    return ('', 'ambiguous' if values else 'unresolved')
'Immutable run bookkeeping; no scientific conclusions inferred from failures.'
from datetime import datetime, timezone
from pathlib import Path
import csv
import math
import hashlib
import json
import platform
import subprocess
import sys
import uuid
import yaml

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def now():
    return datetime.now(timezone.utc).isoformat()

def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)

    def clean(x):
        if isinstance(x, float) and (not math.isfinite(x)):
            return None
        if isinstance(x, dict):
            return {k: clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        return x
    Path(path).write_text(json.dumps(clean(value), indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='raise')
        w.writeheader()
        w.writerows(rows)

def load_config(path):
    path = Path(path).resolve()
    try:
        config = json.loads(path.read_text())
    except json.JSONDecodeError:
        config = yaml.safe_load(path.read_text())
    root = (path.parent / config['project_root']).resolve()
    return (path, root, config)

def new_run(root, config_path, config, run_id=None):
    run_id = run_id or datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8]
    if not run_id or any((c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in run_id)):
        raise ValueError('run_id must contain only letters, digits, underscore or hyphen')
    target = root / 'results' / 'runs' / run_id
    target.mkdir(parents=True, exist_ok=False)
    git = subprocess.run(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True, capture_output=True)
    manifest = {'run_id': run_id, 'started_at': now(), 'finished_at': None, 'timezone': 'UTC', 'git_commit': git.stdout.strip() if git.returncode == 0 else None, 'python_version': platform.python_version(), 'package_versions': {'ppi-target-prioritization': __version__}, 'configuration_file': str(config_path), 'configuration_sha256': sha256(config_path), 'resolved_configuration': config, 'input_files': [], 'database_snapshots': [], 'orthogroup_annotation': {'status': 'not_provided'}, 'stages': {}, 'outputs': [], 'scripts_executed': [sys.argv[0]], 'commands': [sys.argv], 'random_seeds': {}, 'warnings': [] if git.returncode == 0 else ['No Git commit available; workspace is not a Git repository'], 'errors': [], 'final_status': 'running'}
    write_json(target / 'run_manifest.json', manifest)
    return (target, manifest)
'Optional OG annotation; never part of network identity or topology cache keys.'
import csv, hashlib, json, re
from collections import defaultdict, Counter
from pathlib import Path
import pandas as pd

def load(root, config, out=None):
    path = (root / config['input']['path']).resolve()
    with path.open(newline='') as f:
        data = list(csv.reader(f, delimiter='\t'))
    if config['input'].get('header', True):
        header = data.pop(0) if data else []
        if not {'WBGene', 'UniProt'} <= set(header):
            raise ValueError('Input requires WBGene and UniProt columns')
    else:
        header = ['WBGene', 'UniProt']
    rows = []
    for i, r in enumerate(data, 2 if config['input'].get('header', True) else 1):
        if len(r) != len(header):
            raise ValueError(f'Input row {i}: incorrect column count')
        x = dict(zip(header, r))
        g = x['WBGene']
        u = x['UniProt']
        og = x.get('Orthogroup', '')
        if not WB.fullmatch(g) or not UP.fullmatch(u):
            raise ValueError(f'Input row {i}: malformed required identifier')
        if og and (not re.fullmatch('OG\\d+', og)):
            raise ValueError(f'Input row {i}: malformed OG')
        rows.append({'source_row': i, 'WBGene': g, 'UniProt': u, 'Orthogroup': og})
    if not rows:
        raise ValueError('Input is empty')
    genes = {x['WBGene'] for x in rows}
    ogs = defaultdict(set)
    for r in rows:
        if r['Orthogroup']:
            ogs[r['WBGene']].add(r['Orthogroup'])
    ogpath = root / config.get('annotation', {}).get('orthogroup_file', 'data/raw/orthogroup_wbgene_mapping.tsv')
    issues = []
    loaded = 0
    seen = set()
    fatal = []
    if ogpath.exists():
        with ogpath.open(newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            if set(reader.fieldnames or []) != {'Orthogroup', 'WBGene'}:
                fatal.append('OG file requires exactly Orthogroup and WBGene columns')
            else:
                for i, r in enumerate(reader, 2):
                    g = r['WBGene']
                    o = r['Orthogroup']
                    reason = None
                    if not WB.fullmatch(g or '') or not re.fullmatch('OG\\d+', o or ''):
                        reason = 'malformed_or_missing_identifier'
                        fatal.append(f'OG row {i}: {reason}')
                    elif (g, o) in seen:
                        reason = 'duplicate_row'
                    elif g not in genes:
                        reason = 'WBGene_not_in_canonical_mapping'
                    else:
                        ogs[g].add(o)
                        loaded += 1
                    seen.add((g, o))
                    if reason:
                        issues.append({'row': i, 'WBGene': g, 'Orthogroup': o, 'reason': reason})
    for g, values in ogs.items():
        if len(values) > 1:
            issues.append({'WBGene': g, 'Orthogroup': ';'.join(sorted(values)), 'reason': 'multiple_OG_assignments_preserved'})
    for r in rows:
        r['Orthogroup'] = ';'.join(sorted(ogs[r['WBGene']]))
    canonical = pd.DataFrame(rows)[['WBGene', 'UniProt', 'Orthogroup']].drop_duplicates()
    available = ogpath.exists() or any(ogs.values())
    qc = {'orthogroup_annotation': 'available' if available else 'not_provided', 'file_detected': ogpath.exists(), 'file': str(ogpath) if ogpath.exists() else None, 'sha256': sha256(ogpath) if ogpath.exists() else None, 'mappings_loaded': loaded, 'unique_OGs': len(set().union(*ogs.values())) if ogs else 0, 'WBGene_with_OG': sum((bool(ogs[g]) for g in genes)), 'WBGene_without_OG': sum((not ogs[g] for g in genes)), 'issues': issues, 'errors': fatal, 'OG_aware_analyses': 'YES' if available and any(ogs.values()) else 'NO', 'OG_analysis_status': 'SUCCESS' if any(ogs.values()) else 'NOT_RUN', 'reason': '' if any(ogs.values()) else 'orthogroup annotation not provided'}
    identity = sorted({(r['WBGene'], r['UniProt']) for r in rows})
    qc['identity_checksum'] = hashlib.sha256(json.dumps(identity).encode()).hexdigest()
    if out:
        out.mkdir(parents=True, exist_ok=True)
        write_json(out / 'orthogroup_qc.json', qc)
        pd.DataFrame(issues, columns=['row', 'WBGene', 'Orthogroup', 'reason']).to_csv(out / 'orthogroup_issues.tsv', sep='\t', index=False)
    if fatal:
        raise ValueError('; '.join(fatal))
    return (rows, canonical, qc)

def propagate(run, canonical, qc):
    proteins = canonical.groupby('WBGene').UniProt.apply(lambda x: ';'.join(sorted(set(x)))).to_dict()
    node_file = run / '02_network/nodes.tsv'
    if node_file.exists():
        node_table = pd.read_csv(node_file, sep='\t', keep_default_na=False)
        proteins = {**dict(zip(node_table.WBGene, node_table.UniProt)), **proteins}
    og = canonical.groupby('WBGene').Orthogroup.apply(lambda x: ';'.join(sorted({o for v in x for o in str(v).split(';') if o}))).to_dict()
    for stage in ['01_interactions', '02_network', '03_network_analysis', '04_functional_context', '05_biological_evidence', '06_target_profiles']:
        for f in (run / stage).rglob('*.tsv') if (run / stage).exists() else []:
            if 'raw' in f.parts or 'database' in f.parts:
                continue
            frame = pd.read_csv(f, sep='\t', keep_default_na=False)
            if 'WBGene' in frame.columns:
                frame['Orthogroup'] = frame.WBGene.map(og).fillna('')
                if 'UniProt' not in frame:
                    frame['UniProt'] = frame.WBGene.map(proteins).fillna('')
                frame.to_csv(f, sep='\t', index=False)
    metrics = run / '03_network_analysis/node_metrics.tsv'
    if metrics.exists():
        df = pd.read_csv(metrics, sep='\t', keep_default_na=False)
        detail = []
        for r in df.to_dict('records'):
            for o in str(r.get('Orthogroup', '')).split(';'):
                if o:
                    detail.append({'Orthogroup': o, 'WBGene': r['WBGene'], 'community': r['community'], 'degree': r['degree'], 'robust_hub': r['robust_hub']})
        target = metrics.parent / 'og_summary.tsv'
        if detail:
            d = pd.DataFrame(detail)
            d.groupby('Orthogroup').agg(WBGene_count=('WBGene', 'nunique'), communities=('community', 'nunique'), robust_hubs=('robust_hub', 'sum')).reset_index().to_csv(target, sep='\t', index=False)
            d.to_csv(metrics.parent / 'og_gene_membership.tsv', sep='\t', index=False)
        elif target.exists():
            target.unlink()
        write_json(metrics.parent / 'og_status.json', qc)
'Bounded HTTP requests, immutable response evidence and independent source statuses.'
import json, time, traceback
from pathlib import Path
import requests

class RetrievalError(RuntimeError):

    def __init__(self, message, attempts=0, http_status=None):
        super().__init__(message)
        self.attempts = attempts
        self.http_status = http_status

def fetch(url, folder, identifiers=None, params=None):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    start = now()
    last = None
    status = None
    for attempt in range(1, 4):
        try:
            r = requests.get(url, params=params, timeout=(30, 120))
            status = r.status_code
            raw = folder / f'response_{attempt}.bin'
            raw.write_bytes(r.content)
            headers = {k: v for k, v in r.headers.items() if k.lower() not in ('set-cookie', 'cookie', 'authorization')}
            write_json(folder / f'headers_{attempt}.json', headers)
            if r.status_code >= 400:
                last = f'HTTP {r.status_code}: {r.text[:200]}'
                if r.status_code not in (408, 429, 500, 502, 503, 504):
                    break
                delay = min(30, float(r.headers.get('Retry-After', 2 ** (attempt - 1))))
            else:
                meta = {'source': url.split('/')[2], 'database_version': r.headers.get('x-uniprot-release'), 'release_date': r.headers.get('x-uniprot-release-date'), 'retrieved_at': now(), 'timezone': 'UTC', 'endpoint': r.url, 'query_method': 'GET', 'identifiers': identifiers or [], 'parameters': params or {}, 'status': 'SUCCESS', 'http_status': status, 'attempts': attempt, 'raw_response': str(raw.resolve()), 'sha256': sha256(raw), 'software': {'requests': requests.__version__}, 'warnings': [], 'errors': []}
                write_json(folder / 'source_manifest.json', meta)
                return (r, meta)
        except requests.RequestException as e:
            last = str(e)
            delay = 2 ** (attempt - 1)
        if attempt < 3:
            time.sleep(delay)
    error = {'status': 'FAILED', 'error_type': 'HTTPError' if status else 'NetworkError', 'http_status': status, 'message': last, 'attempts': attempt, 'started_at': start, 'failed_at': now(), 'affected_identifiers': identifiers or [], 'endpoint': url}
    write_json(folder / 'error.json', error)
    write_json(folder / 'source_manifest.json', error)
    raise RetrievalError(last, attempt, status)

def independent_sources(names, callbacks, base, identifiers):
    """Every requested callback is attempted even if every earlier source failed."""
    outcomes = []
    for name in names:
        start = now()
        directory = Path(base) / name
        directory.mkdir(parents=True, exist_ok=True)
        try:
            result = callbacks[name](directory)
            status = result.get('status', 'SUCCESS')
            if status not in ('SUCCESS', 'PARTIAL', 'FAILED', 'NOT_RUN'):
                raise ValueError('Invalid evidence status')
            outcome = {'source': name, 'status': status, 'records': result.get('records', 0), 'reason': result.get('reason', ''), 'database_version': result.get('database_version', ''), 'retrieved_at': result.get('retrieved_at', ''), 'started_at': start, 'finished_at': now(), 'error_type': '', 'http_status': '', 'attempts': result.get('attempts', 0)}
        except Exception as e:
            error = {'source': name, 'status': 'FAILED', 'error_type': type(e).__name__, 'http_status': getattr(e, 'http_status', None), 'message': str(e), 'attempts': getattr(e, 'attempts', 0), 'started_at': start, 'failed_at': now(), 'affected_identifiers': list(identifiers), 'traceback': traceback.format_exc()}
            write_json(directory / 'error.json', error)
            write_json(directory / 'source_manifest.json', error)
            outcome = {'source': name, 'status': 'FAILED', 'records': 0, 'reason': str(e), 'database_version': '', 'retrieved_at': '', 'started_at': start, 'finished_at': now(), 'error_type': type(e).__name__, 'http_status': getattr(e, 'http_status', ''), 'attempts': getattr(e, 'attempts', 0)}
        outcomes.append(outcome)
        print(f'[Stage 05][{name}] {outcome['status']}: {outcome['reason'] or str(outcome['records']) + ' records'}', flush=True)
        if outcome['status'] == 'FAILED':
            print('No usable evidence from this source; continuing independent sources.', flush=True)
    return outcomes
'Explicit, versioned source refresh; successful snapshots are never overwritten.'
from copy import deepcopy
from pathlib import Path
import gzip, json, re, zipfile

def download(root, url, directory, name, kind='binary'):
    response, meta = fetch(url, directory / 'request')
    data = response.content
    dest = directory / name
    if kind == 'gzip' and (not data.startswith(b'\x1f\x8b')):
        if data.startswith(b'>'):
            data = gzip.compress(data, mtime=0)
        else:
            raise ValueError('Expected gzip biological resource, got other content')
    dest.write_bytes(data)
    if kind == 'gzip':
        with gzip.open(dest, 'rb') as f:
            f.read(100)
    meta.update(raw_response=str(dest.relative_to(root)), sha256=sha256(dest))
    manifest = directory / 'source_manifest.json'
    write_json(manifest, meta)
    return ({'path': str(dest.relative_to(root)), 'sha256': sha256(dest), 'metadata': str(manifest.relative_to(root)), 'retrieval_metadata': str(manifest.relative_to(root)), **{k: v for k, v in meta.items() if k not in ('raw_response', 'sha256')}}, response)

def refresh(root, cfg, stages, source, run):
    cfg = deepcopy(cfg)
    base = run / 'refreshed_sources'
    base.mkdir(exist_ok=True)
    if '01' in stages:
        directory = base / 'biogrid'
        directory.mkdir()
        spec, response = download(root, 'https://downloads.thebiogrid.org/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip', directory, 'organism.zip')
        with zipfile.ZipFile(root / spec['path']) as archive:
            if archive.testzip():
                raise ValueError('BioGRID archive CRC failure')
            files = [n for n in archive.namelist() if re.fullmatch('BIOGRID-ORGANISM-Caenorhabditis_elegans-[\\d.]+\\.tab3\\.txt', Path(n).name)]
            if len(files) != 1:
                raise ValueError('Cannot identify exactly one C. elegans TAB3 member')
            target = directory / Path(files[0]).name
            target.write_bytes(archive.read(files[0]))
        version = re.search('elegans-([\\d.]+)\\.tab3', target.name).group(1)
        spec.update(path=str(target.relative_to(root)), sha256=sha256(target), database_version=version)
        cfg['biogrid']['tab3'] = spec
        d = base / 'ncbi'
        d.mkdir()
        cfg['biogrid']['gene_info'], _ = download(root, 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Invertebrates/Caenorhabditis_elegans.gene_info.gz', d, 'gene_info.gz', 'gzip')
    if '04' in stages:
        for name, url, fname, kind in [('gaf', 'https://current.geneontology.org/annotations/gaf/CAEEL-mod.gaf.gz', 'annotation.gaf.gz', 'gzip'), ('obo', 'https://current.geneontology.org/ontology/go-basic.obo', 'go-basic.obo', 'binary')]:
            d = base / name
            d.mkdir()
            cfg['go'][name], _ = download(root, url, d, fname, kind)
    cfg['_refresh_sources'] = ([source] if source != 'all' else ['wormbase', 'bgee', 'blast']) if '05' in stages else []
    return cfg

def refresh_evidence(root, cfg, run, source):
    base = run / 'refreshed_sources' / source
    base.mkdir(parents=True, exist_ok=True)
    if source == 'bgee':
        cfg['bgee']['source'], _ = download(root, 'https://www.bgee.org/ftp/pub/current/download/calls/expr_calls/Caenorhabditis_elegans_expr_advanced_all_conditions.tsv.gz', base, 'expression.tsv.gz', 'gzip')
        cfg['bgee']['source']['database_version'] = 'current; numeric release requires confirmation from release metadata'
    if source == 'blast':
        spec, response = download(root, cfg['blast']['target_source'], base, 'human_target.fasta.gz', 'gzip')
        headers = base / 'target.headers'
        headers.write_text('\n'.join((f'{k}: {v}' for k, v in response.headers.items() if k.lower() != 'set-cookie')))
        cfg['blast']['target_gzip'] = spec['path']
        cfg['blast']['target_headers'] = str(headers.relative_to(root))
    return cfg

def table(rows, path, columns=None):
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(path, sep='\t', index=False)
    return df

def resource(root, spec, out, name):
    source = root / spec['path']
    if sha256(source) != spec['sha256']:
        raise ValueError('Resource checksum mismatch: ' + str(source))
    target = out / source.name
    target = out / name
    shutil.copyfile(source, target)
    meta = json.loads((root / spec['metadata']).read_text())
    meta['raw_response'] = str(target.relative_to(root))
    meta['reuse_mode'] = 'explicit snapshot acquired during current project preflight'
    write_json(out / (name + '.source_manifest.json'), meta)
    return (target, meta)
