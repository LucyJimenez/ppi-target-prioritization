"""Stage 1: auditable physical record selection, without graph construction."""
import csv
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from .identifiers import read_mapping, mapping_qc, read_gene_info, resolve
from .provenance import sha256, write_json, write_tsv

REQUIRED = {'#BioGRID Interaction ID', 'Entrez Gene Interactor A', 'Entrez Gene Interactor B',
            'Experimental System Type', 'Experimental System', 'Publication Source',
            'Organism ID Interactor A', 'Organism ID Interactor B'}
EXTRA = ['WBGene_A', 'WBGene_B', 'resolution_A', 'resolution_B', 'input_A', 'input_B', 'qc_status']


def inspect_records(records, mapping, inputs, taxon):
    source_types, organisms, systems = Counter(), Counter(), Counter()
    observed, physical_observed = set(), set()
    candidates, ready, rejected, genetic, cross = [], [], [], [], []
    seen_ids = {}
    for source in records:
        r = dict(source)
        iid = r['#BioGRID Interaction ID']
        if not iid or iid == '-':
            raise ValueError('BioGRID record has no interaction ID')
        if iid in seen_ids and source != seen_ids[iid]:
            raise ValueError(f'Conflicting records for BioGRID ID {iid}')
        seen_ids[iid] = source
        typ = r['Experimental System Type'].strip().lower()
        source_types[typ] += 1
        taxa = (r['Organism ID Interactor A'], r['Organism ID Interactor B'])
        resolved = []
        for side, organism in zip(('A','B'), taxa):
            g, status = resolve(r[f'Entrez Gene Interactor {side}'], mapping) if organism == str(taxon) else ('', 'non_target_taxon')
            r[f'WBGene_{side}'], r[f'resolution_{side}'] = g, status
            r[f'input_{side}'] = g in inputs
            resolved.append(g)
        hits = set(resolved) & inputs
        observed.update(hits)
        if typ == 'physical':
            organisms[taxa] += 1
            systems[r['Experimental System']] += 1
            physical_observed.update(hits)
            if taxa != (str(taxon), str(taxon)):
                c = dict(r, qc_status='non_target_organism_pair')
                cross.append(c)
        if not hits:
            continue
        if typ != 'physical':
            genetic.append(dict(r, qc_status='non_physical'))
            continue
        if taxa != (str(taxon), str(taxon)):
            status = 'non_target_organism_pair'
        elif any(r[f'resolution_{s}'] != 'resolved' for s in ('A','B')):
            status = 'unresolved_or_ambiguous_partner'
        else:
            status = 'ready'
        r['qc_status'] = status
        candidates.append(r)
        (ready if status == 'ready' else rejected).append(r)
    return {'source_types': source_types, 'organisms': organisms, 'systems': systems,
            'observed': observed, 'physical_observed': physical_observed, 'candidates': candidates,
            'ready': ready, 'rejected': rejected, 'non_physical': genetic, 'cross': cross}


def execute(root, config, run_dir, manifest):
    c = config['biogrid']
    if c['interaction_type'] != 'physical' or c['taxon'] != 6239:
        raise ValueError('This approved pipeline requires physical interactions and taxon 6239')
    out = run_dir / '01_interactions'
    out.mkdir()
    input_path = (root / config['input']['path']).resolve()
    manifest['input_files'].append({'path': str(input_path), 'sha256': sha256(input_path)})
    from .annotation import load
    rows,canonical,ogqc=load(root,config,out)
    canonical.to_csv(out/'candidate_annotations.tsv',sep='\t',index=False)
    manifest['orthogroup_annotation']=ogqc
    qc=mapping_qc(rows)
    qc['orthogroup_status']=ogqc['orthogroup_annotation']
    write_json(out/'mapping_qc.json',qc)
    write_tsv(out/'input_mapping.tsv',rows,['source_row','WBGene','UniProt','Orthogroup'])
    inputs = {r['WBGene'] for r in rows}
    snapshot = {}
    for name in ('tab3', 'gene_info'):
        spec = c[name]
        src = (root / spec['path']).resolve()
        digest = sha256(src)
        if digest != spec['sha256']:
            raise ValueError(f'{name} snapshot checksum differs from pinned config: {src}')
        target = out / 'raw' / src.name
        target.parent.mkdir(exist_ok=True)
        shutil.copyfile(src, target)
        metadata = root / spec['retrieval_metadata']
        shutil.copyfile(metadata, target.parent / (name + '.retrieval_metadata.json'))
        snapshot[name] = target
        manifest['database_snapshots'].append({**spec, 'name': name, 'raw_response': str(target.relative_to(root)),
                                                'retrieval_mode': 'reuse_explicit_pinned_snapshot',
                                                'sha256': digest})
    mapping = read_gene_info(snapshot['gene_info'], c['taxon'])
    with snapshot['tab3'].open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        if not REQUIRED <= set(fields or []):
            raise ValueError(f'BioGRID TAB3 missing columns: {REQUIRED - set(fields or [])}')
        result = inspect_records(reader, mapping, inputs, c['taxon'])
    for key, filename in [('candidates','candidate_physical_records.tsv'), ('ready','interactions.tsv'),
                          ('rejected','unresolved_identifiers.tsv'), ('non_physical','input_non_physical_records.tsv'),
                          ('cross','cross_species_interactions.tsv')]:
        write_tsv(out/filename, result[key], fields+EXTRA)
    reverse = defaultdict(set)
    for entrez, genes in mapping.items():
        for gene in genes:
            reverse[gene].add(entrez)
    status_rows = []
    for gene in sorted(inputs):
        ids = sorted(reverse.get(gene, set()))
        status_rows.append({'WBGene':gene, 'NCBI_GeneIDs':';'.join(ids),
                            'mapping_status':'resolved' if len(ids)==1 else ('unresolved' if not ids else 'ambiguous'),
                            'observed_any_interaction':gene in result['observed'],
                            'physical_record_status':'physical_record_found' if gene in result['physical_observed'] else 'no_matching_physical_record'})
    write_tsv(out/'gene_query_status.tsv', status_rows, list(status_rows[0]))
    write_tsv(out/'interaction_type_counts.tsv', [{'type':k,'records':v} for k,v in sorted(result['source_types'].items())], ['type','records'])
    write_tsv(out/'physical_organism_pairs.tsv', [{'taxon_A':k[0],'taxon_B':k[1],'records':v} for k,v in sorted(result['organisms'].items())], ['taxon_A','taxon_B','records'])
    systems = []
    for scope, counter in [('source_physical', result['systems']), ('candidate_physical', Counter(r['Experimental System'] for r in result['candidates']))]:
        systems.extend({'scope':scope,'experimental_system':k,'records':v} for k,v in sorted(counter.items()))
    write_tsv(out/'physical_experimental_systems.tsv', systems, ['scope','experimental_system','records'])
    summary = {**qc, 'source_records':sum(result['source_types'].values()),
               'source_types':dict(result['source_types']), 'input_with_physical':len(result['physical_observed']),
               'input_without_physical':len(inputs-result['physical_observed']),
               'candidate_physical_records':len(result['candidates']),
               'resolved_physical_records':len(result['ready']), 'excluded_physical_records':len(result['rejected']),
               'source_physical_cross_species':len(result['cross']),
               'network_constructed':False, 'lcc_calculated':False}
    write_json(out/'qc_summary.json', summary)
    shutil.copyfile(out/'candidate_physical_records.tsv',out/'interaction_evidence.tsv')
    write_json(out/'source_manifest.json', {'sources':manifest['database_snapshots']})
    manifest['stages']['01'] = {'status':'SUCCESS', 'summary':summary}
    return summary
