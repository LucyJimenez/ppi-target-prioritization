"""Stage 4: go. Run --help for usage."""
from common import sha256, write_json
'Direct positive GO enrichment with a fixed within-community test family.'
from collections import defaultdict, Counter
import csv, gzip, json, shutil
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

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

def parse_obo(path):
    terms = {}
    current = None
    for line in path.read_text().splitlines():
        if line == '[Term]':
            current = {}
        elif line.startswith('['):
            current = None
        elif current is not None and ': ' in line:
            key, val = line.split(': ', 1)
            if key == 'id':
                terms[val] = current
                current['id'] = val
            elif key in ('name', 'namespace', 'is_obsolete'):
                current[key] = val
            elif key == 'alt_id':
                current.setdefault('alt_ids', []).append(val)
    aliases = {a: t for t, d in terms.items() for a in d.get('alt_ids', [])}
    return (terms, aliases)

def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind='stable')
    m = len(p)
    if m == 0:
        return p
    adj = np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    return out

def execute(root, config, run, manifest):
    out = run / '04_functional_context'
    out.mkdir()
    raw = out / 'raw'
    raw.mkdir()
    gaf, gm = resource(root, config['go']['gaf'], raw, 'annotations.gaf.gz')
    obo, om = resource(root, config['go']['obo'], raw, 'go-basic.obo')
    terms, aliases = parse_obo(obo)
    nodes = pd.read_csv(run / '03_network_analysis/node_metrics.tsv', sep='\t')
    universe = set(nodes.WBGene)
    associations = defaultdict(set)
    exclusions = Counter()
    issues = []
    ann = []
    relation_counts = Counter()
    aspects = {'biological_process': 'P', 'molecular_function': 'F', 'cellular_component': 'C'}
    with gzip.open(gaf, 'rt') as f:
        for row in csv.reader((l for l in f if not l.startswith('!')), delimiter='\t'):
            if len(row) < 17:
                raise ValueError('Invalid GAF row width')
            gene, qual, tid, aspect, taxon = (row[1], row[3], row[4], row[8], row[12])
            reason = ''
            if 'NOT' in qual.split('|'):
                reason = 'NOT'
            elif taxon.split('|')[0].split(':')[-1] != '6239':
                reason = 'non_target_taxon'
            t = aliases.get(tid, tid)
            d = terms.get(t)
            if not d:
                reason = reason or 'unknown_term'
            elif d.get('is_obsolete') == 'true':
                reason = reason or 'obsolete_term'
            elif aspects.get(d.get('namespace')) != aspect:
                reason = reason or 'aspect_mismatch'
            relation_counts[qual] += 1
            if reason:
                exclusions[reason] += 1
                issues.append({'WBGene': gene, 'GO_ID': tid, 'reason': reason, 'in_background': gene in universe})
                continue
            if gene in universe:
                associations[t].add(gene)
                ann.append({'WBGene': gene, 'GO_ID': t, 'original_GO_ID': tid, 'aspect': aspect, 'qualifier': qual, 'evidence_code': row[6], 'reference': row[5], 'assigned_by': row[14]})
    pd.DataFrame(issues, columns=['WBGene', 'GO_ID', 'reason', 'in_background']).to_csv(out / 'annotation_exclusions.tsv', sep='\t', index=False)
    problematic = [x for x in issues if x['in_background'] and x['reason'] in ('unknown_term', 'obsolete_term', 'aspect_mismatch')]
    if problematic:
        raise ValueError(f'{len(problematic)} incompatible positive GO annotations in LCC; see exclusions')
    pd.DataFrame(ann).drop_duplicates().to_csv(out / 'lcc_annotation_evidence.tsv', sep='\t', index=False)
    mapped = set().union(*associations.values()) if associations else set()
    pd.DataFrame({'WBGene': sorted(universe), 'mapped': [n in mapped for n in sorted(universe)]}).to_csv(out / 'gene_annotation_coverage.tsv', sep='\t', index=False)
    rows = []
    summary = []
    N = len(universe)
    family = sorted(associations)
    for com, group in nodes.groupby('community'):
        genes = set(group.WBGene)
        n = len(genes)
        local = []
        for t in family:
            K = len(associations[t])
            support = genes & associations[t]
            k = len(support)
            pv = float(fisher_exact([[k, n - k], [K - k, N - n - K + k]], alternative='greater').pvalue)
            local.append({'community': int(com), 'GO_ID': t, 'GO_term': terms[t]['name'], 'aspect': aspects[terms[t]['namespace']], 'observed': k, 'community_nodes': n, 'background_annotated': K, 'background_nodes': N, 'p': pv, 'fold_enrichment': k / n / (K / N), 'supporting_WBGene': ';'.join(sorted(support))})
        for r, fdr in zip(local, bh([r['p'] for r in local])):
            r['FDR'] = float(fdr)
        rows += local
        sig = [r for r in local if r['FDR'] < config['go']['fdr']]
        bp = sorted([r for r in sig if r['aspect'] == 'P'], key=lambda r: (r['FDR'], r['p'], r['GO_ID']))
        summary.append({'community': int(com), 'nodes': n, 'mapped_genes': len(genes & mapped), 'unmapped_genes': len(genes - mapped), 'testable': bool(genes & mapped), 'terms_tested': len(family) if genes & mapped else 0, 'significant_terms_all_aspects': len(sig), 'significant_BP_terms': len(bp), 'top3_GO_IDs': ';'.join((r['GO_ID'] for r in bp[:3])), 'top3_terms': '; '.join((r['GO_term'] for r in bp[:3]))})
    df = pd.DataFrame(rows).sort_values(['community', 'FDR', 'p', 'GO_ID'])
    df['community_testable'] = df.community.map({r['community']: r['testable'] for r in summary})
    df.to_csv(out / 'go_enrichment.tsv', sep='\t', index=False)
    sig = df[(df.FDR < config['go']['fdr']) & df.community_testable]
    sig.to_csv(out / 'enrichment_significant.tsv', sep='\t', index=False)
    bp = sig[sig.aspect == 'P']
    bp.groupby('community').head(3).to_csv(out / 'top3_BP.tsv', sep='\t', index=False)
    pd.DataFrame(summary).to_csv(out / 'community_go_summary.tsv', sep='\t', index=False)
    stats = {'background_genes': N, 'mapped_background_genes': len(mapped), 'family_terms': len(family), 'total_audit_tests': len(df), 'significant_all_aspects': len(sig), 'significant_BP': len(bp), 'communities_with_significant_BP': int(bp.community.nunique()), 'annotation_exclusions': dict(exclusions), 'qualifiers_observed': dict(relation_counts), 'compatibility': 'All positive LCC annotation terms active, resolvable and namespace-consistent; direct annotations, no relation propagation', 'component_dates': {'GAF_retrieved_at': gm.get('retrieved_at'), 'OBO_retrieved_at': om.get('retrieved_at')}}
    write_json(out / 'summary.json', stats)
    write_json(out / 'source_manifest.json', {'sources': [gm, om], 'compatibility': stats['compatibility']})
    pd.DataFrame([{'GO_ID': t, 'GO_term': terms[t]['name'], 'aspect': aspects[terms[t]['namespace']], 'LCC_annotated_genes': len(associations[t])} for t in family]).to_csv(out / 'network_GO_annotation_summary.tsv', sep='\t', index=False)
    manifest['database_snapshots'].extend([gm, om])
    manifest['stages']['04'] = {'status': 'SUCCESS', 'summary': stats}
    return stats
if __name__ == '__main__':
    from pipeline import main
    main(default_stage=4)
