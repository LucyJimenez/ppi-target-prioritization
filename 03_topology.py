"""Stage 3: topology. Run --help for usage."""
from common import table, write_json
'One physical topology/LCC, primary analysis and paired support sensitivity.'
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

def canonical(p):
    groups = defaultdict(list)
    for n, c in p.items():
        groups[c].append(n)
    order = sorted(groups, key=lambda c: (-len(groups[c]), min(groups[c])))
    remap = {c: i + 1 for i, c in enumerate(order)}
    return {n: remap[c] for n, c in p.items()}

def ensemble(g, cfg, out, weight='support', resolution=1.0, prefix='primary'):
    nodes = sorted(g)
    seeds = range(cfg['community_seeds'])
    parts = [canonical(cl.best_partition(g, weight=weight, resolution=resolution, random_state=s)) for s in seeds]
    labels = np.array([[p[n] for n in nodes] for p in parts])
    co = np.zeros((len(nodes), len(nodes)), dtype=float)
    for lab in labels:
        co += (lab[:, None] == lab[None, :]) / len(parts)
    cg = nx.Graph()
    cg.add_nodes_from(nodes)
    ii, jj = np.where(np.triu(co, 1) >= cfg['consensus_threshold'] - 1e-12)
    cg.add_weighted_edges_from(((nodes[i], nodes[j], float(co[i, j])) for i, j in zip(ii, jj)))
    cons = [canonical(cl.best_partition(cg, weight='weight', resolution=1.0, random_state=s)) for s in seeds]
    primary = cons[0]
    ari = np.array([[adjusted_rand_score(x, y) for y in labels] for x in labels])
    means = (ari.sum(axis=1) - 1) / (len(parts) - 1) if len(parts) > 1 else np.ones(1)
    medoid = int(np.argmax(means))
    table([{'seed': s, 'resolution': resolution, 'communities': len(set(p.values())), 'modularity': cl.modularity(p, g, weight=weight), 'ARI_to_consensus': adjusted_rand_score(labels[s], [primary[n] for n in nodes]), 'mean_pairwise_ARI': means[s], 'medoid': s == medoid} for s, p in enumerate(parts)], out / f'{prefix}_seeds.tsv')
    table([{'WBGene': n, **{f'seed_{s}': p[n] for s, p in enumerate(parts)}} for n in nodes], out / f'{prefix}_partitions.tsv')
    table([{'WBGene': n, **{f'seed_{s}': p[n] for s, p in enumerate(cons)}} for n in nodes], out / f'{prefix}_consensus_partitions.tsv')
    np.savez_compressed(out / f'{prefix}_coassignment.npz', nodes=np.array(nodes), coassignment=co, pairwise_ARI=ari)
    summary = {'communities': len(set(primary.values())), 'modularity': cl.modularity(primary, g, weight=weight), 'mean_seed_pairwise_ARI': float(means.mean()), 'medoid_seed': medoid, 'seed_communities_min': min((len(set(p.values())) for p in parts)), 'seed_communities_max': max((len(set(p.values())) for p in parts)), 'consensus_distinct_partitions': len({tuple((p[n] for n in nodes)) for p in cons})}
    return (primary, summary)

def metrics(g, cfg, weighted):
    weight = 'support' if weighted else None
    distance = 'distance' if weighted else None
    methods = {'degree': dict(g.degree()), 'degree_centrality': nx.degree_centrality(g), 'strength': dict(g.degree(weight=weight)), 'betweenness': nx.betweenness_centrality(g, normalized=True, endpoints=False, weight=distance), 'closeness': nx.closeness_centrality(g, distance=distance, wf_improved=True), 'PageRank': nx.pagerank(g, alpha=cfg['pagerank_alpha'], max_iter=cfg['pagerank_max_iter'], tol=cfg['pagerank_tol'], weight=weight), 'eigenvector': nx.eigenvector_centrality(g, max_iter=cfg['eigen_max_iter'], tol=cfg['eigen_tol'], weight=weight), 'kcore': nx.core_number(g)}
    df = pd.DataFrame({'WBGene': sorted(g)})
    for name, values in methods.items():
        df[name] = df.WBGene.map(values)
        df[name + '_rank'] = df[name].rank(ascending=False, method='average')
        df[name + '_percentile'] = df[name].rank(ascending=True, pct=True, method='average')
    return df

def robustness(g, cfg, out):
    nodes = sorted(g)
    edges = sorted((tuple(sorted(e)) for e in g.edges()))
    rng = np.random.default_rng(cfg['seed'])
    values = {mode: dict(g.degree(weight='support' if mode == 'weighted' else None)) for mode in ('primary', 'weighted')}
    thresholds = {mode: float(np.quantile(list(v.values()), cfg['hub_quantile'], method='linear')) for mode, v in values.items()}
    counts = {mode: Counter() for mode in values}
    detail = []
    selections = []
    for rep in range(cfg['hub_replicates']):
        indices = sorted(rng.choice(len(edges), size=max(1, round(cfg['edge_retention'] * len(edges))), replace=False).tolist())
        selections.append({'replicate': rep, 'edge_indices': indices})
        h = nx.Graph()
        h.add_nodes_from(nodes)
        h.add_edges_from(((edges[i][0], edges[i][1], g.edges[edges[i]]) for i in indices))
        for mode in values:
            scores = dict(h.degree(weight='support' if mode == 'weighted' else None))
            threshold = float(np.quantile(list(scores.values()), cfg['hub_quantile'], method='linear'))
            for n in nodes:
                hub = scores[n] >= threshold
                counts[mode][n] += int(hub)
                detail.append({'replicate': rep, 'mode': mode, 'WBGene': n, 'score': scores[n], 'threshold': threshold, 'hub': hub})
    write_json(out / 'perturbation_edge_subsets.json', {'edge_order': edges, 'selections': selections})
    table(detail, out / 'hub_perturbations.tsv')
    frames = {}
    for mode in values:
        frames[mode] = pd.DataFrame([{'WBGene': n, 'hub_threshold': thresholds[mode], 'hub_original': values[mode][n] >= thresholds[mode], 'hub_frequency': counts[mode][n] / cfg['hub_replicates'], 'robust_hub': bool(values[mode][n] >= thresholds[mode] and counts[mode][n] / cfg['hub_replicates'] >= cfg['robust_fraction'])} for n in nodes])
    return (frames, thresholds)

def execute(root, config, run, manifest):
    cfg = config['network']
    out = run / '03_network_analysis'
    out.mkdir()
    g = nx.read_graphml(run / '02_network/network.graphml')
    lcc = nx.read_graphml(run / '02_network/lcc.graphml')
    nodes = sorted(lcc)
    nt = pd.read_csv(run / '02_network/nodes.tsv', sep='\t', keep_default_na=False)
    inputs = set(nt.loc[nt.original_input.astype(str).str.lower() == 'true', 'WBGene'])
    components = list(nx.connected_components(g))
    hubs, thresholds = robustness(lcc, cfg, out)
    results = {}
    summary = {}
    for mode in ('primary', 'weighted'):
        print('Network:', mode, 'metrics and ensemble', flush=True)
        partition, cs = ensemble(lcc, cfg, out, weight='support' if mode == 'weighted' else 'unused_unweighted', prefix=mode)
        df = metrics(lcc, cfg, mode == 'weighted').merge(hubs[mode], on='WBGene', validate='one_to_one').merge(nt, on='WBGene', validate='one_to_one')
        df['community'] = df.WBGene.map(partition)
        df.to_csv(out / ('node_metrics.tsv' if mode == 'primary' else 'weighted_metrics.tsv'), sep='\t', index=False)
        results[mode] = df
        summary[mode] = cs
        summary[mode].update(hub_threshold=thresholds[mode], original_hubs=int(df.hub_original.sum()), robust_hubs=int(df.robust_hub.sum()))
    primary = results['primary']
    primary[primary.robust_hub].sort_values(['degree', 'WBGene'], ascending=[False, True]).to_csv(out / 'hubs.tsv', sep='\t', index=False)
    table([{'community': int(c), 'nodes': len(d)} for c, d in primary.groupby('community')], out / 'communities.tsv')
    compare = []
    k = int(primary.hub_original.sum())
    for metric in ('degree', 'strength', 'betweenness', 'closeness', 'PageRank', 'eigenvector', 'kcore'):
        x = results['primary'].set_index('WBGene').loc[nodes, metric]
        y = results['weighted'].set_index('WBGene').loc[nodes, metric]
        a = set(x.sort_values(ascending=False, kind='stable').head(k).index)
        b = set(y.sort_values(ascending=False, kind='stable').head(k).index)
        at = set(x[x >= x.nlargest(k).min()].index)
        bt = set(y[y >= y.nlargest(k).min()].index)
        compare.append({'metric': metric, 'Spearman': float(spearmanr(x, y).statistic), 'k': k, 'exact_topk_overlap': len(a & b), 'exact_topk_Jaccard': len(a & b) / len(a | b), 'tie_inclusive_primary_n': len(at), 'tie_inclusive_weighted_n': len(bt), 'tie_inclusive_Jaccard': len(at & bt) / len(at | bt)})
    table(compare, out / 'weighted_comparison.tsv')
    ari = float(adjusted_rand_score(primary.community, results['weighted'].community))
    pd.crosstab(primary.community, results['weighted'].community, rownames=['primary'], colnames=['weighted']).to_csv(out / 'community_contingency.tsv', sep='\t')
    sweep = []
    for res in cfg['resolution_sweep']:
        print('Resolution diagnostic:', res, flush=True)
        pp, ss = ensemble(lcc, cfg, out, weight='unused_unweighted', resolution=res, prefix=f'resolution_{res}')
        sweep.append({'resolution': res, **ss, 'ARI_to_primary': adjusted_rand_score(primary.community, [pp[n] for n in nodes])})
    table(sweep, out / 'resolution_sensitivity.tsv')
    deletion = []
    for field in ('publications', 'experimental_systems'):
        categories = sorted({v for _, _, d in lcc.edges(data=True) for v in d[field].split(';')})
        for category in categories:
            h = lcc.copy()
            for a, b, d in list(h.edges(data=True)):
                remaining = set(d[field].split(';')) - {category}
                if not remaining:
                    h.remove_edge(a, b)
                elif field == 'publications':
                    h[a][b]['support'] = len(remaining)
            score = dict(h.degree())
            q = float(np.quantile(list(score.values()), cfg['hub_quantile']))
            wh = dict(h.degree(weight='support'))
            wq = float(np.quantile(list(wh.values()), cfg['hub_quantile']))
            for n in nodes:
                deletion.append({'deleted_field': field, 'deleted_value': category, 'WBGene': n, 'degree': score[n], 'degree_hub': score[n] >= q, 'strength': wh[n], 'strength_hub': wh[n] >= wq, 'remaining_edges': h.number_of_edges()})
    table(deletion, out / 'source_deletion_sensitivity.tsv')
    table([{'quantile': q, 'threshold': float(np.quantile(primary.degree, q)), 'hubs': int((primary.degree >= np.quantile(primary.degree, q)).sum())} for q in cfg['hub_quantile_sweep']], out / 'hub_quantile_sensitivity.tsv')
    support_rows = []
    for n in nodes:
        pubs = set()
        systems = set()
        records = set()
        for _, _, d in lcc.edges(n, data=True):
            pubs.update(d['publications'].split(';'))
            systems.update(d['experimental_systems'].split(';'))
            records.update(d['interaction_ids'].split(';'))
        support_rows.append({'WBGene': n, 'degree': lcc.degree(n), 'publications': len(pubs), 'systems': len(systems), 'records': len(records)})
    bias = table(support_rows, out / 'study_bias.tsv')
    summary.update(full_nodes=len(g), full_edges=g.number_of_edges(), isolates=len(list(nx.isolates(g))), components=len(components), lcc_nodes=len(lcc), lcc_edges=lcc.number_of_edges(), input_in_LCC=len(inputs & set(lcc)), density=nx.density(lcc), mean_degree=2 * lcc.number_of_edges() / len(lcc), average_clustering=nx.average_clustering(lcc), transitivity=nx.transitivity(lcc), diameter=nx.diameter(lcc), average_shortest_path=nx.average_shortest_path_length(lcc), community_ARI_primary_weighted=ari, same_nodes_edges_LCC=True, robust_hub_overlap=len(set(primary[primary.robust_hub].WBGene) & set(results['weighted'][results['weighted'].robust_hub].WBGene)), study_bias_Spearman={col: float(spearmanr(bias.degree, bias[col]).statistic) for col in ('publications', 'systems', 'records')}, optional_Leiden='not_requested_not_run')
    write_json(out / 'summary.json', summary)
    (out / 'network_sensitivity_report.md').write_text('# Weighted sensitivity\n\nSame base nodes, edges and LCC. Support is distinct publications, not affinity.\n\n' + json.dumps(summary, indent=2) + '\n\nSee weighted_comparison.tsv, community_contingency.tsv and paired perturbation outputs. No materiality cutoff imposed.\n')
    table([{'WBGene': n, 'community': c} for n, c in zip(results['weighted'].WBGene, results['weighted'].community)], out / 'weighted_communities.tsv')
    primary[['WBGene', 'hub_original', 'hub_frequency', 'robust_hub', 'hub_threshold']].to_csv(out / 'hub_robustness.tsv', sep='\t', index=False)
    write_json(out / 'sensitivity_summary.json', summary)
    manifest['stages']['03'] = {'status': 'SUCCESS', 'summary': summary}
    return summary
if __name__ == '__main__':
    from pipeline import main
    main(default_stage=3)
