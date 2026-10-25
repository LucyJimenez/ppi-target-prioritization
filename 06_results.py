"""Stage 6: results. Run --help for usage."""
from common import sha256, write_json
'Offline figures; all inputs are saved local analysis outputs.'
import json, textwrap
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

def generate_figures(root, config, run, out):
    figdir = out / 'figures'
    figdir.mkdir()
    md = out / 'figure_manifests'
    md.mkdir()
    installed = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in config['plotting']['fonts'] if f in installed), 'DejaVu Sans')
    plt.rcParams.update({'font.family': font, 'font.size': 8, 'axes.labelsize': 9, 'axes.titlesize': 10, 'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'savefig.dpi': config['plotting']['dpi']})
    nodes = pd.read_csv(run / '03_network_analysis/node_metrics.tsv', sep='\t', keep_default_na=False)
    g = nx.read_graphml(run / '02_network/lcc.graphml')
    cfg = config['plotting']
    width = cfg['width_mm'] / 25.4
    cmap = plt.get_cmap('tab20')
    colors = {c: cmap((int(c) - 1) % 20) for c in nodes.community.unique()}
    manifests = []

    def save(fig, name, caption, inputs):
        fig.savefig(figdir / (name + '.svg'))
        fig.savefig(figdir / (name + '.pdf'))
        fig.savefig(figdir / (name + '.png'), dpi=cfg['dpi'])
        plt.close(fig)
        meta = {'figure': name, 'caption': caption, 'font': font, 'width_mm': cfg['width_mm'], 'dpi': cfg['dpi'], 'formats': ['svg', 'pdf', 'png'], 'inputs': [{'path': str(f.relative_to(root)), 'sha256': sha256(f)} for f in inputs], 'parameters': cfg, 'external_queries': False}
        write_json(md / (name + '.json'), meta)
        manifests.append(meta)
    coords = nx.spring_layout(g, seed=cfg['seed'], k=cfg['spring_k'], iterations=cfg['iterations'], threshold=1e-05, weight=None, scale=1, method='force')
    pd.DataFrame([{'WBGene': n, 'x': coords[n][0], 'y': coords[n][1]} for n in sorted(coords)]).to_csv(figdir / 'network_layout.tsv', sep='\t', index=False)
    ordered = nodes.set_index('WBGene').loc[list(g)]
    hub = ordered.robust_hub.astype(bool)
    fig, ax = plt.subplots(figsize=(width, 7.5))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.9, bottom=0.12)
    nx.draw_networkx_edges(g, coords, ax=ax, edge_color='#a1aeb8', alpha=0.38, width=0.35)
    nx.draw_networkx_nodes(g, coords, ax=ax, node_color=[colors[c] for c in ordered.community], node_size=np.where(hub, 42, 9), edgecolors=np.where(hub, '#101e2c', 'white'), linewidths=np.where(hub, 0.65, 0.15))
    top = nodes[nodes.robust_hub].sort_values(['degree', 'WBGene'], ascending=[False, True]).head(12)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    placed = []
    for r in top.itertuples():
        label = r.Gene_symbol or r.WBGene
        for radius in (12, 22, 34, 48, 64):
            found = False
            for angle in np.linspace(0, 2 * np.pi, 12, endpoint=False):
                offset = (radius * np.cos(angle), radius * np.sin(angle))
                artist = ax.annotate(label, xy=coords[r.WBGene], xytext=offset, textcoords='offset points', fontsize=7, fontfamily=font, ha='center', va='center', bbox={'facecolor': 'white', 'alpha': 0.92, 'edgecolor': 'none', 'pad': 0.8}, arrowprops={'arrowstyle': '-', 'color': '#526b79', 'lw': 0.4})
                bounds = artist.get_window_extent(renderer).expanded(1.05, 1.08)
                if not any((bounds.overlaps(b) for b in placed)):
                    placed.append(bounds)
                    found = True
                    break
                artist.remove()
            if found:
                break
    ax.set_title('A  Physical interaction landscape', loc='left', fontweight='bold', pad=16)
    ax.axis('off')
    handles = [plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=colors[c], markeredgecolor='none', label=f'C{c}', markersize=5) for c in sorted(colors)]
    fig.legend(handles=handles, loc='lower center', ncol=9, frameon=False, bbox_to_anchor=(0.5, 0.035))
    caption = f'Largest connected component: {len(g):,} WBGene nodes and {g.number_of_edges():,} physical edges. Colors show primary unweighted communities; outlined larger nodes are robust hubs. Labels identify {len(placed)} high-degree robust hubs. Layout is deterministic and unweighted.'
    save(fig, 'figure_A_network', caption, [run / '02_network/lcc.graphml', run / '03_network_analysis/node_metrics.tsv'])
    go = pd.read_csv(run / '04_functional_context/community_go_summary.tsv', sep='\t', keep_default_na=False).sort_values('community')
    topgo = pd.read_csv(run / '04_functional_context/top3_BP.tsv', sep='\t')
    wrapped = {int(r.community): [(textwrap.fill(str(t.GO_term), 59), float(t.FDR)) for t in topgo[topgo.community == r.community].itertuples()] for r in go.itertuples()}
    heights = [max(0.38, 0.13 * sum((label.count('\n') + 1 for label, _ in wrapped[int(c)])) + 0.14) for c in go.community]
    y = np.cumsum(heights) - np.array(heights) / 2
    total = sum(heights)
    fig = plt.figure(figsize=(width, total + 0.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 0.8, 4.1], left=0.075, right=0.98, top=0.91, bottom=0.85 / (total + 0.9), wspace=0.18)
    ax1, ax2, ax3 = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax1.barh(y, go.nodes, color=[colors[c] for c in go.community], height=np.minimum(heights, 0.44) * 0.75)
    ax1.set_yticks(y, [f'C{c}' for c in go.community])
    ax1.set_xlabel('Nodes')
    ax2.barh(y, go.significant_BP_terms, color='#237e87', height=np.minimum(heights, 0.44) * 0.75)
    ax2.set_yticks([])
    ax2.set_xlabel('GO BP terms')
    for ax in (ax1, ax2):
        ax.set_ylim(total, 0)
        ax.spines[['top', 'right']].set_visible(False)
    ax3.set_xlim(0, 1)
    ax3.set_ylim(total, 0)
    ax3.axis('off')
    maxlog = max(2.0, float(-np.log10(topgo.FDR.clip(lower=1e-300)).max())) if len(topgo) else 2.0
    norm = matplotlib.colors.Normalize(vmin=0, vmax=maxlog)
    term_cmap = plt.get_cmap('viridis_r')
    for i, r in enumerate(go.itertuples()):
        labels = wrapped[int(r.community)]
        if not labels:
            ax3.text(0, y[i], 'No significant BP term', fontsize=7, color='#6e7f8c', va='center', style='italic')
        else:
            line_count = sum((label.count('\n') + 1 for label, _ in labels))
            at = y[i] - 0.13 * line_count / 2
            for label, fdr in labels:
                ax3.text(0, at, label, fontsize=7, va='top', linespacing=1.1, color=term_cmap(norm(-np.log10(max(fdr, 1e-300)))))
                at += 0.13 * (label.count('\n') + 1)
    fig.suptitle('B  Community size and functional context', x=0.075, ha='left', fontsize=11, fontweight='bold')
    ax3.set_title('Top three BP terms (FDR, p, GO ID order)', fontsize=8, loc='left')
    cax = fig.add_axes([0.66, 0.45 / (total + 0.9), 0.25, 0.06 / (total + 0.9)])
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=term_cmap), cax=cax, orientation='horizontal', label='Term color: -log10(FDR)')
    save(fig, 'figure_B_community_GO', 'Community node counts, significant Biological Process terms (within-community BH FDR <0.05 across all tested aspects), and up to three ranked BP terms. No significant term is an explicit result, not missing data.', [run / '04_functional_context/community_go_summary.tsv', run / '04_functional_context/top3_BP.tsv'])
    weighted = pd.read_csv(run / '03_network_analysis/weighted_metrics.tsv', sep='\t')
    merged = nodes.merge(weighted[['WBGene', 'degree', 'strength', 'community', 'robust_hub']], on='WBGene', suffixes=('', '_weighted'), validate='one_to_one')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(width, 3.2), layout='constrained')
    ax1.scatter(merged.degree, merged.strength_weighted, s=np.where(merged.robust_hub, 22, 6), alpha=0.7, c=['#176e78' if x else '#a7bac3' for x in merged.robust_hub])
    ax1.set(xlabel='Unweighted degree', ylabel='Publication-weighted strength')
    ax1.set_title('C1  Same topology, different support', loc='left', fontsize=9)
    matrix = pd.crosstab(merged.community, merged.community_weighted)
    im = ax2.imshow(matrix, aspect='auto', cmap='Blues')
    ax2.set(xlabel='Weighted community', ylabel='Primary community')
    ax2.set_xticks(np.arange(len(matrix.columns))[::2], matrix.columns[::2])
    ax2.set_yticks(np.arange(len(matrix.index))[::2], matrix.index[::2])
    ax2.set_title('C2  Membership agreement', loc='left', fontsize=9)
    fig.colorbar(im, ax=ax2, shrink=0.65, label='Genes')
    save(fig, 'figure_C_sensitivity', 'Sensitivity uses identical nodes and edges; publication count is documentary support, not biochemical affinity. The contingency matrix compares membership independently of community label numbering.', [run / '03_network_analysis/node_metrics.tsv', run / '03_network_analysis/weighted_metrics.tsv'])
    write_json(md / 'index.json', manifests)
'Available-evidence profiles; no composite importance score or druggability claim.'
import json
import pandas as pd

def integrate(hubs, statuses, summaries):
    result = hubs.copy()
    for source in ('wormbase', 'bgee', 'blast'):
        entry = next((r for r in statuses if r['source'] == source), {'status': 'NOT_RUN', 'reason': 'Not selected'})
        result[source + '_status'] = entry['status']
        result[source + '_status_reason'] = entry.get('reason', '')
        frame = summaries.get(source)
        if entry['status'] in ('SUCCESS', 'PARTIAL') and frame is not None:
            frame = frame.drop(columns=['UniProt', 'Orthogroup', source + '_status'], errors='ignore')
            if frame.WBGene.duplicated().any():
                raise ValueError('Duplicate candidate summary for ' + source)
            result = result.merge(frame, on='WBGene', how='left', validate='one_to_one')
        for col in {'wormbase': ['wormbase_observed_records', 'wormbase_not_observed_records', 'wormbase_essentiality'], 'bgee': ['bgee_conditions', 'bgee_present_conditions'], 'blast': ['blast_hits', 'best_human_UniProt', 'best_human_gene', 'best_evalue', 'best_query_coverage', 'best_subject_coverage']}[source]:
            if col not in result:
                result[col] = pd.NA
    return result

def execute(root, config, run, manifest):
    out = run / '06_target_profiles'
    out.mkdir()
    hubs = pd.read_csv(run / '03_network_analysis/hubs.tsv', sep='\t', keep_default_na=False)
    status = pd.read_csv(run / '05_biological_evidence/evidence_status.tsv', sep='\t', keep_default_na=False).to_dict('records')
    summaries = {}
    for source in ('wormbase', 'bgee', 'blast'):
        file = run / '05_biological_evidence' / source / 'candidate_summary.tsv'
        if file.exists():
            summaries[source] = pd.read_csv(file, sep='\t')
    master = integrate(hubs, status, summaries)
    go = pd.read_csv(run / '04_functional_context/community_go_summary.tsv', sep='\t', keep_default_na=False)
    master = master.merge(go[['community', 'significant_BP_terms', 'top3_GO_IDs', 'top3_terms']], on='community', how='left', validate='many_to_one')
    master.to_csv(out / 'candidate_target_profiles.tsv', sep='\t', index=False, na_rep='NA')
    completeness = master[['WBGene', 'UniProt', 'Orthogroup', 'wormbase_status', 'bgee_status', 'blast_status']].copy()
    completeness.to_csv(out / 'evidence_completeness.tsv', sep='\t', index=False)
    net = json.loads((run / '03_network_analysis/summary.json').read_text())
    gos = json.loads((run / '04_functional_context/summary.json').read_text())
    bg = summaries.get('bgee', pd.DataFrame())
    bl = summaries.get('blast', pd.DataFrame())
    summary = {'candidate_count': len(master), 'network': net, 'GO': gos, 'evidence': status, 'Bgee_hubs_with_records': int((bg.bgee_conditions > 0).sum()) if 'bgee_conditions' in bg else None, 'BLAST_hubs_with_hits': int((bl.blast_hits > 0).sum()) if 'blast_hits' in bl else None, 'BLAST_hubs_without_hits': int((bl.blast_hits == 0).sum()) if 'blast_hits' in bl else None, 'orthogroup': manifest.get('orthogroup_annotation', {}), 'interpretation': 'Candidate target profiles, not validated targets or druggability predictions'}
    write_json(out / 'analysis_summary.json', summary)
    generate_figures(root, config, run, out)
    write_report(root, config, run, out, summary, master)
    manifest['stages']['06'] = {'status': 'SUCCESS', 'summary': {'profiles': len(master), 'evidence_incomplete': any((r['status'] != 'SUCCESS' for r in status))}}
    return summary

def write_report(root, config, run, out, summary, profiles):
    """Render actual saved results, never hard-coded database dates or study conclusions."""
    import base64, html
    blocks = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>PPI results</title>', '<style>body{font:16px Arial;max-width:1100px;margin:40px auto;padding:20px}img{max-width:100%}table{border-collapse:collapse}td,th{padding:8px;border:1px solid #ddd}.table{overflow:auto}</style>', '<h1>Physical PPI pipeline results</h1>', '<p>Network nodes are WBGene IDs. Missing evidence is not negative evidence; sequence similarity does not establish orthology or druggability.</p>']
    blocks.append('<h2>Evidence source status</h2>' + pd.DataFrame(summary['evidence']).to_html(index=False, escape=True))
    blocks.append('<h2>Network and GO summary</h2><pre>' + html.escape(json.dumps({'network': summary['network'], 'GO': summary['GO']}, indent=2)) + '</pre>')
    for name in ['figure_A_network', 'figure_B_community_GO', 'figure_C_sensitivity']:
        image = base64.b64encode((out / 'figures' / (name + '.svg')).read_bytes()).decode()
        blocks.append('<img alt="' + name + '" src="data:image/svg+xml;base64,' + image + '">')
    blocks.append('<h2>Candidate profiles</h2><div class="table">' + profiles.to_html(index=False, escape=True, na_rep='NA') + '</div></html>')
    (out / 'analysis_report.html').write_text('\n'.join(blocks))
if __name__ == '__main__':
    from pipeline import main
    main(default_stage=6)
