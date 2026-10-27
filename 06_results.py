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
    """Readable overview and compact candidate rows with complete accessible details."""
    import base64, html
    def esc(value):
        return html.escape(str(value), quote=True)
    def value(x):
        if x is None or pd.isna(x) or str(x) in ('', 'NA', 'nan'):
            return 'Not available'
        if isinstance(x, (float, np.floating)):
            return f'{x:.4g}'
        return str(x)
    def number(x):
        return f'{int(x):,}' if x is not None else 'Not available'
    def simple_table(headers, rows):
        return '<div class="scroll"><table><thead><tr>'+''.join('<th scope="col">'+esc(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(value(x))+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'
    n=summary['network'];g=summary['GO'];primary=n['primary'];weighted=n['weighted']
    styles='''
    :root{color-scheme:light;--ink:#193344;--muted:#536a79;--line:#dce6eb;--accent:#087d81}*{box-sizing:border-box}
    body{font:16px/1.55 system-ui,-apple-system,Arial,sans-serif;color:var(--ink);background:#f3f7f8;margin:0}
    main{max-width:1240px;margin:auto;padding:36px 24px 70px}h1{font-size:clamp(27px,4vw,40px);line-height:1.2}h2{font-size:25px;margin:0 0 14px}h3{font-size:18px}
    section{background:white;border:1px solid var(--line);border-radius:12px;padding:26px;margin:26px 0;scroll-margin-top:20px}
    p{max-width:980px}a{color:#006b73}nav{display:flex;gap:18px;flex-wrap:wrap;margin:24px 0}.muted{color:var(--muted);font-size:14px}
    .cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px}.card{padding:18px;background:#edf6f6;border-radius:8px}.card strong{display:block;font-size:32px}.card span{font-size:14px}
    .panels{display:grid;grid-template-columns:1fr 1fr;gap:24px}.scroll{overflow:auto;max-width:100%}table{width:100%;border-collapse:separate;border-spacing:0;font-size:14px;text-align:left}
    th{background:#edf3f5;color:#354e5e;font-weight:600}th,td{padding:12px 14px;border-bottom:1px solid var(--line);vertical-align:middle}tbody tr:nth-child(even){background:#f8fafb}
    .profiles{min-width:1080px;table-layout:fixed}.profiles tbody tr{height:64px}.profiles td{height:64px;padding:8px 12px}.profiles th{position:sticky;top:0;z-index:1}
    .profiles th:first-child,.profiles td:first-child{position:sticky;left:0;background:#f0f6f7;z-index:2}.profiles th:first-child{z-index:3}
    .profiles-scroll{max-height:620px;border:1px solid var(--line);border-radius:8px}.gene{display:block;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;font-weight:650}.id{display:block;font-size:12px;color:var(--muted);white-space:nowrap}
    button,.download{display:inline-block;font:inherit;font-size:14px;border:1px solid #7ca5ad;border-radius:6px;background:white;color:#005c66;padding:7px 12px;cursor:pointer;text-decoration:none}button:hover,.download:hover{background:#e5f2f3}button:focus-visible,input:focus-visible,summary:focus-visible{outline:3px solid #e6a12c;outline-offset:3px}
    .toolbar{display:flex;gap:16px;align-items:end;justify-content:space-between;flex-wrap:wrap;margin:18px 0}.toolbar label{display:block;font-size:14px;font-weight:600}input{display:block;width:min(360px,75vw);font:inherit;padding:10px 12px;border:1px solid #9aafb9;border-radius:6px;margin-top:5px}
    dialog{border:1px solid var(--line);border-radius:12px;max-width:850px;width:calc(100% - 32px);max-height:85vh;padding:0;color:var(--ink)}dialog::backdrop{background:#142c40aa}.dialog-head{position:sticky;top:0;background:white;padding:18px 24px;border-bottom:1px solid var(--line);display:flex;justify-content:space-between;gap:16px;align-items:center}.dialog-head h3{margin:0}.detail-body{padding:12px 24px 24px}dl{display:grid;grid-template-columns:minmax(140px,35%) 1fr;gap:0;margin:0}dt,dd{padding:10px 0;border-bottom:1px solid var(--line);overflow-wrap:anywhere}dt{color:var(--muted);padding-right:18px}dd{margin:0}
    details{margin-top:20px}summary{cursor:pointer;color:#006b73;font-weight:600}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f5f7f8;padding:16px;font-size:12px}figure{margin:22px 0}img{max-width:100%;height:auto}figcaption{font-size:14px;color:var(--muted)}.empty{padding:20px}.pill{font-size:12px;display:inline-block;border-radius:12px;padding:3px 8px;background:#e6f2eb}
    @media(max-width:760px){main{padding:20px 12px}section{padding:18px}.cards{grid-template-columns:repeat(2,1fr)}.panels{grid-template-columns:1fr}dl{grid-template-columns:1fr}dt{border:0;padding-bottom:0}dd{padding-top:2px}}
    @media print{body{background:white}main{max-width:none}.profiles-scroll{max-height:none}.toolbar,button{display:none}section{break-inside:avoid}details{display:none}}
    '''
    blocks=['<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PPI results</title><style>'+styles+'</style></head><body><main>',
        '<h1>Physical PPI pipeline results</h1><p>Network nodes are WBGene IDs. Missing evidence is not negative evidence; sequence similarity does not establish orthology or druggability.</p>',
        '<nav aria-label="Report sections"><a href="#overview">Network &amp; GO</a><a href="#figures">Figures</a><a href="#candidates">Candidate profiles</a><a href="#sources">Evidence sources</a></nav>',
        '<section id="overview"><h2>Network and GO summary</h2><p class="muted">The largest connected component (LCC) is the connected part of the network used for community and hub analysis.</p><div class="cards">']
    for label,count in [('LCC genes',n['lcc_nodes']),('Physical interactions',n['lcc_edges']),('Communities',primary['communities']),('Robust hubs',primary['robust_hubs'])]:
        blocks.append('<div class="card"><strong>'+number(count)+'</strong><span>'+label+'</span></div>')
    blocks.append('</div><div class="panels"><div><h3>Network structure</h3>'+simple_table(['Measure','Result'],[
        ('Full network genes',number(n['full_nodes'])),('Full network interactions',number(n['full_edges'])),('Input genes in the LCC',number(n['input_in_LCC'])),('Primary modularity',round(primary['modularity'],3)),('Weighted communities',weighted['communities']),('Hubs shared by both analyses',n.get('robust_hub_overlap'))])+'</div>')
    blocks.append('<div><h3>GO functional context</h3>'+simple_table(['Measure','Result'],[
        ('LCC genes with annotations',f"{g['mapped_background_genes']:,} / {g['background_genes']:,}"),('Significant community–BP associations',g['significant_BP']),('Communities with significant BP terms',g['communities_with_significant_BP']),('Tested GO terms per community',g['family_terms'])])+ '<p class="muted">BP = Biological Process. Associations count a term separately in each enriched community; they are not necessarily unique GO terms.</p></div></div>')
    blocks.append('<details><summary>Show all network and GO statistics</summary><pre>'+esc(json.dumps({'network':n,'GO':g},indent=2))+'</pre></details></section>')
    blocks.append('<section id="figures"><h2>Network, communities and sensitivity</h2>')
    for name,caption in [('figure_A_network','Physical interaction network: colors identify communities; outlined nodes are robust hubs.'),('figure_B_community_GO','Community sizes and significant Biological Process enrichment.'),('figure_C_sensitivity','Comparison of the primary and publication-weighted analyses.')]:
        image=base64.b64encode((out/'figures'/(name+'.svg')).read_bytes()).decode()
        blocks.append('<figure><img alt="'+esc(caption)+'" src="data:image/svg+xml;base64,'+image+'"><figcaption>'+esc(caption)+'</figcaption></figure>')
    blocks.append('</section><section id="candidates"><h2>Candidate profiles</h2><p>Compare the key measures below. Open <b>Details</b> for all identifiers, metrics and source evidence.</p><p class="muted">Phenotype counts are observed records, not an essentiality classification. Human hit counts describe sequence similarity.</p>')
    csv64=base64.b64encode(profiles.to_csv(index=False).encode()).decode()
    blocks.append('<div class="toolbar"><label for="candidate-search">Find a candidate<input id="candidate-search" type="search" placeholder="Gene, WBGene or UniProt"></label><a class="download" download="candidate_profiles.csv" href="data:text/csv;base64,'+csv64+'">Download complete table</a></div><p id="candidate-count" class="muted" aria-live="polite">'+str(len(profiles))+' candidates</p>')
    blocks.append('<div class="scroll profiles-scroll" role="region" aria-label="Candidate profiles; scroll horizontally for more columns" tabindex="0"><table class="profiles"><colgroup><col style="width:210px">'+''.join('<col style="width:'+str(w)+'px">' for w in [95,85,105,130,130,115,100])+'</colgroup><thead><tr>'+''.join('<th scope="col">'+x+'</th>' for x in ['Gene / WBGene','Community','Degree','Hub stability','Observed phenotypes','Expression conditions','Human hits','Full profile'])+'</tr></thead><tbody>')
    dialogs=[]
    for i,row in enumerate(profiles.to_dict('records')):
        gene=value(row.get('Gene_symbol'));gene=row['WBGene'] if gene=='Not available' else gene
        search=' '.join(str(row.get(k,'')) for k in ['Gene_symbol','WBGene','UniProt','Orthogroup']).lower()
        frequency=row.get('hub_frequency');stability='Not available' if frequency is None or pd.isna(frequency) else f'{float(frequency)*100:.0f}%'
        blocks.append('<tr data-search="'+esc(search)+'"><td><span class="gene" title="'+esc(gene)+'">'+esc(gene)+'</span><span class="id">'+esc(row['WBGene'])+'</span></td>'+''.join('<td>'+esc(value(v))+'</td>' for v in [row.get('community'),row.get('degree'),stability,row.get('wormbase_observed_records'),row.get('bgee_conditions'),row.get('blast_hits')])+'<td><button type="button" data-open="candidate-'+str(i)+'" aria-label="View full profile for '+esc(gene)+'">Details</button></td></tr>')
        dialogs.append('<dialog id="candidate-'+str(i)+'" aria-labelledby="candidate-title-'+str(i)+'"><div class="dialog-head"><h3 id="candidate-title-'+str(i)+'">'+esc(gene)+' · '+esc(row['WBGene'])+'</h3><form method="dialog"><button>Close</button></form></div><div class="detail-body"><dl>'+''.join('<dt>'+esc(k.replace('_',' ').capitalize())+'</dt><dd>'+esc(value(v))+'</dd>' for k,v in row.items())+'</dl></div></dialog>')
    blocks.append('</tbody></table></div><p id="no-candidates" class="empty" hidden>No matching candidates. Try another identifier.</p></section>'+''.join(dialogs))
    blocks.append('<section id="sources"><h2>Evidence source status</h2>'+simple_table(['Source','Status','Records','Retrieved'],[(e['source'],e['status'],e.get('records'),str(e.get('retrieved_at',''))[:10]) for e in summary['evidence']]))
    for e in summary['evidence']:
        if e.get('reason'):blocks.append('<p><b>'+esc(e['source'])+':</b> '+esc(e['reason'])+'</p>')
    blocks.append('<details><summary>Source versions and retrieval details</summary>'+simple_table(['Source','Version','Retrieved'],[(e['source'],e.get('database_version'),e.get('retrieved_at')) for e in summary['evidence']])+'</details></section>')
    blocks.append('''</main><script>
    document.querySelectorAll('[data-open]').forEach(b=>b.addEventListener('click',()=>document.getElementById(b.dataset.open).showModal()));
    const search=document.getElementById('candidate-search');
    search.addEventListener('input',()=>{let count=0;document.querySelectorAll('tr[data-search]').forEach(r=>{r.hidden=!r.dataset.search.includes(search.value.trim().toLowerCase());if(!r.hidden)count++;});document.getElementById('candidate-count').textContent=count+' candidates shown';document.getElementById('no-candidates').hidden=count!==0;});
    </script></body></html>''')
    (out/'analysis_report.html').write_text('\n'.join(blocks))

if __name__ == '__main__':
    from pipeline import main
    main(default_stage=6)
