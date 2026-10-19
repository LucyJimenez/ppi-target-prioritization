"""Offline figures; all inputs are saved local analysis outputs."""
import json,textwrap
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from .provenance import write_json,sha256


def generate(root,config,run,out):
    figdir=out/'figures';figdir.mkdir();md=out/'figure_manifests';md.mkdir()
    installed={f.name for f in font_manager.fontManager.ttflist}
    font=next((f for f in config['plotting']['fonts'] if f in installed),'DejaVu Sans')
    plt.rcParams.update({'font.family':font,'font.size':8,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':7,
                         'xtick.labelsize':7,'ytick.labelsize':7,'svg.fonttype':'none','pdf.fonttype':42,'savefig.dpi':config['plotting']['dpi']})
    nodes=pd.read_csv(run/'03_network_analysis/node_metrics.tsv',sep='\t',keep_default_na=False)
    g=nx.read_graphml(run/'02_network/lcc.graphml');cfg=config['plotting'];width=cfg['width_mm']/25.4
    cmap=plt.get_cmap('tab20');colors={c:cmap((int(c)-1)%20) for c in nodes.community.unique()}
    manifests=[]
    def save(fig,name,caption,inputs):
        fig.savefig(figdir/(name+'.svg'));fig.savefig(figdir/(name+'.pdf'));fig.savefig(figdir/(name+'.png'),dpi=cfg['dpi'])
        plt.close(fig)
        meta={'figure':name,'caption':caption,'font':font,'width_mm':cfg['width_mm'],'dpi':cfg['dpi'],'formats':['svg','pdf','png'],
              'inputs':[{'path':str(f.relative_to(root)),'sha256':sha256(f)} for f in inputs],'parameters':cfg,'external_queries':False}
        write_json(md/(name+'.json'),meta);manifests.append(meta)
    coords=nx.spring_layout(g,seed=cfg['seed'],k=cfg['spring_k'],iterations=cfg['iterations'],threshold=1e-5,weight=None,scale=1,method='force')
    pd.DataFrame([{'WBGene':n,'x':coords[n][0],'y':coords[n][1]} for n in sorted(coords)]).to_csv(figdir/'network_layout.tsv',sep='\t',index=False)
    ordered=nodes.set_index('WBGene').loc[list(g)];hub=ordered.robust_hub.astype(bool)
    fig,ax=plt.subplots(figsize=(width,7.5));fig.subplots_adjust(left=.02,right=.98,top=.90,bottom=.12)
    nx.draw_networkx_edges(g,coords,ax=ax,edge_color='#a1aeb8',alpha=.38,width=.35)
    nx.draw_networkx_nodes(g,coords,ax=ax,node_color=[colors[c] for c in ordered.community],node_size=np.where(hub,42,9),edgecolors=np.where(hub,'#101e2c','white'),linewidths=np.where(hub,.65,.15))
    top=nodes[nodes.robust_hub].sort_values(['degree','WBGene'],ascending=[False,True]).head(12)
    fig.canvas.draw();renderer=fig.canvas.get_renderer();placed=[]
    for r in top.itertuples():
        label=r.Gene_symbol or r.WBGene
        for radius in (12,22,34,48,64):
            found=False
            for angle in np.linspace(0,2*np.pi,12,endpoint=False):
                offset=(radius*np.cos(angle),radius*np.sin(angle))
                artist=ax.annotate(label,xy=coords[r.WBGene],xytext=offset,textcoords='offset points',fontsize=7,fontfamily=font,
                    ha='center',va='center',bbox={'facecolor':'white','alpha':.92,'edgecolor':'none','pad':.8},
                    arrowprops={'arrowstyle':'-','color':'#526b79','lw':.4})
                bounds=artist.get_window_extent(renderer).expanded(1.05,1.08)
                if not any(bounds.overlaps(b) for b in placed):
                    placed.append(bounds);found=True;break
                artist.remove()
            if found:break

    ax.set_title('A  Physical interaction landscape',loc='left',fontweight='bold',pad=16);ax.axis('off')
    handles=[plt.Line2D([0],[0],marker='o',color='none',markerfacecolor=colors[c],markeredgecolor='none',label=f'C{c}',markersize=5) for c in sorted(colors)]
    fig.legend(handles=handles,loc='lower center',ncol=9,frameon=False,bbox_to_anchor=(.5,.035))
    caption=f'Largest connected component: {len(g):,} WBGene nodes and {g.number_of_edges():,} physical edges. Colors show primary unweighted communities; outlined larger nodes are robust hubs. Labels identify {len(placed)} high-degree robust hubs. Layout is deterministic and unweighted.'
    save(fig,'figure_A_network',caption,[run/'02_network/lcc.graphml',run/'03_network_analysis/node_metrics.tsv'])
    go=pd.read_csv(run/'04_functional_context/community_go_summary.tsv',sep='\t',keep_default_na=False).sort_values('community')
    topgo=pd.read_csv(run/'04_functional_context/top3_BP.tsv',sep='\t')
    wrapped={int(r.community):[(textwrap.fill(str(t.GO_term),59),float(t.FDR)) for t in topgo[topgo.community==r.community].itertuples()] for r in go.itertuples()}
    heights=[max(.38,.13*sum(label.count('\n')+1 for label,_ in wrapped[int(c)])+.14) for c in go.community]
    y=np.cumsum(heights)-np.array(heights)/2;total=sum(heights)
    fig=plt.figure(figsize=(width,total+.9));gs=fig.add_gridspec(1,3,width_ratios=[1.05,.8,4.1],left=.075,right=.98,top=.91,bottom=.12,wspace=.18)
    ax1,ax2,ax3=[fig.add_subplot(gs[0,i]) for i in range(3)]
    ax1.barh(y,go.nodes,color=[colors[c] for c in go.community],height=np.minimum(heights,.44)*.75);ax1.set_yticks(y,[f'C{c}' for c in go.community]);ax1.set_xlabel('Nodes')
    ax2.barh(y,go.significant_BP_terms,color='#237e87',height=np.minimum(heights,.44)*.75);ax2.set_yticks([]);ax2.set_xlabel('GO BP terms')
    for ax in (ax1,ax2):ax.set_ylim(total,0);ax.spines[['top','right']].set_visible(False)
    ax3.set_xlim(0,1);ax3.set_ylim(total,0);ax3.axis('off')
    maxlog=max(2.,float(-np.log10(topgo.FDR.clip(lower=1e-300)).max())) if len(topgo) else 2.
    norm=matplotlib.colors.Normalize(vmin=0,vmax=maxlog);term_cmap=plt.get_cmap('viridis_r')
    for i,r in enumerate(go.itertuples()):
        labels=wrapped[int(r.community)]
        if not labels:ax3.text(0,y[i],'No significant BP term',fontsize=7,color='#6e7f8c',va='center',style='italic')
        else:
            line_count=sum(label.count('\n')+1 for label,_ in labels);at=y[i]-.13*line_count/2
            for label,fdr in labels:
                ax3.text(0,at,label,fontsize=7,va='top',linespacing=1.1,color=term_cmap(norm(-np.log10(max(fdr,1e-300)))))
                at+=.13*(label.count('\n')+1)
    fig.suptitle('B  Community size and functional context',x=.075,ha='left',fontsize=11,fontweight='bold')
    ax3.set_title('Top three BP terms (FDR, p, GO ID order)',fontsize=8,loc='left')
    cax=fig.add_axes([.66,.055,.25,.013]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=term_cmap),cax=cax,orientation='horizontal',label='Term color: -log10(FDR)')
    save(fig,'figure_B_community_GO','Community node counts, significant Biological Process terms (within-community BH FDR <0.05 across all tested aspects), and up to three ranked BP terms. No significant term is an explicit result, not missing data.',[run/'04_functional_context/community_go_summary.tsv',run/'04_functional_context/top3_BP.tsv'])
    weighted=pd.read_csv(run/'03_network_analysis/weighted_metrics.tsv',sep='\t');merged=nodes.merge(weighted[['WBGene','degree','strength','community','robust_hub']],on='WBGene',suffixes=('','_weighted'),validate='one_to_one')
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(width,3.2),layout='constrained')
    ax1.scatter(merged.degree,merged.strength_weighted,s=np.where(merged.robust_hub,22,6),alpha=.7,c=['#176e78' if x else '#a7bac3' for x in merged.robust_hub])
    ax1.set(xlabel='Unweighted degree',ylabel='Publication-weighted strength');ax1.set_title('C1  Same topology, different support',loc='left',fontsize=9)
    matrix=pd.crosstab(merged.community,merged.community_weighted);im=ax2.imshow(matrix,aspect='auto',cmap='Blues');ax2.set(xlabel='Weighted community',ylabel='Primary community');ax2.set_xticks(np.arange(len(matrix.columns))[::2],matrix.columns[::2]);ax2.set_yticks(np.arange(len(matrix.index))[::2],matrix.index[::2]);ax2.set_title('C2  Membership agreement',loc='left',fontsize=9);fig.colorbar(im,ax=ax2,shrink=.65,label='Genes')
    save(fig,'figure_C_sensitivity','Sensitivity uses identical nodes and edges; publication count is documentary support, not biochemical affinity. The contingency matrix compares membership independently of community label numbering.',[run/'03_network_analysis/node_metrics.tsv',run/'03_network_analysis/weighted_metrics.tsv'])
    write_json(md/'index.json',manifests)
