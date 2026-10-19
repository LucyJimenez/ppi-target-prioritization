# Stage 03 — Network analysis and prioritization

## Question

Which proteins occupy structurally important and robust positions?

## Input

`results/02_network/lcc.graphml`: the exact LCC, with support/distance edge attributes. `nodes.tsv`: WBGene, UniProt, Orthogroup, original_input and component. `network.graphml` supplies full-graph context.

## What this stage does

Compute primary unweighted centralities; consensus Louvain; perturbation-stable degree hubs; paired publication-weighted sensitivity; rank/community/hub comparisons; resolution/source-deletion/quantile diagnostics. Generate OG summaries automatically when annotated.

## Output

Local directory: `results/03_network_analysis/` (a convenience link to an immutable run).

| File | Format / main columns | Meaning |
|---|---|---|
| node_metrics.tsv | WBGene, UniProt, Orthogroup, community, centralities, ranks, hub flags | Primary node-level results |
| communities.tsv | community, nodes | Community sizes; memberships are in node_metrics.tsv |
| hubs.tsv | Node schema filtered to robust_hub=True | Prioritized candidates |
| weighted_metrics.tsv | Same node schema with weighted metrics | Same base nodes/edges/LCC |
| weighted_communities.tsv | WBGene, community | Weighted membership |
| hub_robustness.tsv | WBGene, hub_original, hub_frequency, robust_hub, hub_threshold | Auditable prioritization |
| og_summary.tsv | Orthogroup, WBGene_count, communities, robust_hubs | Only generated when applicable |
| sensitivity_summary.json | Graph, community, hub and agreement summaries | Primary-versus-weighted comparison |

## Run

```bash
python scripts/03_network_analysis.py --config config/config.yaml
```

## Example

Minimal illustrative key fields (synthetic, not a biological result):

```tsv
WBGene	UniProt	Orthogroup	community	degree	robust_hub
WBGene99000001	P12345	OG000001	1	3	True
```

## Notes

The example is synthetic. Production defaults: 100 Louvain seeds, resolution 1, coassignment 0.5; degree p95; 100 replicates retaining 90% edges; robustness ≥90%. Support counts distinct publications; distance=1/support. No composite target score.
