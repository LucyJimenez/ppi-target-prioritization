# Stage 04 — Functional context

## Question

Which biological functions characterize the LCC and its communities?

## Input

`results/03_network_analysis/node_metrics.tsv`: required WBGene/community; associated UniProt/Orthogroup metadata. Current pinned GAF 2.2 annotations and OBO ontology are configured in `config/config.yaml`.

## What this stage does

Validate active terms and ontology aspects; preserve direct positive annotations; exclude NOT; report LCC coverage; run community Fisher enrichment with all LCC genes as background and a fixed P/F/C test family; BH within community including zero-overlap terms.

## Output

Local directory: `results/04_functional_context/` (a convenience link to an immutable run).

| File | Format / main columns | Meaning |
|---|---|---|
| go_enrichment.tsv | community, GO_ID, GO_term, aspect, observed, p, FDR, fold_enrichment, supporting_WBGene | Complete test family, not just significant terms |
| community_go_summary.tsv | community, nodes, mapped/unmapped genes, significant BP terms, top3 | Community-level functional context |
| gene_annotation_coverage.tsv | WBGene, mapped, UniProt, Orthogroup | Whole-LCC mapping coverage |
| lcc_annotation_evidence.tsv | WBGene, GO_ID, qualifier, evidence_code, reference | Detailed direct annotation evidence |
| source_manifest.json | Annotation/ontology retrieval metadata | External provenance |
| summary.json | Coverage, tests, significance counts and semantic validation | Functional summary |

## Run

```bash
python scripts/04_functional_context.py --config config/config.yaml
```

## Example

Minimal illustrative key fields (synthetic, not a biological result):

```tsv
WBGene	community
WBGene99000001	1
WBGene99000002	1
```

## Notes

GO is not restricted to hubs. Whole-LCC coverage/annotation summaries are descriptive; no whole-LCC enrichment against a new species-wide background is invented. BP display uses FDR<0.05; ranking is FDR, p, GO ID. No propagation or redundancy pruning.
