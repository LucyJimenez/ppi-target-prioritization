# Stage 06 — Integrated target profiles

## Question

What available evidence profile supports each prioritized candidate?

## Input

`hubs.tsv` from Stage 03; `community_go_summary.tsv` from Stage 04; `evidence_status.tsv` and available successful `candidate_summary.tsv` files from each Stage 05 source. Optional OG metadata is automatically propagated.

## What this stage does

Join one row per WBGene with network robustness, community function and available candidate evidence. Preserve missing/failed/deferred states. Render offline publication figures and a self-contained English HTML methods/results report.

## Output

Local directory: `results/06_target_profiles/` (a convenience link to an immutable run).

| File | Format / main columns | Meaning |
|---|---|---|
| candidate_target_profiles.tsv | WBGene, UniProt, Orthogroup, community, metrics, robustness, GO context, source statuses/evidence | Interpretable candidate profiles; no global target score |
| evidence_completeness.tsv | WBGene and wormbase_status/bgee_status/blast_status | Evidence availability |
| figures/ | PDF, SVG, PNG and deterministic layout coordinates | Local publication outputs |
| figure_manifests/ | Inputs/checksums, layout/font and rendering parameters | Offline figure provenance |
| analysis_summary.json | Network, GO, evidence and OG summaries | Machine-readable results |
| analysis_report.html | Self-contained English HTML with embedded figures | Local paper-style methods, results and limitations |

## Run

```bash
python scripts/06_target_profiles.py --config config/config.yaml
```

## Example

Minimal illustrative key fields (synthetic, not a biological result):

```tsv
WBGene	wormbase_status	bgee_status	blast_status
WBGene99000001	FAILED	SUCCESS	SUCCESS
```

## Notes

FAILED or NOT_RUN WormBase produces NA/unavailable fields, never phenotype=absent or essentiality=false. The stage can finish even when all Stage 05 sources failed. The report and figures remain local and are ignored by Git.
