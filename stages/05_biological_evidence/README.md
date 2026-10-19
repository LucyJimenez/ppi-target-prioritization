# Stage 05 — Biological evidence

## Question

What phenotype, expression and human-similarity evidence exists for each prioritized candidate?

## Input

`results/03_network_analysis/hubs.tsv`: TSV with required WBGene/UniProt, optional Orthogroup and network metrics. Bgee gzip snapshot and metadata; NCBI gene_info for gene resolution; current UniProt query responses and human target FASTA; BLAST+ executables.

## What this stage does

Run WormBase, Bgee and BLAST independently. Save detailed source results, candidate summaries, raw responses, manifests and explicit failure records. Generate an exact candidate list and manual WormBase evidence template when integration is deferred.

## Output

Local directory: `results/05_biological_evidence/` (a convenience link to an immutable run).

| File | Format / main columns | Meaning |
|---|---|---|
| evidence_status.tsv | source, status, records, reason, database_version, retrieved_at, error_type, attempts | One source-specific state |
| evidence_summary.json | Stage status and source outcomes | Partial success is valid |
| wormbase/phenotypes.tsv | WBGene, phenotype ID/label, observed, evidence JSON, references | Generated only for usable phenotype responses |
| wormbase/candidate_summary.tsv | WBGene, observed/not-observed counts, record status | No automatic essentiality assertion |
| bgee/expression.tsv | Bgee gene, anatomy/stage, call, quality, rank/score and evidence fields | All retrieved candidate conditions |
| bgee/candidate_summary.tsv | WBGene, condition/anatomy/stage counts and record status | Descriptive expression context |
| blast/human_similarity.tsv | WBGene, UniProt, Orthogroup, query_id, human_hit_id/gene/UniProt, evalue, bit_score, percent_identity, alignment_length, query_coverage, subject_coverage | All qualifying query–human subject results |
| blast/candidate_summary.tsv | WBGene, query status, hits, best-hit summary | Sequence similarity, not orthology |
| */source_manifest.json and */error.json | Versions, parameters, dates, commands, checksums or errors | Permanent provenance |

## Run

```bash
python scripts/05_biological_evidence.py --source all --config config/config.yaml
```

## Example

Minimal illustrative key fields (synthetic, not a biological result):

```tsv
WBGene	UniProt	Orthogroup
WBGene99000001	P12345	
```

## Notes

WormBase phenotype classification IDs must be reviewed before essentiality categories are asserted. Missing phenotype evidence is not a negative phenotype. Bgee propagated conditions are not independent experiments. No arbitrary high-expression threshold; BLAST hits are not automatic orthologs.

### Independent sources and failure handling

```bash
python scripts/05_biological_evidence.py --source wormbase --config config/config.yaml
python scripts/05_biological_evidence.py --source bgee --config config/config.yaml
python scripts/05_biological_evidence.py --source blast --config config/config.yaml
python scripts/05_biological_evidence.py --source all --config config/config.yaml
```

`--source wormbase` explicitly attempts the service; default `all` respects configured deferral. Compatible successful sources are reused; add `--refresh` to acquire a new snapshot for the selected source. Each source catches its own exceptions, so a WormBase failure does not suppress Bgee or BLAST. `SUCCESS`, `PARTIAL`, `FAILED`, `NOT_RUN` are distinct controlled states with reasons. A successful empty response means zero matching biological records; a failed empty response supports no biological inference.

Connection/read timeouts are 30/120 seconds; transient HTTP/network failures receive at most three attempts with backoff, capped Retry-After delay. Permanent 403 is not retried indefinitely. Error records include type, HTTP code, attempts, timestamps, affected identifiers and response location. Source-specific output from another successful module is preserved. Stage 06 proceeds with explicit unavailable fields, even if all sources fail.

## Explicit historical WormBase snapshot

Set `wormbase.snapshot` to a directory, retrieval timestamp, database-version label and a filename-to-SHA256 mapping. The importer checks gene identity and file checksums, preserves observed and not-observed records separately, and marks hubs without a file as `not_in_snapshot` with missing counts. It never infers non-essentiality. The bundled [100-input demonstration](../../examples/demo100/README.md) uses this mode; it does not change the main study's deferred WormBase configuration. When snapshot mode is configured, it takes precedence over live REST requests, including `--refresh`: external databases may refresh, but WormBase remains frozen to the explicit snapshot.
