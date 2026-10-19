# Stage 01 — Physical interaction discovery

## Question

Which physical protein–protein interactions are associated with the input candidate genes?

## Input

`data/raw/candidates.tsv`: TSV, required WBGene/UniProt, optional Orthogroup. `config/config.yaml`: YAML with input header policy and pinned BioGRID TAB3/NCBI gene_info sources. Optional `data/raw/orthogroup_wbgene_mapping.tsv`: TSV with Orthogroup/WBGene.

## What this stage does

Validate canonical identities, automatically attach optional OG metadata, inspect all physical organism/experimental-system combinations, retain physical 6239–6239 records incident to input genes, and preserve unresolved/cross-species evidence.

## Output

Local directory: `results/01_interactions/` (a convenience link to an immutable run).

| File | Format / main columns | Meaning |
|---|---|---|
| candidate_annotations.tsv | WBGene, UniProt, Orthogroup | Canonical three-column annotation; OG may be empty |
| interactions.tsv | Original BioGRID fields, WBGene_A/B, input_A/B, qc_status | Validated evidence for topology; self-records still present |
| interaction_evidence.tsv | Complete input-associated physical records | Preserves raw evidence multiplicity |
| unresolved_identifiers.tsv | Original fields and resolution/status | Explicit exclusions, not fabricated identities |
| cross_species_interactions.tsv | Original fields and both taxons | All physical cross-species QC |
| qc_summary.json | Counts and input mapping QC | Dynamic counts |
| source_manifest.json | Source/version/date/hash/path | Pinned retrieval provenance |

## Run

```bash
python scripts/01_interactions.py --config config/config.yaml
```

## Example

Minimal illustrative key fields (synthetic, not a biological result):

```tsv
WBGene	UniProt	Orthogroup
WBGene99000001	P12345	OG000001
WBGene99000002	Q12345	
```

## Notes

Both interactor taxa are checked explicitly. Sharing a UniProt accession never merges two WBGene nodes. No symbol-only fallback nodes are created. Refresh explicitly to obtain new database versions.
