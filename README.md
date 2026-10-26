# C. elegans physical PPI pipeline

A step-by-step protocol: **input genes → physical interactions → network → communities and hubs → GO → biological evidence → tables and figures**.

Each numbered Python file contains the analysis for that step. `pipeline.py` runs the steps in order; `common.py` provides shared downloads, identifier checks and file handling. There is no package tree to navigate, no screenshots to interpret and no pre-downloaded study documents required.

## 1. Install

Use Python 3.12 or later:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install NCBI BLAST+ separately and make sure `blastp` and `makeblastdb` are on your PATH. Network access is required for database retrieval. Missing BLAST software or an unavailable biological-evidence service is reported explicitly; other independent evidence sources continue.

## 2. Prepare your input

Create a tab-separated file called `input.tsv` with these headers:

```tsv
WBGene	UniProt	Orthogroup
WBGene00014086	Q23597	OG0004167
WBGene00011831	Q22562	OG0004728
```

These two rows illustrate the input format, not a sufficient analysis dataset. `WBGene` and `UniProt` are required; `Orthogroup` is optional and may be omitted or left empty. Multiple accessions for one WBGene are preserved, but the network has one node per gene. All counts are calculated from your file.

## 3. Run the complete pipeline

```bash
python pipeline.py --input input.tsv --out results
```

The normal workflow downloads machine-readable BioGRID, NCBI, GO, Bgee and UniProt data and queries WormBase. It does not use our historical study snapshots. Inputs, retrieval metadata, parameters, checksums and outputs are recorded in your local results folder.

If a biological-evidence database is unavailable, the terminal and `05_biological_evidence/evidence_status.tsv` identify the source and failure. A WormBase failure does not stop Bgee or BLAST, and results can still be generated. Missing records are never interpreted as proof of non-essentiality or absence of homology. Failure of a required upstream source such as BioGRID stops the dependent analysis with an error.

## Or run one step at a time

Use the **same output directory** for successive steps:

```bash
python 01_interactions.py --input input.tsv --out results
python 02_network.py --out results
python 03_topology.py --out results
python 04_go.py --out results
python 05_evidence.py --out results
python 06_results.py --out results
```

| Script | What it does | Main output inside your results folder |
|---|---|---|
| `01_interactions.py` | Download BioGRID and resolve identifiers; select physical interactions with both partners in C. elegans | `01_interactions/interactions.tsv` |
| `02_network.py` | Remove self-loops, collapse repeated pairs, retain evidence and extract one largest connected component | `02_network/lcc.graphml` |
| `03_topology.py` | Centrality, Louvain communities, robust hubs and publication-weighted sensitivity | `03_network_analysis/node_metrics.tsv` and `hubs.tsv` |
| `04_go.py` | GO enrichment of LCC communities | `04_functional_context/go_enrichment.tsv` |
| `05_evidence.py` | WormBase phenotypes, Bgee expression and human BLAST similarity for the hubs | `05_biological_evidence/evidence_status.tsv` and source-specific tables |
| `06_results.py` | Join evidence, draw the network/community/sensitivity figures and write an HTML report | `06_target_profiles/candidate_target_profiles.tsv` and `analysis_report.html` |

Open the HTML report in a browser. Figures are generated from result tables; they are not screenshots. The final step makes no database queries.

A completed step is not overwritten. Choose a new `--out` directory for a new analysis or after a failed step. Later steps verify the recorded upstream checksums. Sequential commands reload the saved parameters; use a new output folder to change them.

## Scientific settings

Defaults are in the single `parameters.yaml` file. For another settings file, pass `--parameters my_parameters.yaml` when starting the run.

- Physical interactions only; both interactor taxonomy IDs must be **6239 (Caenorhabditis elegans)**. The graph includes edges incident to at least one input gene, including external partners; it is not the complete organism interactome.
- Primary analysis is unweighted. Publication-weighted sensitivity uses the same graph and LCC.
- Louvain: resolution 1, 100 seeds; consensus threshold 0.5. Hubs: degree ≥95th percentile; retained as robust when present in ≥90% of 100 edge-retention replicates (90% of edges retained).
- GO: direct positive annotations, all LCC genes as background, one-sided Fisher tests, BH correction within each community; FDR <0.05. Complete results are retained.
- BLASTP: E-value ≤1e-5, BLOSUM62, no extra identity/coverage cutoff or reciprocal search. Similarity is not orthology or proof of host toxicity.
- WormBase preserves observed and explicitly not-observed phenotypes, evidence and references. It does not automatically classify genes as essential/non-essential.

Database versions and dates come from each run's retrieval records. New releases can change results. Candidate profiles are not validated therapeutic targets or a druggability score.

## Worked example: from input genes to results

A demonstration was previously run with **100 records sampled from 165 historical inputs**, corresponding to 76 unique WBGene IDs. Its WormBase evidence used the **5 August 2026 snapshot**. The example dataset is only a sample of data planned for full publication later.

[Download the example HTML report](https://github.com/LucyJimenez/ppi-target-prioritization/raw/refs/heads/main/example_report.html) and open it in your browser to explore the network, communities, GO enrichment and integrated hub evidence. The figures are embedded, so the report works offline. GitHub displays HTML source rather than rendering the report.

The default scripts do **not** use that snapshot. An explicit `--example-config /path/to/frozen-example.local.json` option exists only to validate a locally available frozen example. It requires local pinned resources and checksum-verified evidence matching the recalculated hubs; it never silently substitutes live data. Only the self-contained example HTML is included here; the complete example files and study results remain local.

This GitHub version contains the protocol, runnable scripts and one example HTML report. The previous expanded demonstration remains identifiable in release history; the current branch is the simplified workflow.
