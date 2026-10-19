# Protein–Protein Interaction Network Analysis for Antiparasitic Target Prioritization

A modular and reproducible workflow that combines physical protein–protein interaction networks, network analysis, functional context, phenotype evidence, expression data, and host-protein sequence similarity to support antiparasitic target prioritization.

The endpoint is **candidate target profiles**, not drug discovery, demonstrated druggability or experimental validation. WBGene identifies a network node; UniProt is an associated protein; Orthogroup is optional annotation.

## Code public, data local

Pipeline code, documentation, example configuration, synthetic tests and the explicitly selected five-gene historical example below are intended for GitHub. Real inputs, database snapshots, FASTA, BLAST databases, results, figures, HTML reports and local configuration are excluded by `.gitignore`. No remote publishing command is part of this workflow. Do not use `git add -f` on excluded study files.

## Six stages

| Stage | Question | Principal local output |
|---|---|---|
| [01 Interactions](stages/01_interactions/README.md) | Which physical interactions are associated with input genes? | Validated interactions and evidence |
| [02 Network](stages/02_network/README.md) | How is the interaction environment organized? | One graph and one LCC |
| [03 Analysis](stages/03_network_analysis/README.md) | Which positions are structurally important and robust? | Metrics, communities, robust hubs and sensitivity |
| [04 Functional context](stages/04_functional_context/README.md) | Which functions characterize the LCC and communities? | GO coverage and community enrichment |
| [05 Biological evidence](stages/05_biological_evidence/README.md) | What evidence describes each prioritized candidate? | Independent WormBase, Bgee and BLAST results/statuses |
| [06 Target profiles](stages/06_target_profiles/README.md) | What available evidence supports each candidate? | Profiles, figures and English HTML report |

GO concerns the LCC and its communities, not only hubs. Candidate-specific evidence uses primary robust hubs. Unweighted analysis is primary; publication-weighted sensitivity uses the exact same nodes, edges and LCC. There is no arbitrary global target score.

## Installation

Python 3.12 and BLAST+ are required. The tested dependency versions are pinned in `requirements.txt`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
blastp -version
makeblastdb -version
```

If BLAST+ is missing on macOS with Homebrew: `brew install blast`. Alternatively use `conda env create -f environment.yml` and `conda activate ppi-target`. Run commands from the repository root. Installation does not download scientific data.

## Input and configuration

```bash
cp config/config.example.yaml config/config.yaml
mkdir -p data/raw
```

Place your own input at `data/raw/candidates.tsv`, or edit `input.path` in the local configuration:

```tsv
WBGene	UniProt	Orthogroup
WBGene99000001	P12345	OG000001
WBGene99000002	Q12345	OG000001
WBGene99000003	Q67890	
```

The example identifiers above are synthetic format illustrations. Required fields are WBGene and UniProt. Orthogroup may be absent or empty. A headerless two-column input is supported with `input.header: false`. One accession may map to multiple WBGene nodes; these nodes are not merged.

Orthogroup annotations are optional. When `data/raw/orthogroup_wbgene_mapping.tsv` is present, the pipeline automatically incorporates OG annotations and generates OG-aware outputs. When it is absent, all analyses that do not depend on Orthogroup information continue normally.

```tsv
Orthogroup	WBGene
OG000001	WBGene99000001
OG000001	WBGene99000002
```

Malformed OG files fail annotation validation with a QC report. Duplicates, unknown genes and multiple OG assignments are reported; multiple assignments are preserved explicitly. One OG can annotate several genes without collapsing nodes. Adding OG later does not change scientific cache keys: rerun normally to reuse compatible science and update annotations/profiles.

## First acquisition and execution

A fresh clone has no biological snapshots. After supplying input and reviewing parameters, explicitly acquire current sources:

```bash
python run_pipeline.py --config config/config.yaml --refresh
```

Every acquisition is versioned. A successful explicit refresh updates local snapshot pins in `config/snapshot_pins.local.json`; older snapshots remain available. A normal run reuses compatible snapshots and completed stages:

```bash
python run_pipeline.py --config config/config.yaml
```

The shipped example defers WormBase, reflecting the local study's deliberate pending integration. Set `wormbase.deferred: false` to include it in `--source all`; `--source wormbase` explicitly attempts WormBase even when deferral is configured.

## Individual stages and restart

```bash
python scripts/01_interactions.py --config config/config.yaml
python scripts/02_network.py --config config/config.yaml
python scripts/03_network_analysis.py --config config/config.yaml
python scripts/04_functional_context.py --config config/config.yaml
python scripts/05_biological_evidence.py --source all --config config/config.yaml
python scripts/06_target_profiles.py --config config/config.yaml

python run_pipeline.py --stage 3
python run_pipeline.py --from-stage 4
python run_pipeline.py --to-stage 3
python run_pipeline.py --resume RUN_ID --from-stage 4
```

Individual stages require compatible upstream results and fail clearly if those results are unavailable. `--resume` selects a previous run as the reuse source; it creates a new immutable run rather than overwriting the old one. Cache validation includes canonical WBGene/UniProt identity, scientific configuration, snapshot checksums, upstream cache identities, package versions and output checksums. Methods have an explicit version in the cache signature; update it when scientific semantics change. Source code copies are saved with each executed stage.

`--force` recomputes selected local stages without refreshing external sources. For Stage 05, forcing a source requests a new query/search for that source. `--refresh --source bgee --stage 5` explicitly acquires a new Bgee snapshot; other compatible successful evidence remains reusable. No successful snapshot is overwritten.

## Independent evidence failures

Biological evidence sources are queried independently. Failure or temporary unavailability of one source does not prevent the remaining evidence modules from running, and source-specific failures are explicitly recorded rather than interpreted as absence of biological evidence.

```bash
python scripts/05_biological_evidence.py --source wormbase --config config/config.yaml
python scripts/05_biological_evidence.py --source bgee --config config/config.yaml
python scripts/05_biological_evidence.py --source blast --config config/config.yaml
```

Each source is `SUCCESS`, `PARTIAL`, `FAILED` or `NOT_RUN`, with a reason. HTTP requests use 30-second connection and 120-second read timeouts, at most three attempts for transient failures, and bounded backoff. Permanent HTTP errors such as 403 stop retries. A failed source produces an error record, not fake zero-valued biological results. Stage 06 still runs with unavailable fields and explicit source status, even when all evidence sources fail.

`SUCCESS` with zero records means a successful empty response. `FAILED` with zero records means no biological conclusion is possible. Missing lethal phenotypes never establish non-essentiality. BLAST hits are sequence similarity, not automatic orthology; no hits do not establish host selectivity.

## Results and reports remain local

Each run writes `results/runs/RUN_ID/`, with six named stage directories, configuration and provenance in `run_manifest.json`, and source-code snapshots. `results/01_interactions/` through `results/06_target_profiles/` are convenience links to the latest corresponding outputs; immutable run directories are authoritative.

Stage 06 creates `candidate_target_profiles.tsv`, `evidence_completeness.tsv`, `analysis_summary.json`, figures in PDF/SVG/400-dpi PNG, figure manifests and `analysis_report.html`. The English HTML is self-contained for figures and works offline; linked TSV files remain alongside it. Full study products remain local; the only real-data documentation exception is the five-gene historical excerpt below.

## Small historical input/output example

This five-gene excerpt comes from the earlier `PPI_network_analysis` project. It illustrates input traceability and result presentation, including an available **WormBase REST snapshot retrieved on 2026-08-05** (release recorded as “live REST snapshot”; no numbered release established). It is not a new run of this pipeline and does not resolve the current run's pending WormBase integration. Only this explicitly selected excerpt is distributed; neither the full 165-row historical input nor the current study dataset is included.

Selection: the first five entries in the historical input mapping with both archived physical-network metrics and observed WormBase phenotype records. This is a demonstration subset, not a random or representative sample.

### Example input

Copy this tab-separated block into a local input file with `input.header: true`. The OG values are recovered orthogroup identifiers, not generated row labels.

```tsv
WBGene	UniProt	Orthogroup
WBGene00014086	Q23597	OG0004167
WBGene00011831	Q22562	OG0004728
WBGene00007934	G5ECE7	OG0004924
WBGene00011273	Q22021	OG0005398
WBGene00009739	O45525	OG0005773
```

### Example network output

Metrics below belong to the archived **physical-only** network, not the old project's mixed-interaction analysis. Degrees and community labels were calculated in that larger historical network, **not in a network of these five genes**. A new run using only this subset will not reproduce these numbers. Historical gene symbols and community numbering are preserved.

| WBGene | UniProt | Historical symbol | Degree | Historical community |
|---|---|---|---:|---:|
| WBGene00014086 | Q23597 | ZK809.3 | 37 | 1 |
| WBGene00011831 | Q22562 | T19B10.2 | 11 | 7 |
| WBGene00007934 | G5ECE7 | C34E7.4 | 8 | 6 |
| WBGene00011273 | Q22021 | R53.4 | 13 | 4 |
| WBGene00009739 | O45525 | F45H10.2 | 52 | 1 |

### Example WormBase evidence output

One observed annotation is shown per gene. For a compact illustration, embryonic-lethal records are selected first, then other lethal records, then phenotype ID and evidence identifier order. This selection is not a complete phenotype summary or an essentiality classification. Missing lethal records must never be interpreted as non-essentiality.

| WBGene | Observed phenotype | Phenotype ID | Evidence identifiers | Reference |
|---|---|---|---|---|
| WBGene00014086 | embryonic lethal | WBPhenotype:0000050 | WBPaper00005654;WBRNAi00009497 | WBPaper00005654 |
| WBGene00011831 | larval lethal | WBPhenotype:0000054 | WBPaper00029258;WBRNAi00074756 | WBPaper00029258 |
| WBGene00007934 | small | WBPhenotype:0000229 | WBPaper00006395;WBRNAi00024735 | WBPaper00006395 |
| WBGene00011273 | embryonic lethal | WBPhenotype:0000050 | WBPaper00004540;WBRNAi00000198 | WBPaper00004540 |
| WBGene00009739 | slow growth | WBPhenotype:0000031 | WBPaper00005654;WBRNAi00008806 | WBPaper00005654 |

All five displayed records are marked observed and carry RNAi/paper evidence in the archive. Archived evidence is shown for documentation only; it is not silently imported into current candidate profiles. [Example provenance](docs/historical_example_provenance.md) records the source files and checksums.

## Validation

```bash
python -m unittest discover -s tests -v
python run_pipeline.py --help
```

Tests cover identifier preservation, both-interactor taxon filtering, no overwriting, optional OG joins/propagation, malformed OG, one-OG/many-gene relations, five independent evidence-failure scenarios, failure versus zero evidence, coverage calculations, multiple-testing correction and synthetic integration with and without OG. No real analysis results are bundled with tests.

## Scientific limitations

The graph includes edges incident to input genes and is not the complete organism interactome. Publication weighting measures documentation, not binding affinity. Physical assay categories do not all imply direct binary binding. GO uses direct positive annotations and an LCC background, with annotation coverage reported. Expression counts include propagated conditions and are not independent experiments. Phenotype interpretation requires evidence and references; human similarity requires coverage-aware interpretation. Current unavailable resources are reported, never silently substituted.

Consult stage READMEs for schemas and exact outputs. MIT licensing applies to code; cite the software version/commit and each external biological resource separately.
