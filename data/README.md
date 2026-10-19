# Data policy

Only small **synthetic format examples** in `example/` are public. They are not biological mappings and are not expected to match BioGRID. Executable integration fixtures are constructed in temporary directories by tests.

Put your own required WBGene/UniProt table in `raw/candidates.tsv`. An optional Orthogroup column may be empty. `raw/orthogroup_wbgene_mapping.tsv` is automatically discovered and joined by WBGene.

Real inputs, snapshots, sequences, local BLAST databases, logs, figures, reports and results stay local and are ignored by Git. Do not force-add these directories. The software license does not replace source-database licenses.
