# Stage 02 — Network construction

## Question

How is the physical interaction environment organized?

## Input

`results/01_interactions/interactions.tsv`: TSV requiring BioGRID interaction ID, experimental type, both organism IDs, WBGene_A/B and input_A/B. `candidate_annotations.tsv` and `input_mapping.tsv` preserve input identity and protein mapping.

## What this stage does

Validate Stage 01 constraints; preserve self-evidence; collapse unordered gene pairs; retain publication/system/record metadata; retain input isolates; identify components and one deterministic LCC.

## Output

Local directory: `results/02_network/` (a convenience link to an immutable run).

| File | Format / main columns | Meaning |
|---|---|---|
| nodes.tsv | WBGene, UniProt, Orthogroup, original_input, in_LCC, component | All input genes plus resolved partners |
| edges.tsv | WBGene_A/B, evidence_count, publication_count, experimental_system_count, support, distance | Unique non-self pairs |
| edge_evidence.tsv | Pair identifiers and BioGRID evidence | Detailed evidence behind each edge |
| connected_components.tsv | component, nodes; edges where available | Component inventory |
| network.graphml | GraphML nodes and edges | Full simple undirected graph |
| lcc_nodes.tsv / lcc_edges.tsv | Node/edge schemas above | The single largest connected component |
| lcc.graphml | GraphML | Exact LCC used for both analyses |
| network_summary.json | Nodes, edges, isolates and components | Structural counts |

## Run

```bash
python scripts/02_network.py --config config/config.yaml
```

## Example

Minimal illustrative key fields (synthetic, not a biological result):

```tsv
WBGene_A	WBGene_B	Experimental System Type	Organism ID Interactor A	Organism ID Interactor B
WBGene99000001	WBGene99000002	physical	6239	6239
```

## Notes

The example shows key fields only; use the complete Stage 01 table rather than a manually abbreviated input. No recursive partner expansion. Weights do not define LCC membership.
