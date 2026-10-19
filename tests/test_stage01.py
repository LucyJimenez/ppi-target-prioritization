import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from ppi_target.biogrid import inspect_records
from ppi_target.identifiers import read_mapping, mapping_qc, read_gene_info
from ppi_target.provenance import new_run


def record(i,a,b,ta='6239',tb='6239',kind='physical'):
    return {'#BioGRID Interaction ID':str(i), 'Entrez Gene Interactor A':str(a),
            'Entrez Gene Interactor B':str(b), 'Organism ID Interactor A':ta,
            'Organism ID Interactor B':tb, 'Experimental System Type':kind,
            'Experimental System':'Two-hybrid', 'Publication Source':'PUBMED:123'}


class IdentityTests(unittest.TestCase):
    def test_shared_accession_preserves_both_genes_and_first_row(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'map.tsv'
            p.write_text('WBGene99007208\tP12345\nWBGene99007206\tP12345\n')
            rows,errors=read_mapping(p)
            self.assertFalse(errors)
            qc=mapping_qc(rows)
            self.assertEqual(qc['unique_WBGene'],2)
            self.assertEqual(qc['unique_UniProt'],1)
            self.assertEqual(qc['duplicate_pair_extra_rows'],0)

    def test_both_taxa_and_stable_identity_required(self):
        g,h='WBGene99007208','WBGene99007206'
        m={'1':{g},'2':{h},'3':{'WBGene00000003','WBGene00000004'}}
        rows=[record(1,1,2), record(2,1,2,tb='9606'), record(3,1,2,ta='10090'),
              record(4,1,3),record(5,1,99),record(6,1,2,kind='genetic'),record(7,1,1)]
        r=inspect_records(rows,m,{g,h},6239)
        self.assertEqual(len(r['ready']),2) # self-record stays evidence for later cleaning
        self.assertEqual(len(r['rejected']),4)
        self.assertEqual(len(r['cross']),2)
        self.assertEqual(len(r['non_physical']),1)
        self.assertEqual(r['ready'][0]['WBGene_A'],g)
        self.assertEqual(r['ready'][0]['WBGene_B'],h)

    def test_conflicting_interaction_id_fails(self):
        with self.assertRaisesRegex(ValueError,'Conflicting'):
            inspect_records([record(1,1,2),record(1,1,3)],{'1':{'WBGene00000001'}},{'WBGene00000001'},6239)

    def test_gene_info_deduplicates_same_wbgene_xrefs_not_different_ids(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'gene.gz'
            with gzip.open(p,'wt') as f:
                f.write('#tax_id\tGeneID\tdbXrefs\n6239\t1\tWormBase:WBGene00000001|AllianceGenome:WB:WBGene00000001\n9606\t2\tWormBase:WBGene00000002\n')
            self.assertEqual(read_gene_info(p,6239),{'1':{'WBGene00000001'}})

    def test_run_does_not_overwrite(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); cfg=p/'config.json'; cfg.write_text('{}')
            new_run(p,cfg,{},'same-run')
            with self.assertRaises(FileExistsError):
                new_run(p,cfg,{},'same-run')

if __name__=='__main__':
    unittest.main()
