import json,tempfile,unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import pandas as pd
from ppi_target.wormbase import from_snapshot
from ppi_target.provenance import sha256

class ArchivedPhenotypes(unittest.TestCase):
    def test_missing_is_not_negative_and_observed_is_separate(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'archive';source.mkdir();out=root/'output';out.mkdir()
            gene='WBGene99000001';file=source/(gene+'.json')
            ph={'phenotype':{'id':'WBPhenotype:0000050','label':'embryonic lethal'},'evidence':{'paper':'WBPaper00000001'}}
            file.write_text(json.dumps({'name':gene,'fields':{'phenotype':{'data':[ph]},'phenotype_not_observed':{'data':[ph]}}}))
            spec={'directory':'archive','files':{file.name:sha256(file)},'retrieved_at':'2026-08-05','database_version':'snapshot'}
            result=from_snapshot(root,spec,pd.DataFrame({'WBGene':[gene,'WBGene99000002']}),out)
            self.assertEqual(result['status'],'PARTIAL');self.assertEqual(result['records'],2)
            summary=pd.read_csv(out/'candidate_summary.tsv',sep='\t')
            self.assertEqual(summary.iloc[0].wormbase_observed_records,1)
            self.assertEqual(summary.iloc[0].wormbase_not_observed_records,1)
            self.assertTrue(pd.isna(summary.iloc[1].wormbase_observed_records))
            self.assertEqual(summary.iloc[1].wormbase_record_status,'not_in_snapshot')
    def test_snapshot_checksum_rejects_modified_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'archive').mkdir();(root/'out').mkdir();gene='WBGene99000001'
            (root/'archive'/f'{gene}.json').write_text('{}')
            spec={'directory':'archive','files':{gene+'.json':'wrong'},'retrieved_at':'2026-08-05','database_version':'snapshot'}
            with self.assertRaisesRegex(ValueError,'checksum mismatch'):
                from_snapshot(root,spec,pd.DataFrame({'WBGene':[gene]}),root/'out')

if __name__=='__main__':unittest.main()
