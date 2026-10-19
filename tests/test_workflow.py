import copy,csv,gzip,json,tempfile,unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import pandas as pd
from ppi_target.annotation import load,propagate
from ppi_target.sources import independent_sources
from ppi_target.profiles import integrate
from ppi_target.enrichment import bh
from ppi_target.human_homology import union_length
from ppi_target.provenance import load_config,sha256

class OptionalOGTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        (self.root/'input.tsv').write_text('WBGene\tUniProt\nWBGene99000001\tP12345\nWBGene99000002\tQ12345\nWBGene99000003\tQ67890\n')
        self.cfg={'input':{'path':'input.tsv','header':True},'annotation':{'orthogroup_file':'og.tsv'}}
    def test_absence_does_not_block_and_schema_is_unified(self):
        rows,c,q=load(self.root,self.cfg)
        self.assertEqual(list(c.columns),['WBGene','UniProt','Orthogroup']);self.assertEqual(q['orthogroup_annotation'],'not_provided');self.assertEqual(q['WBGene_without_OG'],3)
    def test_join_one_og_many_genes_and_propagation(self):
        (self.root/'og.tsv').write_text('Orthogroup\tWBGene\nOG000001\tWBGene99000001\nOG000001\tWBGene99000002\n')
        _,c,q=load(self.root,self.cfg);self.assertEqual(q['unique_OGs'],1);self.assertEqual(q['WBGene_with_OG'],2)
        d=self.root/'03_network_analysis';d.mkdir()
        pd.DataFrame({'WBGene':c.WBGene,'degree':[2,1,1],'community':[1,1,1],'robust_hub':[True,False,False]}).to_csv(d/'node_metrics.tsv',sep='\t',index=False)
        propagate(self.root,c,q);v=pd.read_csv(d/'node_metrics.tsv',sep='\t');self.assertEqual(v.degree.tolist(),[2,1,1]);self.assertEqual(v.Orthogroup[:2].tolist(),['OG000001']*2);self.assertTrue((d/'og_summary.tsv').exists())
    def test_malformed_file_fails_with_qc(self):
        (self.root/'og.tsv').write_text('Orthogroup\tWBGene\nBAD\tWBGene99000001\n')
        with self.assertRaises(ValueError):load(self.root,self.cfg,self.root/'qc')
        self.assertTrue((self.root/'qc/orthogroup_qc.json').exists())
    def test_duplicates_unknown_and_multiple_ogs_not_silently_lost(self):
        (self.root/'og.tsv').write_text('Orthogroup\tWBGene\nOG1\tWBGene99000001\nOG1\tWBGene99000001\nOG2\tWBGene99000001\nOG3\tWBGene99000999\n')
        _,c,q=load(self.root,self.cfg);self.assertEqual(c.iloc[0].Orthogroup,'OG1;OG2');self.assertEqual(len(q['issues']),3)
    def test_annotation_does_not_change_identity_checksum(self):
        _,_,a=load(self.root,self.cfg);(self.root/'og.tsv').write_text('Orthogroup\tWBGene\nOG1\tWBGene99000001\n');_,_,b=load(self.root,self.cfg)
        self.assertEqual(a['identity_checksum'],b['identity_checksum'])

class FailureIsolationTests(unittest.TestCase):
    def test_all_requested_scenarios(self):
        for failed in [{'wormbase'},{'bgee'},{'blast'},set(),{'wormbase','bgee','blast'}]:
            with self.subTest(failed=failed),tempfile.TemporaryDirectory() as td:
                called=[]
                def callback(name):
                    def run(folder):
                        called.append(name)
                        if name in failed:raise TimeoutError('simulated timeout')
                        return {'status':'SUCCESS','records':0,'reason':'successful empty biological response'}
                    return run
                names=['wormbase','bgee','blast'];result=independent_sources(names,{n:callback(n) for n in names},Path(td),['WBGene99000001'])
                self.assertEqual(called,names)
                self.assertEqual({r['source'] for r in result if r['status']=='FAILED'},failed)
                for n in failed:self.assertTrue((Path(td)/n/'error.json').exists())
    def test_failure_never_becomes_negative_essentiality(self):
        hubs=pd.DataFrame({'WBGene':['WBGene99000001'],'UniProt':['P12345'],'Orthogroup':['']})
        states=[{'source':'wormbase','status':'FAILED','reason':'timeout'},{'source':'bgee','status':'SUCCESS'},{'source':'blast','status':'SUCCESS'}]
        x=integrate(hubs,states,{'bgee':pd.DataFrame({'WBGene':['WBGene99000001'],'bgee_conditions':[0]})})
        self.assertEqual(x.wormbase_status.iloc[0],'FAILED');self.assertTrue(pd.isna(x.wormbase_essentiality.iloc[0]));self.assertTrue(pd.isna(x.wormbase_observed_records.iloc[0]));self.assertEqual(x.bgee_conditions.iloc[0],0)

class ScientificUtilities(unittest.TestCase):
    def test_union_coverage_not_double_counted(self):self.assertEqual(union_length([(1,10),(5,20),(30,35)]),26)
    def test_bh_includes_null_tests(self):
        self.assertAlmostEqual(bh([.01,1,1])[0],.03)
    def test_scientific_numeric_json_is_not_yaml_string(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'config.yaml';p.write_text('{"project_root":".","tolerance":1e-10}')
            self.assertIsInstance(load_config(p)[2]['tolerance'],float)

class FixtureIntegration(unittest.TestCase):
    def test_01_02_03_04_with_and_without_og(self):
        from ppi_target import biogrid,network,analysis,enrichment
        project=Path(__file__).resolve().parents[1]
        example=__import__('yaml').safe_load((project/'config/config.example.yaml').read_text())
        for has_og in (False,True):
            with self.subTest(OG=has_og),tempfile.TemporaryDirectory() as td:
                root=Path(td);cfg=copy.deepcopy(example);(root/'input.tsv').write_text('WBGene\tUniProt\nWBGene99000001\tP12345\nWBGene99000002\tQ12345\n')
                cfg['input']={'path':'input.tsv','header':True};cfg['annotation']={'orthogroup_file':'og.tsv','optional':True}
                if has_og:(root/'og.tsv').write_text('Orthogroup\tWBGene\nOG1\tWBGene99000001\nOG1\tWBGene99000002\n')
                with gzip.open(root/'gene.gz','wt') as f:
                    f.write('#tax_id\tGeneID\tdbXrefs\n')
                    for i in range(1,7):f.write(f'6239\t{i}\tWormBase:WBGene9900000{i}\n')
                from test_stage01 import record
                records=[record(i+1,a,b) for i,(a,b) in enumerate([(1,2),(1,3),(1,4),(2,3),(2,5),(2,6),(1,1)])]
                with (root/'bio.tsv').open('w') as f:
                    w=csv.DictWriter(f,fieldnames=records[0],delimiter='\t');w.writeheader();w.writerows(records)
                (root/'meta.json').write_text('{"database_version":"synthetic fixture","status":"SUCCESS"}')
                for key,name in [('tab3','bio.tsv'),('gene_info','gene.gz')]:cfg['biogrid'][key]={'path':name,'sha256':sha256(root/name),'retrieval_metadata':'meta.json'}
                cfg['network'].update(community_seeds=2,hub_replicates=4,resolution_sweep=[1.])
                (root/'go.obo').write_text('[Term]\nid: GO:9990001\nname: synthetic process\nnamespace: biological_process\n')
                with gzip.open(root/'go.gz','wt') as f:
                    for i in (1,2,3):f.write('\t'.join(['WB',f'WBGene9900000{i}','g','involved_in','GO:9990001','TEST:1','IEA','','P','','','protein','taxon:6239','20260101','TEST','',''])+'\n')
                for key,name in [('gaf','go.gz'),('obo','go.obo')]:cfg['go'][key]={'path':name,'sha256':sha256(root/name),'metadata':'meta.json'}
                run=root/'run';run.mkdir();manifest={'input_files':[],'database_snapshots':[],'stages':{},'warnings':[]}
                biogrid.execute(root,cfg,run,manifest);network.execute(root,cfg,run,manifest);analysis.execute(root,cfg,run,manifest);enrichment.execute(root,cfg,run,manifest)
                _,canon,qc=load(root,cfg);propagate(run,canon,qc)
                n=pd.read_csv(run/'02_network/nodes.tsv',sep='\t');self.assertEqual(len(n),6)
                metrics=pd.read_csv(run/'03_network_analysis/node_metrics.tsv',sep='\t');self.assertEqual(len(metrics),6)
                self.assertTrue((run/'04_functional_context/go_enrichment.tsv').exists())
                self.assertEqual(metrics.Orthogroup.notna().sum(),2 if has_og else 0)
                self.assertTrue(manifest['stages']['03']['summary']['same_nodes_edges_LCC'])

if __name__=='__main__':unittest.main()
