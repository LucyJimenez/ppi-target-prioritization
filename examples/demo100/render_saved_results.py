#!/usr/bin/env python3
"""Regenerate example figures and HTML offline from the distributed tables."""
import argparse,json,shutil,sys
from pathlib import Path
repo=Path(__file__).resolve().parents[2];sys.path.insert(0,str(repo/'src'))
from ppi_target.profiles import execute
from ppi_target.provenance import write_json

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',default=str(repo/'local/demo100-rendered'));a=p.parse_args()
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
    here=Path(__file__).resolve().parent
    for source in (here/'results').iterdir():
        if source.is_dir() and source.name!='06_target_profiles':shutil.copytree(source,out/source.name)
    cfg=json.loads((here/'display_config.json').read_text())
    manifest=json.loads((here/'results/run_manifest.json').read_text())
    execute(out,cfg,out,manifest)
    write_json(out/'run_manifest.json',manifest)
    print(out/'06_target_profiles/analysis_report.html')

if __name__=='__main__':main()
