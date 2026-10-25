"""Run all six steps, or invoke a numbered script for a single step."""
import argparse, importlib.util, json, subprocess, sys, traceback
from pathlib import Path
import yaml
from common import load, propagate, refresh, sha256, write_json, now
STEPS = {1: '01_interactions', 2: '02_network', 3: '03_topology', 4: '04_go', 5: '05_evidence', 6: '06_results'}
OUTPUTS = {1: '01_interactions', 2: '02_network', 3: '03_network_analysis', 4: '04_functional_context', 5: '05_biological_evidence', 6: '06_target_profiles'}
DEPENDENCIES = {1: [], 2: [1], 3: [2], 4: [2, 3], 5: [3], 6: [3, 4, 5]}

def module(number):
    spec = importlib.util.spec_from_file_location(STEPS[number], Path(__file__).with_name(STEPS[number] + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def main(default_stage=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, help='TSV with WBGene, UniProt and optional Orthogroup headers; required for step 1')
    p.add_argument('--out', type=Path, default=Path('results'), help='Local output directory; use the same directory for successive steps')
    p.add_argument('--parameters', type=Path, default=Path(__file__).with_name('parameters.yaml'))
    p.add_argument('--stage', type=int, choices=range(1, 7), default=default_stage)
    p.add_argument('--example-config', type=Path, help='Local frozen-example configuration ONLY; never used by the normal live pipeline')
    a = p.parse_args()
    root = a.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    mp = root / 'run_manifest.json'
    manifest = json.loads(mp.read_text()) if mp.exists() else {'started_at': now(), 'stages': {}, 'input_files': [], 'database_snapshots': [], 'errors': []}
    stages = [a.stage] if a.stage else list(STEPS)
    configfile = root / 'resolved_parameters.json'
    if configfile.exists():
        cfg = json.loads(configfile.read_text())
    else:
        cfg = yaml.safe_load(a.parameters.read_text())
        cfg['wormbase'] = {'deferred': False}
        cfg['annotation'] = {'optional': True, 'orthogroup_file': 'optional_orthogroup_mapping.tsv'}
        cfg['_source'] = 'all'
        cfg['_reuse_sources'] = []
        if a.example_config:
            cfg = json.loads(a.example_config.read_text())
            cfg['_example'] = True
            cfg['_source'] = 'all'
            cfg['_reuse_sources'] = []
        if not a.input:
            p.error('--input is required when starting a new output directory')
        cfg['input'] = {'path': str(a.input.resolve()), 'header': True}
    if a.input and str(a.input.resolve()) != cfg['input']['path']:
        p.error('This output directory belongs to another input; choose a new --out')
    if bool(a.example_config) != bool(cfg.get('_example')):
        p.error('Use --example-config consistently for a frozen-example output directory')
    inputfile = Path(cfg['input']['path'])
    digest = sha256(inputfile)
    if manifest.get('input_files') and manifest['input_files'][0]['sha256'] != digest:
        p.error('Input changed; choose a new --out')
    manifest['input_files'] = [{'path': str(inputfile), 'sha256': digest}]
    git = subprocess.run(['git', '-C', str(Path(__file__).parent), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    manifest.update(software_version='0.4.0', git_commit=git.stdout.strip() if git.returncode == 0 else None, command=sys.argv)
    try:
        _, canonical, qc = load(root, cfg, root / 'annotation')
        manifest['orthogroup_annotation'] = qc
        for stage in stages:
            key = f'{stage:02d}'
            dest = root / OUTPUTS[stage]
            if dest.exists():
                raise ValueError(f'Step {stage} already has output. Preserve it and use a new --out for a new run.')
            for required in DEPENDENCIES[stage]:
                if not (root / OUTPUTS[required]).exists():
                    raise ValueError(f'Run step {required} first with the same --out')
                for record in manifest.get('outputs', []):
                    if record['path'].startswith(OUTPUTS[required] + '/') and sha256(root / record['path']) != record['sha256']:
                        raise ValueError('Upstream output modified: ' + record['path'])
            print(f'Step {stage}: {STEPS[stage]}', flush=True)
            if not cfg.get('_example'):
                if stage in (1, 4):
                    cfg = refresh(root, cfg, {key}, 'all', root)
                if stage == 5:
                    cfg['_refresh_sources'] = ['bgee', 'blast']
            elif stage == 5:
                saved = cfg.get('example_evidence')
                if not saved:
                    raise ValueError('Frozen example requires example_evidence; no live query will be substituted')
                import shutil, pandas as pd
                source = Path(saved['path'])
                for name, expected in saved['sha256'].items():
                    if sha256(source / name) != expected:
                        raise ValueError('Example evidence checksum mismatch: ' + name)
                hubs = set(pd.read_csv(root / '03_network_analysis/hubs.tsv', sep='\t').WBGene)
                old = set(pd.read_csv(source / 'wormbase_manual_genes.tsv', sep='\t').WBGene)
                if hubs != old:
                    raise ValueError('Frozen evidence does not match the recalculated hub set')
                shutil.copytree(source, dest)
                summary = json.loads((dest / 'evidence_summary.json').read_text())
                manifest['stages'][key] = {'status': summary['status'], 'summary': summary, 'mode': 'frozen_example_only'}
            if not (cfg.get('_example') and stage == 5):
                module(stage).execute(root, cfg, root, manifest)
            propagate(root, canonical, qc)
            write_json(configfile, cfg)
            manifest['stages'][key]['finished_at'] = now()
            manifest['outputs'] = [{'path': str(f.relative_to(root)), 'sha256': sha256(f)} for d in OUTPUTS.values() for f in (root / d).rglob('*') if f.is_file()]
            write_json(mp, manifest)
        state = manifest['stages'].get('05', {}).get('status')
        manifest['status'] = 'PARTIAL' if state in ('FAILED', 'PARTIAL') else 'SUCCESS'
        print('Finished:', manifest['status'], root, flush=True)
    except Exception as exc:
        manifest['status'] = 'FAILED'
        manifest['errors'].append({'at': now(), 'message': str(exc)})
        print('ERROR:', exc, file=sys.stderr)
        raise
    finally:
        manifest['updated_at'] = now()
        write_json(mp, manifest)
if __name__ == '__main__':
    main()
