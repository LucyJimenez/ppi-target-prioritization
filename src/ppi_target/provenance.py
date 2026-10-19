"""Immutable run bookkeeping; no scientific conclusions inferred from failures."""
from datetime import datetime, timezone
from pathlib import Path
import csv
import math
import hashlib
import json
import platform
import subprocess
import sys
import uuid
import yaml
from . import __version__


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    def clean(x):
        if isinstance(x,float) and not math.isfinite(x):return None
        if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [clean(v) for v in x]
        return x
    Path(path).write_text(json.dumps(clean(value), indent=2, ensure_ascii=False,allow_nan=False) + '\n')


def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='raise')
        w.writeheader()
        w.writerows(rows)


def load_config(path):
    path = Path(path).resolve()
    # Accept conventional YAML as well as its JSON subset.
    try:
        config = json.loads(path.read_text())
    except json.JSONDecodeError:
        config = yaml.safe_load(path.read_text())
    root = (path.parent / config['project_root']).resolve()
    return path, root, config


def new_run(root, config_path, config, run_id=None):
    run_id = run_id or datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8]
    if not run_id or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in run_id):
        raise ValueError('run_id must contain only letters, digits, underscore or hyphen')
    target = root / 'results' / 'runs' / run_id
    target.mkdir(parents=True, exist_ok=False)
    git = subprocess.run(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True, capture_output=True)
    manifest = {'run_id': run_id, 'started_at': now(), 'finished_at': None,
                'timezone': 'UTC', 'git_commit': git.stdout.strip() if git.returncode == 0 else None,
                'python_version': platform.python_version(), 'package_versions': {'ppi-target-prioritization': __version__},
                'configuration_file': str(config_path), 'configuration_sha256': sha256(config_path),
                'resolved_configuration': config, 'input_files': [], 'database_snapshots': [],
                'orthogroup_annotation': {'status': 'not_provided'}, 'stages': {}, 'outputs': [],
                'scripts_executed': [sys.argv[0]], 'commands': [sys.argv],
                'random_seeds': {}, 'warnings': ([] if git.returncode == 0 else ['No Git commit available; workspace is not a Git repository']),
                'errors': [], 'final_status': 'running'}
    write_json(target / 'run_manifest.json', manifest)
    return target, manifest
