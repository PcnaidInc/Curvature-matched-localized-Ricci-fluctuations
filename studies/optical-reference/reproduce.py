"""Usage: python reproduce.py verify | checks --out NEW_DIRECTORY
Verify hashes only, or rerun all bounded scientific checks in a new directory.
No prior SPV-1 archive is required for these reference benchmarks.
"""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('mode',choices=['verify','checks']);p.add_argument('--out',type=Path);a=p.parse_args()
if a.mode=='verify':
    manifest=json.loads((ROOT/'FILE_HASHES.json').read_text())
    for rel,expected in manifest.items():
        file=ROOT/rel
        if not file.is_file() or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
            raise SystemExit(f'FAIL: {rel}')
    print(f'PASS: {len(manifest)} hashes (no scientific calculation).')
else:
    if a.out is None: p.error('--out is required for checks')
    out=a.out.resolve()
    if out.exists() and any(out.iterdir()):p.error('Use a new or empty directory to preserve earlier results.')
    out.mkdir(parents=True,exist_ok=True)
    for script,folder in [('optical_benchmark.py','results'),('check_curved_optics.py','verification'),('check_normalization.py','verification'),('check_tail.py','verification')]:
        command=[sys.executable,str(ROOT/'code'/script),'--out',str(out/folder)]
        run=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        (out/(script+'.log')).write_text(run.stdout)
        if run.returncode:raise SystemExit(f'FAIL: {script}; inspect its saved log')
    print('PASS: bounded reference/derivation checks. No SPV-1 quantum-mode production.')
