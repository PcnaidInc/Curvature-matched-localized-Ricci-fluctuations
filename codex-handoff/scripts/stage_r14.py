#!/usr/bin/env python3
"""Safely verify and stage the immutable R14 release under a gitignored directory.
No network access; no scientific calculation. Does not overwrite any prior stage.
"""
import argparse, hashlib, json, os, shutil, stat, tempfile, zipfile
from pathlib import Path, PurePosixPath
EXPECTED_NAME='BlackHoles_R93R14_Complete_2026-10-12.zip'
EXPECTED_SHA='2212ca6e2592de38121636a90efe5c92ae256968bf920e689e0acdc16b1711d1'
ROOT='BlackHoles_R93R14_2026-10-12'
REPO=Path(__file__).resolve().parents[2]
def file_sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()
def safe_name(n):
    p=PurePosixPath(n)
    if not n or n.startswith('/') or '\\' in n or '..' in p.parts or len(p.parts)<2 or p.parts[0]!=ROOT:
        raise ValueError(f'Unsafe/unexpected archive path: {n!r}')
    return p
def stage(archive:Path,dest:Path):
    archive=archive.expanduser().resolve()
    if not archive.is_file() or archive.name!=EXPECTED_NAME:raise FileNotFoundError('Supply exact R14 complete ZIP')
    if file_sha(archive)!=EXPECTED_SHA:raise ValueError('R14 archive SHA-256 mismatch')
    if dest.exists():raise FileExistsError(f'Existing stage is immutable: {dest}')
    if REPO not in dest.parents:raise ValueError('Destination must be inside repository .codex-local')
    names=set();meta=None
    with zipfile.ZipFile(archive) as z:
        if z.testzip():raise ValueError('ZIP CRC failed')
        if len(z.infolist())>3000 or sum(i.file_size for i in z.infolist())>350_000_000:raise ValueError('Archive size limit exceeded')
        for i in z.infolist():
            safe_name(i.filename)
            if i.filename in names:raise ValueError('Duplicate archive entry')
            names.add(i.filename)
            if (i.external_attr>>16)&0o170000==stat.S_IFLNK:raise ValueError('Symlink in ZIP')
        manifest_name=ROOT+'/MANIFEST.json'
        if manifest_name not in names:raise ValueError('Manifest missing')
        meta=json.loads(z.read(manifest_name))
        expected={ROOT+'/'+f['path'] for f in meta['files']}
        actual={n for n in names if not n.endswith('/') and n!=manifest_name}
        if expected!=actual:raise ValueError('Unlisted or missing payload(s)')
        for f in meta['files']:
            data=z.read(ROOT+'/'+f['path'])
            if len(data)!=f['bytes'] or hashlib.sha256(data).hexdigest()!=f['sha256']:
                raise ValueError('File hash mismatch: '+f['path'])
        dest.parent.mkdir(parents=True,exist_ok=True)
        temp=Path(tempfile.mkdtemp(prefix='r14-stage-',dir=dest.parent))
        try:
            z.extractall(temp)
            temp.rename(dest)
        except BaseException:
            if temp.exists():shutil.rmtree(temp)
            raise
    return {'destination':str(dest/ROOT),'archive_sha256':EXPECTED_SHA,'verified_payloads':len(meta['files']),'physics_reexecuted':False}
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--archive',required=True,type=Path);a.add_argument('--destination',type=Path,default=REPO/'.codex-local'/'r14');args=a.parse_args()
    print(json.dumps(stage(args.archive,args.destination.expanduser().resolve()),indent=2))
