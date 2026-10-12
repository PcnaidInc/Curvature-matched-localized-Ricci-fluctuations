#!/usr/bin/env python3
"""Verify public branch's pinned source snapshot. No scientific computation."""
from pathlib import Path
import hashlib,json
D=Path(__file__).resolve().parents[1];S=D/'r14-source'
m=json.loads((S/'SOURCE_MANIFEST.json').read_text())
for path,data in m['source_files'].items():
    f=S/path
    assert f.is_file() and f.stat().st_size==data['bytes'],path
    assert hashlib.sha256(f.read_bytes()).hexdigest()==data['sha256'],path
for path in ['START_HERE.md','CURRENT_STATE.md','RESEARCH_INTEGRITY.md','NEXT_WORK.md','REPRODUCIBILITY.md','CODEX_PROMPT.md','DATA_MANIFEST.json']:
    assert (D/path).is_file(),path
print(json.dumps({'screened_source_files_checked':len(m['source_files']),'state':'source-only','physics_reexecuted':False}))
