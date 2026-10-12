# Local reproducibility and command map

The public branch contains a lightweight, screened text-only R14 source snapshot at `codex-handoff/r14-source/`. It is **not** a stand-alone copy of all interval cell arrays and production evidence. For actual proof and reproducibility, obtain the original **61,043,196-byte** R14 release ZIP with SHA-256 `2212ca6e2592de38121636a90efe5c92ae256968bf920e689e0acdc16b1711d1` and run the staging script:

```shell
python codex-handoff/scripts/stage_r14.py --archive /path/to/BlackHoles_R93R14_Complete_2026-10-12.zip
cd .codex-local/r14/BlackHoles_R93R14_2026-10-12
python -m pip install -r requirements.txt
python reproduce.py verify
python reproduce.py certificate-check
python reproduce.py reduce-saved --out /tmp/blackholes-r14-reduced.json
```

On Windows PowerShell use a real destination outside the evidence directory, for example `$env:TEMP\blackholes-r14-reduced.json`. `verify` is a hash/integrity check, `certificate-check` checks **saved** rational/interval enclosures, and `reduce-saved` assembles cached output. **None re-evolves the field.** A new interval run requires `python reproduce.py preparation --out <new external directory> --time-div 64 --directions 8 --integral-div 8192 --dps 40`; it should not be called by default because it is new numerical work. Never send outputs into an existing source/evidence path.

The full final preparation setting, if actually needed, is `--time-div 256 --directions 32 --integral-div 32768 --dps 70`, evaluating 25,216 interval cells. It was already executed twice in the archived R14 evidence at different precisions and should not be repeated absent a predeclared goal.

The R14 archive contains `code/`, `vendor/`, `results/`, `verification/`, `inputs/`, `reports/METHODS.md`, `MANIFEST.json`, the immutable R13C release nested in `source_release/`, and the current research tracker. The file and SHA inventory in `DATA_MANIFEST.json` distinguishes screened Git source from the full private numerical bundle.

**Git discipline:** keep `.codex-local/` ignored. Use commits for code/docs and compact reviewed results; keep raw blocks local/private. Before proposing to merge any work, reproduce source hashes, inspect diffs for private data, and obtain scientific review/author authorization. Do not modify manuscript files under `manuscript/` for the next error-control batch.
