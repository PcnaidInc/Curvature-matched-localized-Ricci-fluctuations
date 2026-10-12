# Codex CLI continuation — BlackHoles (12 October 2026)

**Repository:** https://github.com/PcnaidInc/Curvature-matched-localized-Ricci-fluctuations
**Prepared branch:** `codex/blackholes-r14-handoff-2026-10-12` (do not merge automatically)
**Base:** `main` at `f282dae2d2de62acdd336283a4c8409e87c47e2f`
**Current scientific release:** `BlackHoles_R93R14_Complete_2026-10-12.zip`, SHA-256 `2212ca6e2592de38121636a90efe5c92ae256968bf920e689e0acdc16b1711d1`.
**Parent scientific work:** PCNBR-81. Last bounded execution: PCNBR-111. Scientific state: **unfinished**.

This repository's `main` branch contains an **older public manuscript/reproducibility set**. It is **not already a complete checkout of the latest R14 working campaign**. This branch adds a screened execution handoff and text-only source snapshots, not the 61 MB private R14 evidence archive. The complete source, large interval arrays, old checkpoints, updated 56-task tracker and failures are in the separate local R14 ZIP supplied with the private handoff download.

## Recommended startup

1. Clone the repository (or pull changes) and check out `codex/blackholes-r14-handoff-2026-10-12`. Read root `AGENTS.md`.
2. Place `BlackHoles_R93R14_Complete_2026-10-12.zip` from the private handoff download in a local accessible path. Do **not** stage it or any Jira material for Git.
3. From repository root, run `python codex-handoff/scripts/stage_r14.py --archive /path/to/BlackHoles_R93R14_Complete_2026-10-12.zip`. This verifies SHA-256, ZIP integrity, each listed payload, safe filenames, then extracts to gitignored `.codex-local/r14`.
4. Run the documented R14 `reproduce.py verify` and `certificate-check` first; these do not rerun the full field campaign. Record actual exit codes; do not claim tests pass without executing them.
5. Proceed with `NEXT_WORK.md` starting with the **comparison-width retarded response vs R12 local-order-nine trial** and historical geometry/field numerical-error audit. Save a first verified stage before any long numerical job.
6. Update the 56-task tracker **append-only**, create a linked Jira ledger entry if tools are connected, verify readback, and preserve separate report, science status, and unresolved errors.

**Read in this order:** `CURRENT_STATE.md`, `RESEARCH_INTEGRITY.md`, `NEXT_WORK.md`, `REPRODUCIBILITY.md`, `CODEX_PROMPT.md`. Archive-only inputs and program-specific examples are in `r14-source/`, which is public **source reference**, NOT the full release.
