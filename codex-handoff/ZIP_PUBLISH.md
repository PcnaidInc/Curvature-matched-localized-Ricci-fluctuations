# Publish the verified R14 ZIP to this handoff branch

The GitHub connector used to create this branch cannot ingest a large local binary by file reference. The R14 ZIP is below GitHub's 100 MiB per-file hard limit, so an authenticated local Git client can commit it normally.

Expected archive:

- `BlackHoles_R93R14_Complete_2026-10-12.zip`
- bytes: `61043196`
- SHA-256: `2212ca6e2592de38121636a90efe5c92ae256968bf920e689e0acdc16b1711d1`
- target: `codex-handoff/BlackHoles_R93R14_Complete_2026-10-12.zip`
- branch: `codex/blackholes-r14-handoff-2026-10-12`

Windows PowerShell:

```powershell
.\codex-handoff\scripts\publish_r14_zip.ps1 -ArchivePath "C:\path\BlackHoles_R93R14_Complete_2026-10-12.zip"
```

Linux / WSL:

```bash
bash codex-handoff/scripts/publish_r14_zip.sh /path/BlackHoles_R93R14_Complete_2026-10-12.zip
```

Both scripts refuse a wrong branch, byte size, SHA-256, or extra staged files. They commit and push only the verified ZIP. Because this repository is public, running either command publishes the archive publicly on that branch.
