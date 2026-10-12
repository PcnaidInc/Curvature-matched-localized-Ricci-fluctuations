param(
  [Parameter(Mandatory=$true)][string]$ArchivePath,
  [string]$Remote = "origin",
  [string]$Branch = "codex/blackholes-r14-handoff-2026-10-12"
)
$ErrorActionPreference = "Stop"
$ExpectedName = "BlackHoles_R93R14_Complete_2026-10-12.zip"
$ExpectedSize = 61043196
$ExpectedSha = "2212ca6e2592de38121636a90efe5c92ae256968bf920e689e0acdc16b1711d1"
$RepoPath = "codex-handoff/$ExpectedName"

if (-not (Test-Path -LiteralPath $ArchivePath -PathType Leaf)) { throw "Archive not found: $ArchivePath" }
$item = Get-Item -LiteralPath $ArchivePath
if ($item.Length -ne $ExpectedSize) { throw "Wrong archive size: $($item.Length); expected $ExpectedSize" }
$sha = (Get-FileHash -Algorithm SHA256 -LiteralPath $ArchivePath).Hash.ToLowerInvariant()
if ($sha -ne $ExpectedSha) { throw "SHA-256 mismatch: $sha; expected $ExpectedSha" }

$root = (git rev-parse --show-toplevel).Trim()
if (-not $root) { throw "Not inside a Git repository" }
Set-Location $root
$current = (git branch --show-current).Trim()
if ($current -ne $Branch) { throw "Wrong branch: '$current'. Checkout '$Branch' first." }

$dest = Join-Path $root ($RepoPath -replace '/', [IO.Path]::DirectorySeparatorChar)
Copy-Item -LiteralPath $ArchivePath -Destination $dest -Force
$copiedSha = (Get-FileHash -Algorithm SHA256 -LiteralPath $dest).Hash.ToLowerInvariant()
if ($copiedSha -ne $ExpectedSha) { throw "Copied archive hash mismatch" }

git add -- $RepoPath
if ($LASTEXITCODE -ne 0) { throw "git add failed" }
git diff --cached --stat
git diff --cached --name-status
$staged = (git diff --cached --name-only)
if (($staged | Measure-Object).Count -ne 1 -or $staged -ne $RepoPath) { throw "Refusing commit: staged paths are not exactly $RepoPath" }

git commit -m "data(research): add verified R14 handoff archive"
if ($LASTEXITCODE -ne 0) { throw "git commit failed" }
git push $Remote "HEAD:refs/heads/$Branch"
if ($LASTEXITCODE -ne 0) { throw "git push failed" }

Write-Host "Pushed $RepoPath"
Write-Host "SHA256 $ExpectedSha"
Write-Host "Branch $Branch"
