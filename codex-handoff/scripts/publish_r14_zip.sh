#!/usr/bin/env bash
set -euo pipefail
ARCHIVE=${1:?Usage: publish_r14_zip.sh /path/to/BlackHoles_R93R14_Complete_2026-10-12.zip}
REMOTE=${REMOTE:-origin}
BRANCH=${BRANCH:-codex/blackholes-r14-handoff-2026-10-12}
NAME=BlackHoles_R93R14_Complete_2026-10-12.zip
EXPECTED_SIZE=61043196
EXPECTED_SHA=2212ca6e2592de38121636a90efe5c92ae256968bf920e689e0acdc16b1711d1
REPO_PATH="codex-handoff/$NAME"
[[ -f "$ARCHIVE" ]] || { echo "Archive not found: $ARCHIVE" >&2; exit 2; }
[[ $(stat -c%s "$ARCHIVE") == "$EXPECTED_SIZE" ]] || { echo "Archive size mismatch" >&2; exit 3; }
[[ $(sha256sum "$ARCHIVE" | awk '{print $1}') == "$EXPECTED_SHA" ]] || { echo "Archive SHA-256 mismatch" >&2; exit 4; }
ROOT=$(git rev-parse --show-toplevel)
cd "$ROOT"
[[ $(git branch --show-current) == "$BRANCH" ]] || { echo "Checkout $BRANCH first" >&2; exit 5; }
cp -- "$ARCHIVE" "$REPO_PATH"
[[ $(sha256sum "$REPO_PATH" | awk '{print $1}') == "$EXPECTED_SHA" ]] || { echo "Copied SHA mismatch" >&2; exit 6; }
git add -- "$REPO_PATH"
git diff --cached --stat
git diff --cached --name-status
mapfile -t staged < <(git diff --cached --name-only)
[[ ${#staged[@]} -eq 1 && "${staged[0]}" == "$REPO_PATH" ]] || { echo "Refusing commit: unexpected staged paths" >&2; exit 7; }
git commit -m "data(research): add verified R14 handoff archive"
git push "$REMOTE" "HEAD:refs/heads/$BRANCH"
echo "Pushed $REPO_PATH"
echo "SHA256 $EXPECTED_SHA"
echo "Branch $BRANCH"
