#!/usr/bin/env bash
# Publishes wiki/*.md (tracked in this repo) to the GitHub wiki
# (github.com/duceppemo/QIIME2_ITS.wiki.git, a separate git repo GitHub maintains
# alongside this one). The wiki/ folder here is the source of truth -- edit it, commit
# it to this repo like any other file, then run this script to push those pages live.
#
# Usage: scripts/sync_wiki.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WIKI_SRC="$REPO_ROOT/wiki"
WIKI_REMOTE="https://github.com/duceppemo/QIIME2_ITS.wiki.git"
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

git clone --quiet "$WIKI_REMOTE" "$WORKDIR"

# Mirror wiki/*.md onto the wiki clone, removing pages that no longer exist here.
find "$WORKDIR" -maxdepth 1 -name '*.md' -exec rm -f {} +
cp "$WIKI_SRC"/*.md "$WORKDIR"/

cd "$WORKDIR"
git add -A
if git diff --cached --quiet; then
    echo "Wiki already up to date, nothing to push."
    exit 0
fi

SRC_COMMIT="$(git -C "$REPO_ROOT" rev-parse --short HEAD)"
git commit --quiet -m "Sync from ${SRC_COMMIT} (duceppemo/QIIME2_ITS wiki/)"
git push --quiet origin master
echo "Wiki updated from $REPO_ROOT/wiki (source commit ${SRC_COMMIT})."
