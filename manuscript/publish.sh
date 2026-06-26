#!/usr/bin/env bash
#
# Build the manuscript with press and publish the static bundle to the
# `gh-pages` branch for GitHub Pages.
#
# The gh-pages branch holds ONLY the generated site. It is force-pushed as a
# single commit each time and never carries your source history — `main` is
# never touched. Run from anywhere:  ./manuscript/publish.sh
#
set -euo pipefail

MANUSCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIST="$MANUSCRIPT_DIR/dist"
BRANCH="gh-pages"
REMOTE_URL="$(git -C "$MANUSCRIPT_DIR" remote get-url origin)"

if ! command -v press >/dev/null 2>&1; then
  echo "error: the 'press' command is not on your PATH." >&2
  echo "       install it (see manuscript/README.md) and try again." >&2
  exit 1
fi

echo "==> Building bundle with press"
press build "$MANUSCRIPT_DIR"

echo "==> Publishing $DIST to '$BRANCH' on $REMOTE_URL"
cd "$DIST"
rm -rf .git
git init -q
git add -A
git -c user.name="press-deploy" -c user.email="deploy@local" \
  commit -q -m "Deploy manuscript $(date -u +%Y-%m-%dT%H:%M:%SZ)"
git branch -M "$BRANCH"
git push -f "$REMOTE_URL" "$BRANCH"
rm -rf .git

echo "==> Done. Live at https://blab.github.io/fitness-flux/ once Pages is enabled (see README)."
