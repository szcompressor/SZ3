#!/bin/bash
# Build the FZ docs site. Outputs to docs/site/html/.
#
# Run from anywhere:  bash docs/site/build.sh
#                or:  cd docs/site && bash build.sh
#
# Layout:
#   docs/
#     claude-skills/             installable Claude skills (parallel concern)
#     site/                      ALL website inputs + generated output
#       build.sh                  this script
#       Doxyfile                  Doxygen config
#       DoxygenLayout.xml         topbar / sidebar layout overrides
#       header.html, footer.html  Doxygen template overrides
#       mainpage.md               Doxygen home page (USE_MDFILE_AS_MAINPAGE)
#       pagefind.css              site-wide CSS (search bar, drawer, chat FAB)
#       modules.json              T2.2 catalog metadata (single source of truth)
#       module_page_template.html template for per-entry catalog pages
#       scripts/                  python build helpers
#       static/                   hand-written assets copied verbatim into html/
#       html/                     GENERATED output (gitignored), serve with
#                                   `cd docs/site/html && python3 -m http.server`
set -e

# Make all paths below resolve from the repo root regardless of cwd.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

# Pull SZ3 version from CMakeLists.txt so docs always show the real version.
# Doxyfile uses $(FZ_VERSION) which Doxygen expands from this env var.
export FZ_VERSION=$(grep -oE 'project\(SZ3 VERSION [0-9.]+\)' CMakeLists.txt | grep -oE '[0-9]+\.[0-9.]+')
echo "Building docs for FZ ${FZ_VERSION}"

# Clean stale top-level files so removed sources don't leave orphan HTML.
# Subdirs (pagefind/, modules/, assets/) are preserved.
find docs/site/html -maxdepth 1 -type f -delete 2>/dev/null || true

doxygen docs/site/Doxyfile

# Copy hand-written static assets (chat.js / chat.css / webgpu-test.html)
# into the generated tree. Keeps them under version control in docs/site/static/
# while making them servable from docs/site/html/.
mkdir -p docs/site/html/assets/chat
cp -R docs/site/static/. docs/site/html/

python3 docs/site/scripts/build_module_pages.py docs/site/modules.json docs/site/html/modules
python3 docs/site/scripts/build_chat_assets.py docs/site/modules.json docs/site/html/assets/chat/system-prompt.txt
npx -y pagefind --site docs/site/html
