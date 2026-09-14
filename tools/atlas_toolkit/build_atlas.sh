#!/bin/bash
# Full chain: ORCHIDEE outputs -> standalone atlas.html.
#   1) edit atlas_config.py  (simulation, variables, layers, panels)
#   2) ./build_atlas.sh
set -e
cd "$(dirname "$0")"
python3 extract_atlas_cube.py
python3 make_borders.py atlas_meta.json borders_domain.json || true
python3 assemble_atlas.py . atlas.html
echo "Open atlas.html in a browser (a double-click is enough)."
