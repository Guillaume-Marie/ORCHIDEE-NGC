#!/bin/bash
# Chaîne complète : extraction des sorties ORCHIDEE -> atlas.html autonome.
# 1) Éditer la CONFIGURATION en tête de extract_atlas_cube.py (SIM_DIR, RUN, OUT)
# 2) ./build_atlas.sh
set -e
cd "$(dirname "$0")"
python3 extract_atlas_cube.py
python3 assemble_atlas.py . atlas.html
echo "Ouvrir atlas.html dans un navigateur (double-clic suffit)."
