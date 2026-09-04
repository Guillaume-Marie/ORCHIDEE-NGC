# Atlas ORCHIDEE — boîte à outils portable

Produit un **fichier HTML autonome** (`atlas.html`) : carte de l'Europe cliquable
des sorties ORCHIDEE (variables au choix, curseur d'année), chaque pixel ouvrant
une « planche pixel » (séries par classe d'âge, filtre par essence, âges porté /
pondéré aire / de flux, perturbations, lisière). Aucun serveur, aucun compte :
le fichier s'ouvre par double-clic, hors ligne, et se partage par mail ou clé USB
(données embarquées, quantifiées int16/int8 + gzip).

## Prérequis
- Python 3 avec `numpy` et `netCDF4` (pour l'extraction seulement).
- Un navigateur récent (Firefox/Chrome/Edge) pour la lecture — rien d'autre.

## Utilisation
1. Éditer la **CONFIGURATION** en tête de `extract_atlas_cube.py` :
   - `SIM_DIR` : dossier des sorties mensuelles (`.../SBG/Output/MO`)
   - `RUN` : préfixe des fichiers (`<RUN>_YYYY0101_YYYY1231_1M_stomate_history.nc`),
     sans underscore
   - `OUT` : où écrire les blobs intermédiaires (`.` par défaut)
2. `./build_atlas.sh` (ou les deux scripts python à la main)
3. Ouvrir `atlas.html`.

## Ce que la simulation doit fournir
Champ **requis** : `VEGET_MAX` par PFT (mensuel). Tous les autres sont
optionnels — un champ absent (p. ex. `AGE_STAND_AREA`, `AC12_F_REAL`,
`EDGE_LENGTH_DYN`, `DA_*` sur un trunk standard) donne simplement des panneaux
vides, sans faire échouer l'extraction. Champs exploités :
`DIA_INV, HEIGHT, AGE_STAND_BM, RDI, IND, TOTAL_M_c, AC12_F_REAL,
AGE_STAND_AREA, FOREST_MANAGED, CLASS0_VEG, FRAC_CLASS0, DA_FIRE, DA_STORM,
DA_PEST, DA_HARVEST, EDGE_LENGTH_DYN, Areas, CONTFRAC`.

## Adapter à une autre configuration PFT
Le mapping PFT → essence × classe d'âge est codé dans `extract_atlas_cube.py` :
- `TREE_BLOCKS` : liste de blocs `(début, fin)` **0-based** de PFT contigus, un
  bloc = une essence × ses classes d'âge dans l'ordre 1→4 ;
- `ESSENCES` : les noms affichés, dans le même ordre.
La valeur fournie correspond à la config 51 PFT / 4 classes (11 essences).
Pour 1 classe d'âge, le concept de « classes » perd son sens : l'atlas reste
utilisable mais les panneaux par classe n'auront qu'une courbe.

## Limites connues
- Grille régulière lat/lon attendue (testé à 1°) ; domaine Europe pour le fond
  de carte (frontières Natural Earth 50m embarquées, projection Lambert
  azimutale équivalente 52°N/10°E).
- Les taux `DA_*` sont interprétés comme des accumulateurs m²/pas de temps
  moyennés au mois (convention AED_FEEDBACK).
- Page limitée en pratique par le poids des données embarquées : ~150 ans × 1000
  pixels ≈ 15 Mo, au-delà réduire la fenêtre ou le domaine.

Origine : développé pour les simulations euage/euage2 (NextGenCarbon T4.3/T5.2),
Guillaume Marie / Science-Partners, avec l'assistance de Claude (Anthropic), 2026.
