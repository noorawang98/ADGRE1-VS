#!/bin/bash
# Demo of the ligand-based part of the toolkit: recover trustworthy SMILES, build a dataset,
# run the geometry gate, then fit a surrogate model and test its enrichment.
#
# Needs: RDKit (for recover_smiles.py / build_dataset.py) and numpy+scipy (for the model).
#   RUN_RD="python -c 'import rdkit'"  RUN_ML=python  bash demos/run_demo.sh
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUT=${OUT:-$REPO/demo_out}
mkdir -p "$OUT" "$OUT/logs"
PY_RD=${PY_RD:-python}       # interpreter with RDKit
PY_ML=${PY_ML:-python}       # interpreter with numpy + scipy
step () { echo; echo "=== $1 ==="; }

step "environment"
bash "$REPO/environment/check_env.sh" 2>/dev/null | tail -12 || true
echo "RDKit interpreter : $PY_RD"
echo "model interpreter : $PY_ML"
echo "output directory  : $OUT"

step "1) recover trustworthy SMILES from an SDF with bond orders"
[ -s "$OUT/ligand_smiles.tsv" ] && echo "already present, skipping" || \
  $PY_RD "$REPO/demos/recover_smiles.py" --out "$OUT/ligand_smiles.tsv"

step "2) build the modelling dataset (scores + descriptors + ligand efficiency + alerts)"
$PY_RD "$REPO/demos/build_dataset.py" --out "$OUT/dataset.tsv" --ligand-table "$OUT/ligand_smiles.tsv"

step "3) geometry gate and ligand-efficiency ranking"
$PY_RD "$REPO/workflow/select_poses.py" --table "${AFFINITY_CSV:?set AFFINITY_CSV to an id,affinity table}" \
  --out "$OUT/selection.tsv"

step "4) surrogate model and enrichment test"
$PY_ML "$REPO/demos/ai_screen_surrogate.py" --data "$OUT/dataset.tsv"

step "5) multi-task demo (affinity for two receptors + selectivity)"
$PY_ML "$REPO/demos/ai_screen_multitask.py" --out "$OUT"

echo; echo "outputs in $OUT"
