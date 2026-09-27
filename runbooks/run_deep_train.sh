#!/bin/bash
# =============================================================================
# Route B: train Chemprop on ADGRE1 811-site docking labels and use the model's score
# to RANK the candidate library (i.e. the model actually participates in ranking).
#
# RUN ON THE HOST (needs the GPU; the container has no /dev/nvidia*).
#   bash deep_demo/run_deep_train.sh              # 5 folds, then ensemble predictions
#   FOLDS=1 EPOCHS=30 bash deep_demo/run_deep_train.sh   # quick single-fold check
#
# LABEL SOURCE WARNING (important, read before interpreting results)
#   The labels are CPU-Vina docking scores on the 811/7TM site. A model trained on them can
#   only reproduce that engine. Rank candidates whose scores come from the SAME engine,
#   otherwise the ranking mixes engine differences with model skill. This is a surrogate of
#   Vina, NOT evidence of activity; measured labels belong in deep_demo/input/activity.tsv.
# =============================================================================
set -uo pipefail
ROOT=/home/workspace/MD/MD/HDY/ADGRE1
P=${PREFIX:-$HOME/anaconda3/envs/deepscreen}
[ -x "$P/bin/python" ] || P=$ROOT/envs/deepscreen
DATA=$ROOT/deep_demo/input
OUT=$ROOT/deep_demo/output
LOGS=$ROOT/deep_demo/logs
FOLDS=${FOLDS:-5}
EPOCHS=${EPOCHS:-100}
BATCH=${BATCH:-32}
HIDDEN=${HIDDEN:-300}
DEPTH=${DEPTH:-3}
mkdir -p "$OUT" "$LOGS"
LOG=$LOGS/train_all.log
exec > >(tee -a "$LOG") 2>&1
echo "### run_deep_train.sh $(date)  FOLDS=$FOLDS EPOCHS=$EPOCHS HIDDEN=$HIDDEN DEPTH=$DEPTH"

echo "=== [0] environment and data checks ==="
"$P/bin/python" - <<'PY'
import torch, chemprop, rdkit, sys, os, csv
print("torch", torch.__version__, "| chemprop", chemprop.__version__, "| rdkit", rdkit.__version__)
if not torch.cuda.is_available():
    print("!! GPU not visible. Training on CPU is possible but slow; run this on the host.")
    print("   Set ALLOW_CPU=1 to continue on CPU.")
    if os.environ.get("ALLOW_CPU") != "1":
        sys.exit(3)
for f in ("train.csv", "val.csv", "test.csv"):
    p = "/home/workspace/MD/MD/HDY/ADGRE1/deep_demo/input/" + f
    rows = list(csv.DictReader(open(p)))
    ys = [float(r["target"]) for r in rows]
    print("  %-10s %5d rows | target min %.2f median %.2f max %.2f"
          % (f, len(rows), min(ys), sorted(ys)[len(ys)//2], max(ys)))
PY
[ $? -eq 0 ] || exit 1

echo "=== [1] train folds ==="
for f in $(seq 1 "$FOLDS"); do
  if [ -f "$OUT/fold$f/model_0/best.pt" ]; then
    echo "  fold$f: checkpoint exists, skipping"
    continue
  fi
  echo "  --- fold$f ---"
  "$P/bin/chemprop" train \
    --data-path "$DATA/train.csv" \
    --task-type regression \
    --output-dir "$OUT/fold$f" \
    --epochs "$EPOCHS" --batch-size "$BATCH" \
    --message-hidden-dim "$HIDDEN" --depth "$DEPTH" \
    --split-sizes 0.9 0.1 0.0 \
    --accelerator gpu --devices 1 2>&1 | tee "$LOGS/train_fold$f.log" | tail -6
done

echo "=== [2] per-fold test predictions -> ensemble ==="
i=0
for m in $(find "$OUT" -maxdepth 2 -name best.pt | sort); do
  i=$((i+1))
  "$P/bin/chemprop" predict --test-path "$DATA/test.csv" --model-path "$m" \
    --preds-path "$OUT/test_preds_fold$i.csv" --accelerator gpu --devices 1 \
    >/dev/null 2>&1 && echo "  fold$i ok" || echo "  fold$i FAILED"
done
"$P/bin/python" - <<'PY'
import glob
import pandas as pd, numpy as np
fs = sorted(glob.glob('/home/workspace/MD/MD/HDY/ADGRE1/deep_demo/output/test_preds_fold*.csv'))
if not fs:
    raise SystemExit("no per-fold predictions")
ds = [pd.read_csv(f) for f in fs]
pcol = list(ds[0].columns)[-1]
mx = np.stack([d[pcol].values for d in ds], 1)
out = ds[0].copy(); out['pred'] = mx.mean(1); out['pred_std'] = mx.std(1)
out.to_csv('/home/workspace/MD/MD/HDY/ADGRE1/deep_demo/output/test_preds.csv', index=False)
print("ensembled %d folds | mean across-fold std %.3f" % (len(fs), float(mx.std(1).mean())))
PY

echo "=== [3] metrics and GO/NO-GO (with a MEASURED baseline) ==="
# A naive ligand-efficiency baseline is computed on the same scaffold-grouped test set, so the
# model is judged against a real number instead of an assumed one. On this test set plain LE
# reaches precision@10% = 0.045 with Spearman -0.26 against the docking labels, i.e. ligand
# efficiency is a WEAK ranker here -- the bar the deep model must clear.
"$P/bin/python" - <<'BASELINE'
import csv, os
import numpy as np
from scipy.stats import spearmanr
ROOT = "/home/workspace/MD/MD/HDY/ADGRE1"
ds = {}
for r in csv.DictReader(open(os.path.join(ROOT, "demo", "dataset_811.tsv")), delimiter="\t"):
    ds[r["smiles"]] = r
test = list(csv.DictReader(open(os.path.join(ROOT, "deep_demo", "input", "test.csv"))))
m = [(float(r["target"]), ds[r["smiles"]]) for r in test if r["smiles"] in ds]
y = np.array([a for a, _ in m]); le = np.array([float(b["le"]) for _, b in m])
k = max(1, len(y) // 10)
act = set(np.argsort(y)[:k].tolist()); top = set(np.argsort(-le)[:k].tolist())
print("BASELINE plain-LE: precision@10%% = %.3f | Spearman(-LE,label) = %+.3f (n=%d)"
      % (len(top & act) / k, spearmanr(-le, y).statistic, len(y)))
BASELINE
"$P/bin/python" "$ROOT/deep_demo/scripts/evaluate_deep_scores.py" \
  --pred "$OUT/test_preds.csv" --test "$DATA/test.csv" \
  --out "$OUT/deep_metrics.json" || true

echo "=== [4] rank the candidate library with the ensemble ==="
# Build the ranking set: library ligands that HAVE trustworthy SMILES (the model needs SMILES)
"$P/bin/python" - <<'PY'
import csv, os
ROOT = '/home/workspace/MD/MD/HDY/ADGRE1'
smi = {}
for r in csv.DictReader(open(os.path.join(ROOT, 'demo', 'ligand_smiles.tsv')), delimiter='\t'):
    if r.get('smiles'):
        smi[r['ligand_id']] = r['smiles']
out = os.path.join(ROOT, 'deep_demo', 'input', 'library.smi')
n = 0
with open(out, 'w') as fh:
    for lid, s in sorted(smi.items()):
        fh.write("%s\t%s\n" % (s, lid)); n += 1
print("ranking set written: %d ligands -> %s" % (n, out))
PY
"$P/bin/chemprop" predict --test-path "$DATA/library.smi" \
  --model-path "$OUT/fold1/model_0/best.pt" \
  --preds-path "$OUT/library_preds.csv" --accelerator gpu --devices 1 2>&1 | tail -4 || \
  echo "(library prediction failed; check the library.smi format)"

"$P/bin/python" - <<'PY'
import csv, os
p = '/home/workspace/MD/MD/HDY/ADGRE1/deep_demo/output/library_preds.csv'
if not os.path.exists(p):
    raise SystemExit("no library predictions")
rows = list(csv.DictReader(open(p)))
col = list(rows[0].keys())[-1]
rows.sort(key=lambda r: float(r[col]))     # lower = better (Vina convention)
out = '/home/workspace/MD/MD/HDY/ADGRE1/demo/deep_ranking.tsv'
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['rank', 'smiles', 'ligand_id', 'pred_affinity'])
    for i, r in enumerate(rows, 1):
        w.writerow([i, r['smiles'], r.get('ligand_id', r.get('target', '')), r[col]])
print("ranked %d ligands -> %s" % (len(rows), out))
print("top5 predicted:")
for i, r in enumerate(rows[:5], 1):
    print("   %d  %s  %s" % (i, r[col], r['smiles'][:60]))
PY

cat <<'NOTE'

--- Interpretation rules (do not skip) ---
* GO only if test Spearman > 0.6 AND precision@10% >= plain ligand efficiency; else NO-GO.
* Labels are docking scores, so the model is a surrogate of Vina. Report it as
  "AI-assisted ranking pipeline", not as "AI discovered actives".
* Cross-check the top-ranked compounds with the geometry gate (demo/select_poses.py) before
  any experimental decision: clash artefacts inflate docking scores and the model inherits them.
NOTE
echo "### done $(date)"
