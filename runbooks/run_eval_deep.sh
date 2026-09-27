#!/bin/bash
# Wait for the training folds to finish -> ensemble test predictions -> metrics and GO/NO-GO.
# Idempotent: only folds that already have best.pt are processed.
set -uo pipefail
ROOT=/home/workspace/MD/MD/HDY/ADGRE1
P=${PREFIX:-$HOME/anaconda3/envs/deepscreen}
[ -x "$P/bin/python" ] || P=$ROOT/envs/deepscreen
OUT=$ROOT/deep_demo/output; LOGS=$ROOT/deep_demo/logs
ACC=${ACC:-cpu}            # use ACC=gpu when a GPU is visible
LOG=$LOGS/eval_deep.log; mkdir -p "$LOGS"
exec > >(tee -a "$LOG") 2>&1
echo "### run_eval_deep.sh $(date)  ACC=$ACC"

WAIT=${WAIT:-1}
if [ "$WAIT" = 1 ]; then
  while pgrep -f "chemprop train" >/dev/null 2>&1; do
    echo "  [$(date +%H:%M)] training still running: folds complete $(find $OUT -maxdepth 2 -name best.pt 2>/dev/null | wc -l)"
    sleep 60
  done
fi
FOLDS=$(find "$OUT" -maxdepth 2 -name best.pt -printf '%h\n' 2>/dev/null | grep -c . || true)
echo "folds available: $FOLDS"
[ "$FOLDS" -ge 1 ] || { echo "no finished fold yet; re-run later"; exit 1; }

echo "--- per-fold test predictions (ensemble input)"
i=0
for m in $(find "$OUT" -maxdepth 2 -name best.pt | sort); do
  i=$((i+1))
  "$P/bin/chemprop" predict --test-path "$ROOT/deep_demo/input/test.csv" \
    --model-path "$m" --preds-path "$OUT/test_preds_fold$i.csv" \
    --accelerator "$ACC" --devices 1 >/dev/null 2>&1 \
    && echo "  fold$i OK -> test_preds_fold$i.csv" || echo "  fold$i prediction failed"
done

echo "--- ensemble"
"$P/bin/python" - <<'PY'
import glob
import pandas as pd, numpy as np
fs = sorted(glob.glob('/home/workspace/MD/MD/HDY/ADGRE1/deep_demo/output/test_preds_fold*.csv'))
if not fs:
    raise SystemExit("no prediction files")
ds = [pd.read_csv(f) for f in fs]
pcol = list(ds[0].columns)[-1]
mx = np.stack([d[pcol].values for d in ds], 1)
out = ds[0].copy(); out['pred'] = mx.mean(1); out['pred_std'] = mx.std(1)
out.to_csv('/home/workspace/MD/MD/HDY/ADGRE1/deep_demo/output/test_preds.csv', index=False)
print("ensembled %d folds | prediction column=%s | mean across-fold std=%.3f"
      % (len(fs), pcol, float(mx.std(1).mean())))
PY

echo "--- metrics and verdict"
"$P/bin/python" "$ROOT/deep_demo/scripts/evaluate_deep_scores.py" \
  --pred "$OUT/test_preds.csv" --test "$ROOT/deep_demo/input/test.csv" \
  --out "$OUT/deep_metrics.json" || echo "evaluation failed (check the test_preds.csv format)"
echo "### done $(date)"
