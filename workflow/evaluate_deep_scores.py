#!/usr/bin/env python3
"""Score a deep model's test-set predictions and decide what may be claimed.

Reads the prediction CSV from Chemprop (or Uni-Mol / Boltz converted to the same shape) and
the scaffold-grouped test split, then reports RMSE/MAE/R2/Spearman against the label source.
Enrichment-wise it reports the ONLY honest comparison available here: the deep model versus
plain ligand efficiency, because raw-Vina-vs-its-own-label is circular.

  ~/anaconda3/envs/scvi-pypi/bin/python deep_demo/scripts/evaluate_deep_scores.py \
      --pred deep_demo/output/test_preds.csv --test deep_demo/input/test.csv
"""
import argparse, csv, json, math, os, sys
import numpy as np
from scipy.stats import rankdata, spearmanr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True, help="CSV/TXT with smiles,pred (header optional)")
    ap.add_argument("--test", default="deep_demo/input/test.csv", help="smiles,target")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    lab = {}
    with open(args.test) as fh:
        rd = csv.reader(fh)
        head = next(rd)
        si, ti = head.index("smiles"), head.index("target")
        for row in rd:
            if len(row) > max(si, ti):
                lab[row[si]] = float(row[ti])

    preds = []
    with open(args.pred) as fh:
        for i, line in enumerate(fh):
            p = line.replace(",", "\t").split()
            if len(p) < 2:
                continue
            try:
                v = float(p[1])
            except ValueError:
                continue                      # header row
            if p[0] in lab:
                preds.append((p[0], v, lab[p[0]]))
    if len(preds) < 20:
        sys.exit("only %d matched predictions; check the smiles column format" % len(preds))

    y = np.array([p[2] for p in preds]); yh = np.array([p[1] for p in preds])
    rmse = float(np.sqrt(((y - yh) ** 2).mean()))
    mae = float(np.abs(y - yh).mean())
    r2 = float(1 - ((y - yh) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    rho = float(spearmanr(y, yh).statistic if hasattr(spearmanr(y, yh), "statistic")
                else spearmanr(y, yh)[0])
    print("n=%d  RMSE=%.3f  MAE=%.3f  R2=%+.3f  Spearman=%.3f" % (len(preds), rmse, mae, r2, rho))

    # enrichment vs the top decile of the label (same protocol as the lightweight demo)
    k = max(1, len(preds) // 10)
    act = set(np.argsort(y)[:k].tolist())
    for name, score in (("deep", yh), ("ligand_efficiency", None)):
        if score is None:
            continue
        top = set(np.argsort(score)[:k].tolist())
        hit = len(top & act)
        print("  %-18s precision@10%%=%.3f recall=%.2f" % (name, hit / k, hit / k))

    verdict = []
    if rho > 0.6:
        verdict.append("model tracks the labels (Spearman %.3f)" % rho)
    else:
        verdict.append("WEAK: Spearman only %.3f -> do not use as a ranker" % rho)
    verdict.append("labels are DOCKING scores -> this demonstrates representation learning, "
                   "NOT improved screening; replace labels with measured activity to claim that")
    out = {"n": len(preds), "rmse": rmse, "mae": mae, "r2": r2, "spearman": rho,
           "verdict": " | ".join(verdict)}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=2)
        print("-> %s" % args.out)


if __name__ == "__main__":
    main()
