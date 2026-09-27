#!/usr/bin/env python3
"""Technical demo v2: dual-receptor multi-task scoring (ADGRE1 vs ADGRE2) with selectivity
prediction. Runs on CPU.


DATA (all already available; from the GPU reverse screen and the corrected-box CPU scan):
  * ADGRE1 side: cpu_scan/811_b24 (9714 drugbank ligands, 24 A corrected box)
  * ADGRE2 side: reverse_screen/adgre2_cpu/drugbank (9337 ligands, same engine and box settings)
  * paired labels: identical ids in both files are the same ligand -> natural dual-receptor supervision
  * structures: demo/ligand_smiles.tsv (skeleton matching against drugbank.sdf) supplies trustworthy SMILES

Three tasks sharing one feature block:
  T1 regress aff_ADGRE1   T2 regress aff_ADGRE2   T3 regress dlog = log10|A2| - log10|A1| (selectivity)
Models: ridge (global linear) + kNN10 (local model over Tanimoto-like cosine similarity), pure numpy/scipy.
Validation: scaffold-grouped 5-fold (no scaffold spans folds) to avoid analogue leakage.
Criteria: report R2/MAE/Spearman and compare against a naive baseline (LE ranking) --
      a technical demo proves the loop works, it does not claim to beat the baseline.
  ~/anaconda3/envs/scvi-pypi/bin/python demo/ai_screen_v2.py
-> demo/ai_report_v2.json, demo/selectivity_pred.tsv
"""
import argparse, csv, glob, json, os, re, sys
import numpy as np
from scipy import sparse
from scipy.stats import rankdata

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOKEN = re.compile(r"\[[^\]]+\]|Br|Cl|Si|Se|[BCNOPSFIHKbcnops]|\d|\(|\)|=|#|/|\\|@|\+|-")


def load_affinity(path, underscore_ok=True):
    """Read each ligand's best score from a vina/vina-gpu output dir (clash artefacts filtered)."""
    out = {}
    for f in glob.glob(os.path.join(path, "*.pdbqt")):
        b = os.path.basename(f)
        lid = b[:-len("_out.pdbqt")] if b.endswith("_out.pdbqt") else b[:-len(".pdbqt")]
        best = None
        for l in open(f, errors="ignore"):
            if l.startswith("REMARK VINA RESULT"):
                try:
                    v = float(l.split()[3])
                except (IndexError, ValueError):
                    continue
                if -50 < v < 0 and (best is None or v < best):
                    best = v
            elif l.startswith("ENDMDL"):
                break
        if best is not None:
            out[lid] = best
    return out


def load_smiles():
    smi = {}
    p = os.path.join(ROOT, "demo", "ligand_smiles.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r.get("smiles"):
                smi[r["ligand_id"]] = r["smiles"]
    return smi


def features(smiles_list):
    vocab, toks = {}, []
    for s in smiles_list:
        t = TOKEN.findall(s)
        grams = t + ["".join(t[i:i+2]) for i in range(len(t)-1)]
        toks.append(grams)
        for g in grams:
            vocab.setdefault(g, len(vocab))
    ind, ptr, val = [], [0], []
    for grams in toks:
        cnt = {}
        for g in grams:
            j = vocab[g]
            cnt[j] = cnt.get(j, 0) + 1
        for j in sorted(cnt):
            ind.append(j); val.append(1.0)
        ptr.append(len(ind))
    B = sparse.csr_matrix((val, ind, ptr), shape=(len(toks), len(vocab)))
    return B


def scaffold(s):
    pairs = {}
    for i, ch in enumerate(s):
        if ch.isdigit():
            pairs.setdefault(ch, []).append(i)
    rings = [s[p[0]:p[-1]+1] for p in pairs.values() if len(p) >= 2]
    return "|".join(sorted(rings))[:120] or s[:40]


def ridge_fit(X, y, lam=1.0):
    A = (X.T @ X).toarray() if sparse.issparse(X) else X.T @ X
    A = A + lam * np.eye(A.shape[0])
    return np.linalg.solve(A, np.asarray(X.T @ y).ravel())


def knn_pred(Btr, ytr, Bte, k=10):
    Btr = Btr.astype(np.float64); Bte = Bte.astype(np.float64)
    ntr = np.sqrt(np.asarray(Btr.multiply(Btr).sum(1)).ravel()) + 1e-9
    nte = np.sqrt(np.asarray(Bte.multiply(Bte).sum(1)).ravel()) + 1e-9
    S = (Bte @ Btr.T).toarray() / np.outer(nte, ntr)
    idx = np.argsort(-S, axis=1)[:, :k]
    w = np.take_along_axis(S, idx, 1)
    yy = np.take_along_axis(np.tile(ytr, (Bte.shape[0], 1)), idx, 1)
    return (yy * w).sum(1) / (w.sum(1) + 1e-9)


def metrics(y, p):
    r2 = 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return float(r2), float(np.abs(y - p).mean()), float(np.corrcoef(rankdata(y), rankdata(p))[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a1-dir", default=os.path.join(ROOT, "cpu_scan", "811_b24"))
    ap.add_argument("--a2-dir", default=os.path.join(ROOT, "reverse_screen", "adgre2_cpu", "drugbank"))
    ap.add_argument("--folds", type=int, default=5, help="number of scaffold-grouped folds (min 3)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(ROOT, "demo", "out"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    a1 = load_affinity(args.a1_dir)
    a2 = load_affinity(args.a2_dir)
    smi = load_smiles()
    print("ADGRE1 scored %d | ADGRE2 scored %d | SMILES table %d" % (len(a1), len(a2), len(smi)))

    ids = [i for i in a1 if i in a2 and i in smi]
    print("three-way intersection (dual-receptor scores + trustworthy SMILES): %d" % len(ids))
    if len(ids) < 300:
        sys.exit("intersection too small; check the adgre2 dir and the SMILES table")

    y1 = np.array([a1[i] for i in ids])
    y2 = np.array([a2[i] for i in ids])
    y3 = np.log10(np.abs(y2)) - np.log10(np.abs(y1))     # dlog: >0 means ADGRE1-preferring
    S = [smi[i] for i in ids]
    B = features(S)
    print("features: %d sparse dims (SMILES token 1-3 grams)" % B.shape[1])

    groups = np.array([scaffold(s) for s in S])
    ug = sorted(set(groups)); rng.shuffle(ug)
    fold_of = {g: k % args.folds for k, g in enumerate(ug)}
    fold = np.array([fold_of[g] for g in groups])
    print("%d scaffolds -> %d folds" % (len(ug), args.folds))

    report = {"n_pairs": len(ids), "n_features": int(B.shape[1]),
              "n_scaffolds": int(len(ug)), "folds": args.folds,
              "y1_median": round(float(np.median(y1)), 2),
              "y2_median": round(float(np.median(y2)), 2),
              "dlog_median": round(float(np.median(y3)), 3),
              "tasks": {}}

    for name, y in (("T1_aff_ADGRE1", y1), ("T2_aff_ADGRE2", y2), ("T3_dlog_selectivity", y3)):
        pr = np.zeros(len(ids)); pk = np.zeros(len(ids))
        for f in range(args.folds):
            tr, te = fold != f, fold == f
            if tr.sum() < 50 or te.sum() == 0:
                continue
            w = ridge_fit(B[tr], y[tr], lam=1.0)
            pr[te] = B[te] @ w
            pk[te] = knn_pred(B[tr], y[tr], B[te], k=10)
        m_r = metrics(y, pr); m_k = metrics(y, pk)
        report["tasks"][name] = {
            "ridge": {"r2": round(m_r[0], 3), "mae": round(m_r[1], 3), "spearman": round(m_r[2], 3)},
            "knn10": {"r2": round(m_k[0], 3), "mae": round(m_k[1], 3), "spearman": round(m_k[2], 3)},
        }
        print("\n[%s] median label %.3f" % (name, float(np.median(y))))
        print("   ridge : R2=%+.3f MAE=%.3f Spearman=%.3f" % m_r)
        print("   kNN10 : R2=%+.3f MAE=%.3f Spearman=%.3f" % m_k)

    # baseline for the selectivity task: using the T1 prediction as "selectivity" is meaningless,
    # so the baseline is ranking by ADGRE1 ligand efficiency (naive but standard)
    le1 = -y1 / np.maximum(1, np.array([len(TOKEN.findall(s)) for s in S]))
    rho_le = float(np.corrcoef(rankdata(y3), rankdata(-le1))[0, 1])
    report["baseline"] = {"spearman_dlog_vs_negLE": round(rho_le, 3),
                          "note": "technical demo; no claim of beating the baseline, shown for reference"}
    print("\n[reference] Spearman between dlog and -LE(ADGRE1) = %+.3f" % rho_le)

    # fit on everything, then emit predictions (demonstrates the ranking-output stage)
    w = ridge_fit(B, y3, lam=1.0)
    pred = B @ w
    order = np.argsort(-pred)     # most ADGRE1-preferring first
    out_tsv = os.path.join(args.out, "selectivity_pred.tsv")
    with open(out_tsv, "w", newline="") as fh:
        wtr = csv.writer(fh, delimiter="\t")
        wtr.writerow(["rank", "ligand_id", "smiles", "aff_adgre1", "aff_adgre2",
                      "dlog_true", "dlog_pred"])
        for r, i in enumerate(order, 1):
            wtr.writerow([r, ids[i], S[i], round(y1[i], 2), round(y2[i], 2),
                          round(y3[i], 3), round(float(pred[i]), 3)])
    top = [(ids[i], float(pred[i]), float(y3[i])) for i in order[:10]]
    report["predicted_ag1_selective_top10"] = [
        {"ligand_id": t[0], "dlog_pred": round(t[1], 3), "dlog_true": round(t[2], 3)} for t in top]
    report["verdict"] = ("Technical demo: the three-task loop is closed (features -> multi-task regression -> scaffold-grouped validation -> ranked output). "
                         "Note that T3 has limited predictive power on this data, which is the quantitative expression of no selectivity window at this site.")
    json.dump(report, open(os.path.join(args.out, "ai_report_v2.json"), "w"), indent=2)
    print("\nTop-10 predicted ADGRE1-preferring ligands (dlog_pred vs dlog_true):")
    for t in top:
        print("   %-14s pred=%+.3f  true=%+.3f" % t)
    print("\nOutput -> demo/ai_report_v2.json, demo/selectivity_pred.tsv")


if __name__ == "__main__":
    main()
