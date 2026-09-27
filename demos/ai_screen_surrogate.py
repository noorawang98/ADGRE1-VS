#!/usr/bin/env python3
"""AI screening demo on DrugBank/FDA-approved small molecules (no sklearn/rdkit needed).

Environment note: this box has no env carrying both RDKit and scikit-learn (adfr has RDKit
but no sklearn; scvi-pypi has sklearn but no RDKit; pip installs fail on the read-only
site-packages), so the model is implemented with numpy/scipy only and the chemistry comes
from the already-computed dataset (demo/dataset_811.tsv, built from SDF-matched SMILES).

Two complementary QSAR-style predictors, both standard in chemoinformatics:
  * global : ridge regression on a sparse chemistry + property block
  * local  : k-nearest-neighbour Tanimoto-style similarity over the same sparse block
Validation is scaffold-grouped (no scaffold in both train and test) and the demo reports
whether the models actually beat plain Vina / ligand-efficiency ranking.

  ~/anaconda3/envs/scvi-pypi/bin/python demo/ai_screen.py --data demo/dataset_811.tsv
-> demo/ai_report.json, demo/shortlist.tsv
"""
import argparse, csv, json, os, re, sys
import numpy as np
from scipy import sparse
from scipy.stats import rankdata

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOKEN = re.compile(r"\[[^\]]+\]|Br|Cl|Si|Se|Na|Li|Mg|Ca|Zn|Fe|Cu|Mn|Sn|Ag|Au|Hg|Pb|Pt|"
                   r"[BCNOPSFIHKbcnops]|\d|\(|\)|=|#|/|\\|@|\+|%-?\d*")


def tokenize(smiles):
    """SMILES -> token multiset. A count-based token block is the bond-order-free analogue of
    a Morgan fingerprint: it needs no ring perception, so it works on SDF-recovered SMILES
    even when the original PDBQT was bond-order blind."""
    return TOKEN.findall(smiles)


def features(rows):
    """bag-of-tokens (substring n-grams of the token stream) + 6 physchem descriptors"""
    vocab, toks = {}, []
    for r in rows:
        t = tokenize(r["smiles"])
        grams = t + ["".join(t[i:i + 2]) for i in range(len(t) - 1)]
        grams += ["".join(t[i:i + 3]) for i in range(len(t) - 2)]
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
    n = len(rows)
    B = sparse.csr_matrix((val, ind, ptr), shape=(n, len(vocab)))
    P = np.array([[float(r["mw"]), float(r["clogp"]), float(r["tpsa"]), float(r["hbd"]),
                   float(r["hba"]), float(r["rotb"])] for r in rows])
    P = (P - P.mean(0)) / (P.std(0) + 1e-9)
    return sparse.hstack([B, sparse.csr_matrix(P)]).tocsr(), B


def ridge_fit(X, y, lam=1.0):
    A = (X.T @ X).toarray() if sparse.issparse(X) else X.T @ X
    A = A + lam * np.eye(A.shape[0])
    b = np.asarray(X.T @ y).ravel()
    return np.linalg.solve(A, b)


def scaffold_key(smiles):
    """Bemis-Murcko-style ring-system key without RDKit: keep ring-closure pairs and the
    atoms between them, drop everything else. Used only to group folds so the same scaffold
    cannot appear in both train and test."""
    rings, stack = [], []
    for i, ch in enumerate(smiles):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack:
                stack.pop()
    pairs = {}
    for i, ch in enumerate(smiles):
        if ch.isdigit():
            pairs.setdefault(ch, []).append(i)
    for k, pos in pairs.items():
        if len(pos) >= 2:
            rings.append(smiles[pos[0]:pos[-1] + 1])
    return "|".join(sorted(rings))[:120] or smiles[:40]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(ROOT, "demo", "dataset_811.tsv"))
    ap.add_argument("--top", type=int, default=100)
    ap.add_argument("--active-frac", type=float, default=0.05)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    rows = [r for r in csv.DictReader(open(args.data), delimiter="\t") if r.get("smiles")]
    for r in rows:
        r["affinity"] = float(r["affinity"]); r["hac"] = int(r["hac"]); r["le"] = float(r["le"])
    print("dataset: %d ligands with trustworthy SMILES" % len(rows))
    if len(rows) < 200:
        sys.exit("need >=200 rows; let the matcher and the 24 A rescan progress further")

    X, B = features(rows)
    y = np.array([r["affinity"] for r in rows])
    groups = np.array([scaffold_key(r["smiles"]) for r in rows])
    print("features: %d (%d token grams) | scaffolds: %d" % (X.shape[1], B.shape[1], len(set(groups))))

    # scaffold-grouped folds
    ug = list(set(groups)); rng.shuffle(ug)
    fold_of = {g: i % args.folds for i, g in enumerate(ug)}
    fold = np.array([fold_of[g] for g in groups])

    pred_ridge = np.zeros(len(rows)); pred_knn = np.zeros(len(rows))
    for f in range(args.folds):
        tr, te = fold != f, fold == f
        if tr.sum() < 50 or te.sum() == 0:
            continue
        w = ridge_fit(X[tr], y[tr], lam=1.0)
        pred_ridge[te] = X[te] @ w
        # local model: similarity-weighted mean of the k most similar training ligands
        Btr = B[tr].astype(np.float64); Bte = B[te].astype(np.float64)
        nrm_tr = np.sqrt(np.asarray(Btr.multiply(Btr).sum(1)).ravel()) + 1e-9
        nrm_te = np.sqrt(np.asarray(Bte.multiply(Bte).sum(1)).ravel()) + 1e-9
        S = (Bte @ Btr.T).toarray() / np.outer(nrm_te, nrm_tr)     # cosine ~ Tanimoto proxy
        k = 10
        idx = np.argsort(-S, axis=1)[:, :k]
        wts = np.take_along_axis(S, idx, 1)
        pred_knn[te] = (np.take_along_axis(np.tile(y[tr], (te.sum(), 1)), idx, 1) * wts).sum(1) / (wts.sum(1) + 1e-9)

    def metrics(p):
        ss = 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        return float(ss), float(np.abs(y - p).mean()), float(np.corrcoef(rankdata(y), rankdata(p))[0, 1])

    r2r, maer, spr = metrics(pred_ridge)
    r2k, maek, spk = metrics(pred_knn)
    print("\n[1] affinity prediction (scaffold-grouped %d-fold)" % args.folds)
    print("    ridge : R2=%+.3f MAE=%.3f Spearman=%.3f" % (r2r, maer, spr))
    print("    kNN10 : R2=%+.3f MAE=%.3f Spearman=%.3f" % (r2k, maek, spk))

    # ---- enrichment over size-matched decoys ----
    order = np.argsort(y)
    n_act = max(10, int(len(y) * args.active_frac))
    act = set(order[:n_act].tolist())
    lo, hi = np.percentile([rows[i]["hac"] for i in act], [5, 95])
    pool = [i for i in range(len(y)) if i not in act and lo <= rows[i]["hac"] <= hi]
    dec = set(rng.choice(pool, size=min(len(pool), 4 * n_act), replace=False).tolist())
    sel = np.array(sorted(act | dec))
    ybin = np.array([1 if i in act else 0 for i in sel])

    def auc(score):
        # AUC of "actives score better (lower) than decoys" = AUC of the negated score
        r = rankdata(-np.asarray(score))
        n1, n0 = ybin.sum(), len(ybin) - ybin.sum()
        return float((r[ybin == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

    def ef1(score):
        # scores follow the Vina convention (more negative = better), so take the smallest
        k = max(1, int(0.01 * len(sel)))
        top = sel[np.argsort(score)[:k]]
        return float((sum(1 for i in top if i in act) / k) / (n1_frac := (ybin.sum() / len(ybin))))

    # every score here is in the same orientation: LOWER = better
    scores = {"ridge": pred_ridge[sel], "knn10": pred_knn[sel],
              "raw_vina": y[sel], "ligand_efficiency": -np.array([rows[i]["le"] for i in sel])}
    print("\n[2] enrichment (actives = best %.0f%% by Vina, size-matched decoys)" % (100 * args.active_frac))
    res = {}
    for k, v in scores.items():
        a, e = auc(v), ef1(v)
        res[k] = {"auc": round(a, 3), "ef1pct": round(e, 2)}
        print("    %-18s AUC=%.3f EF1%%=%.2f" % (k, a, e))

    # ---- transfer ranking + shortlist ----
    w_all = ridge_fit(X, y, lam=1.0)
    pred_all = X @ w_all
    win = {}
    thr = os.path.join(ROOT, "filters", "approved_drug_thresholds.json")
    if os.path.exists(thr):
        win = json.load(open(thr))["gates"]["approved_like_p5_p95"]
    out = []
    for i, r in enumerate(rows):
        gate = "pass"
        for kk in ("mw", "clogp", "tpsa", "hbd", "hba", "rotb"):
            if win and kk in win and not (win[kk][0] <= float(r[kk]) <= win[kk][1]):
                gate = "fail"; break
        if gate == "pass" and r.get("alerts"):
            gate = "flag"
        # ai_score: lower = better; combines predicted affinity, efficiency and an alert penalty
        out.append({**r, "pred_affinity": round(float(pred_all[i]), 3), "gate": gate,
                    "ai_score": round(float(pred_all[i]) + 2.0 * (0.25 - r["le"]) +
                                      (0.5 if r.get("alerts") else 0.0), 3)})
    out.sort(key=lambda r: r["ai_score"])
    with open(os.path.join(ROOT, "demo", "shortlist.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader()
        for r in out:
            w.writerow(r)

    # IMPORTANT caveat: the actives ARE defined by raw Vina score, so raw_vina's AUC is
    # circular and meaningless as a benchmark. The only honest benchmark is ligand efficiency
    # (size-decoupled, not derived from the label definition).
    best_ml_auc = max(res["ridge"]["auc"], res["knn10"]["auc"])
    best_ml_ef = max(res["ridge"]["ef1pct"], res["knn10"]["ef1pct"])
    le_auc, le_ef = res["ligand_efficiency"]["auc"], res["ligand_efficiency"]["ef1pct"]
    if best_ml_ef > le_ef + 0.5 and best_ml_auc > le_auc + 0.02:
        verdict = ("ML adds signal: AUC %.3f vs LE %.3f, EF1%% %.2f vs %.2f"
                   % (best_ml_auc, le_auc, best_ml_ef, le_ef))
    else:
        verdict = ("NO-GO as a ranker: ML (AUC %.3f / EF1%% %.2f) does not beat plain ligand "
                   "efficiency (AUC %.3f / EF1%% %.2f) on this dataset. ML is still useful as a "
                   "QSPR surrogate (R2 %.3f) but must not be presented as improving enrichment. "
                   "raw_vina AUC=1.0 is circular and is NOT a benchmark."
                   % (best_ml_auc, best_ml_ef, le_auc, le_ef, max(r2r, r2k)))
    report = {"n_ligands": len(rows), "n_features": int(X.shape[1]), "n_scaffolds": len(set(groups)),
              "prediction": {"ridge": {"r2": round(r2r, 3), "mae": round(maer, 3), "spearman": round(spr, 3)},
                             "knn10": {"r2": round(r2k, 3), "mae": round(maek, 3), "spearman": round(spk, 3)}},
              "enrichment": res,
              "caveat": "actives are defined by raw Vina score, so the raw_vina row is circular; "
                        "compare ML only against ligand_efficiency",
              "verdict": verdict,
              "gates": {g: sum(1 for r in out if r["gate"] == g) for g in ("pass", "flag", "fail")},
              "top10": [{k: r[k] for k in ("ligand_id", "affinity", "le", "pred_affinity", "mw",
                                           "clogp", "qed", "gate", "ai_score")} for r in out[:10]]}
    json.dump(report, open(os.path.join(ROOT, "demo", "ai_report.json"), "w"), indent=2)
    print("\n[3] verdict: %s" % verdict)
    print("    approved-drug gates: %s" % report["gates"])
    print("out -> demo/ai_report.json, demo/shortlist.tsv")


if __name__ == "__main__":
    main()
