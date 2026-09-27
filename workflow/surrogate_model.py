#!/usr/bin/env python3
"""CPU-only affinity surrogate + decoy-aware size correction for the 811 pocket.

Trains on the existing 2102-molecule fda screen (or any pocket's out/ dir) using
RDKit Morgan fingerprints (SVD-compressed), predicts a best-pose affinity for the
whole ZINC15 library, and reports which existing ligands are merely size/desolvation
favoured (score explained by descriptors alone) rather than pocket-specific.

Run inside the adfr or AmberTools24 env (RDKit + sklearn):
  ~/anaconda3/envs/AmberTools24/bin/python train_surrogate.py --pocket 811 \
      --zinc-root /home/workspace/MD/MD/zinc --out 811_triage
"""
import argparse, glob, json, os, sys
import numpy as np

FP_BITS, FP_RADIUS = 2048, 2
LEAD_RULES = dict(mw=(280, 500), clogp=(1.0, 4.0), tpsa=(0, 95), hbd=(0, 3), rotb=(0, 6))


def props(mol, Descriptors, Crippen, rdMolDescriptors):
    return dict(
        mw=Descriptors.MolWt(mol),
        clogp=Crippen.MolLogP(mol),
        tpsa=rdMolDescriptors.CalcTPSA(mol),
        hbd=rdMolDescriptors.CalcNumHBD(mol),
        hba=rdMolDescriptors.CalcNumHBA(mol),
        rotb=rdMolDescriptors.CalcNumRotatableBonds(mol),
        rings=rdMolDescriptors.CalcNumRings(mol),
        hac=mol.GetNumHeavyAtoms(),
    )


def to_smiles(pdbqt):
    """Extract SMILES-free fallback: use RDKit on the pdbqt converted through openbabel."""
    import subprocess, tempfile
    tmp = tempfile.mktemp(suffix=".smi")
    try:
        subprocess.run(["obabel", pdbqt, "-osmi", "-O", tmp], capture_output=True, timeout=30)
        with open(tmp) as fh:
            line = fh.readline().split()
        return line[0] if line else None
    except Exception:
        return None
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pocket", default="811")
    ap.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--zinc-root", default="/home/workspace/MD/MD/zinc")
    ap.add_argument("--out", default=None)
    ap.add_argument("--top", type=int, default=500)
    ap.add_argument("--max-train", type=int, default=2200)
    args = ap.parse_args()
    outdir = args.out or os.path.join(args.root, args.pocket + "_triage")
    os.makedirs(outdir, exist_ok=True)

    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, AllChem
    from rdkit.Chem import DataStructs
    from sklearn.decomposition import TruncatedSVD
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.model_selection import cross_val_score
    RDLogger.DisableLog("rdApp.*")

    # ---- 1. training set: existing fda screen for this pocket ----
    rows = []
    for f in sorted(glob.glob(os.path.join(args.root, args.pocket, "out", "*.pdbqt")))[:args.max_train]:
        best = None
        with open(f, errors="ignore") as fh:
            for line in fh:
                if line.startswith("REMARK VINA RESULT"):
                    try: v = float(line.split()[3])
                    except (IndexError, ValueError): continue
                    best = v if best is None else min(best, v)
                elif line.startswith("ENDMDL"):
                    break
        if best is None or best > 0 or best < -50:   # drop clash artefacts (+1e6, -107)
            continue
        smi = to_smiles(f)
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        rows.append((os.path.basename(f).replace("_out.pdbqt", ""), mol, best))
    print("training ligands (clash-filtered): %d" % len(rows))
    if len(rows) < 100:
        sys.exit("too few training rows; check obabel and %s/%s/out" % (args.root, args.pocket))

    gen = AllChem.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)
    X_fp = np.zeros((len(rows), FP_BITS), dtype=np.float32)
    X_pr, y = [], []
    for i, (lid, mol, aff) in enumerate(rows):
        DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), X_fp[i])
        p = props(mol, Descriptors, Crippen, rdMolDescriptors)
        X_pr.append([p["mw"], p["clogp"], p["tpsa"], p["hbd"], p["hba"], p["rotb"], p["rings"], p["hac"]])
        y.append(aff)
    X_pr, y = np.array(X_pr), np.array(y)

    svd = TruncatedSVD(n_components=64, random_state=0)
    Xf = svd.fit_transform(X_fp)
    X = np.hstack([Xf, X_pr])
    print("feature matrix: %s (fp64 SVD + 8 props)" % (X.shape,))

    # ---- 2. two models: full vs property-only (size/desolvation baseline) ----
    full = GradientBoostingRegressor(random_state=0, n_estimators=400, max_depth=3, learning_rate=0.05)
    base = GradientBoostingRegressor(random_state=0, n_estimators=400, max_depth=3, learning_rate=0.05)
    cv_full = cross_val_score(full, X, y, cv=5, scoring="r2")
    cv_base = cross_val_score(base, X_pr, y, cv=5, scoring="r2")
    print("5-fold R2  full=%.3f+-%.3f   property-only=%.3f+-%.3f"
          % (cv_full.mean(), cv_full.std(), cv_base.mean(), cv_base.std()))
    full.fit(X, y); base.fit(X_pr, y)

    # ---- 3. score the lead-like subset of the ZINC15 library by SMILES ----
    # ligand pdbqt dirs contain no SMILES, so convert a capped sample per shard
    cands = []
    shards = sorted(glob.glob(os.path.join(args.zinc_root, "ZINC15_*")))
    for sh in shards:
        for f in sorted(glob.glob(os.path.join(sh, "*.pdbqt")))[:args.top]:
            smi = to_smiles(f)
            if not smi: continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None: continue
            p = props(mol, Descriptors, Crippen, rdMolDescriptors)
            ok = (LEAD_RULES["mw"][0] <= p["mw"] <= LEAD_RULES["mw"][1]
                  and LEAD_RULES["clogp"][0] <= p["clogp"] <= LEAD_RULES["clogp"][1]
                  and p["tpsa"] <= LEAD_RULES["tpsa"][1]
                  and p["hbd"] <= LEAD_RULES["hbd"][1]
                  and p["rotb"] <= LEAD_RULES["rotb"][1])
            if not ok: continue
            cands.append((os.path.basename(f).replace(".pdbqt", ""), mol, p))
        if len(cands) >= 20000:
            break
    print("lead-like library candidates collected: %d" % len(cands))
    if not cands:
        sys.exit("no candidates passed the property gate")

    F = np.zeros((len(cands), FP_BITS), dtype=np.float32)
    PR = np.zeros((len(cands), 8), dtype=np.float32)
    for i, (lid, mol, p) in enumerate(cands):
        DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), F[i])
        PR[i] = [p["mw"], p["clogp"], p["tpsa"], p["hbd"], p["hba"], p["rotb"], p["rings"], p["hac"]]
    XC = np.hstack([svd.transform(F), PR])
    pred = full.predict(XC)
    pred_base = base.predict(PR)
    le = -pred / PR[:, 7]                       # ligand efficiency from predicted affinity

    order = np.argsort(pred)
    with open(os.path.join(outdir, "surrogate_shortlist.tsv"), "w") as fh:
        fh.write("rank\tligand_id\tpred_aff\tpred_prop_only\tspecificity_delta\tpred_LE\tmw\tclogp\ttpsa\thbd\trotb\thac\n")
        for r, i in enumerate(order[:args.top], 1):
            lid, _, p = cands[i]
            fh.write("%d\t%s\t%.2f\t%.2f\t%.2f\t%.3f\t%.1f\t%.2f\t%.1f\t%d\t%d\t%d\n"
                     % (r, lid, pred[i], pred_base[i], pred[i] - pred_base[i], le[i],
                        p["mw"], p["clogp"], p["tpsa"], p["hbd"], p["rotb"], p["hac"]))
    json.dump({"n_train": len(rows), "cv_r2_full": float(cv_full.mean()),
               "cv_r2_property_only": float(cv_base.mean()), "n_candidates": len(cands)},
              open(os.path.join(outdir, "surrogate_metrics.json"), "w"), indent=2)

    # ---- 4. flag training ligands whose score is pure size effect ----
    tr_pred = full.predict(X); tr_base = base.predict(X_pr)
    with open(os.path.join(outdir, "training_specificity.tsv"), "w") as fh:
        fh.write("ligand_id\taff_true\taff_pred\taff_prop_only\tspecificity_delta\n")
        for (lid, _, aff), a, b in zip(rows, tr_pred, tr_base):
            fh.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\n" % (lid, aff, a, b, a - b))
    print("outputs -> %s/{surrogate_shortlist,training_specificity}.tsv, surrogate_metrics.json" % outdir)


if __name__ == "__main__":
    main()
