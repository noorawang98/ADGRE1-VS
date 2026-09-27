#!/usr/bin/env python3
"""Prepare train/test splits (and later a library file) for the deep-scoring demo.

Two label sources, and the difference matters for what you can claim:
  A) docking labels  : demo/dataset_811.tsv (2630 rows, Vina score on pocket 811/the 600-852
     surface site). Training on these makes the model a *surrogate of Vina* -- a legitimate
     ML exercise, but it cannot be presented as improving on docking, because the target
     simply is the docking score.
  B) activity labels : requires measured IC50/Kd/EC50. If you have any, put them in
     deep_demo/input/activity.tsv with columns (ligand_id, smiles, value, units, assay) and
     this script will prefer them. ONLY (B) makes "AI screening" a defensible claim.

Split policy: scaffold-grouped 5-fold (Bemis-Murcko when RDKit is present, else a
ring-system key), so no scaffold appears in both train and test. Writes Chemprop- and
Uni-Mol-friendly files (SMILES<TAB>target).

  ~/anaconda3/envs/AmberTools24/bin/python deep_demo/scripts/prepare_deep_inputs.py
-> deep_demo/input/{train.csv,val.csv,test.csv,folds.csv,library.smi}
"""
import argparse, csv, os, random, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def scaffold(smiles, use_rdkit=True):
    if use_rdkit:
        try:
            from rdkit import Chem, RDLogger
            from rdkit.Chem.Scaffolds import MurckoScaffold
            RDLogger.DisableLog("rdApp.*")
            m = Chem.MolFromSmiles(smiles)
            if m is not None:
                s = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
                if s:
                    return s
        except ImportError:
            pass
    # fallback: ring-system key (ring-closure digits + enclosed atoms)
    pairs = {}
    for i, ch in enumerate(smiles):
        if ch.isdigit():
            pairs.setdefault(ch, []).append(i)
    rings = [smiles[p[0]:p[-1] + 1] for p in pairs.values() if len(p) >= 2]
    return "|".join(sorted(rings))[:120] or smiles[:40]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=os.path.join(ROOT, "demo", "dataset_811.tsv"))
    ap.add_argument("--activity", default=os.path.join(ROOT, "deep_demo", "input", "activity.tsv"))
    ap.add_argument("--library", default=os.path.join(ROOT, "deep_demo", "input", "library.tsv"),
                    help="optional: id<TAB>smiles table to write as library.smi for inference")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default=os.path.join(ROOT, "deep_demo", "input"))
    args = ap.parse_args()
    random.seed(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    rows, source = [], "docking"
    if os.path.exists(args.activity) and os.path.getsize(args.activity) > 0:
        with open(args.activity) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("smiles") and r.get("value"):
                    rows.append({"ligand_id": r["ligand_id"], "smiles": r["smiles"],
                                 "target": -float(r["value"])})   # lower = better, Vina convention
        source = "activity"
        print("using MEASURED activity labels: %d rows (this is what makes 'AI screening' defensible)"
              % len(rows))
    if not rows:
        with open(args.dataset) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("smiles"):
                    rows.append({"ligand_id": r["ligand_id"], "smiles": r["smiles"],
                                 "target": float(r["affinity"])})
        print("using DOCKING labels: %d rows -- the model can only reproduce Vina; "
              "do not claim it beats docking" % len(rows))
    if len(rows) < 100:
        sys.exit("need >=100 labelled rows")

    # scaffold-grouped folds
    groups = {}
    for r in rows:
        groups.setdefault(scaffold(r["smiles"]), []).append(r)
    keys = sorted(groups)
    random.shuffle(keys)
    fold_of = {k: i % args.folds for i, k in enumerate(keys)}
    for r in rows:
        r["fold"] = fold_of[scaffold(r["smiles"])]
    print("scaffolds: %d -> folds: %s" % (len(keys),
          {f: sum(1 for r in rows if r["fold"] == f) for f in range(args.folds)}))

    # fold 0 = test; within the remaining, hold out 10% of scaffolds as validation
    test = [r for r in rows if r["fold"] == 0]
    rest = [r for r in rows if r["fold"] != 0]
    val_keys = [k for k in keys if fold_of[k] != 0][:max(1, len(keys) // 10)]
    val = [r for r in rest if scaffold(r["smiles"]) in set(val_keys)]
    train = [r for r in rest if scaffold(r["smiles"]) not in set(val_keys)]

    def write(name, data, extra=False):
        p = os.path.join(args.outdir, name)
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh)
            cols = ["smiles", "target"] + (["ligand_id", "fold"] if extra else [])
            w.writerow(cols)
            for r in data:
                w.writerow([r["smiles"], r["target"]] + ([r["ligand_id"], r["fold"]] if extra else []))
        print("  %-14s %5d rows -> %s" % (name, len(data), os.path.relpath(p, ROOT)))
    print("split (scaffold-grouped):")
    write("train.csv", train)
    write("val.csv", val)
    write("test.csv", test)
    write("folds.csv", rows, extra=True)

    # library for inference (no labels)
    lib_src = args.library if os.path.exists(args.library) else None
    if lib_src:
        with open(lib_src) as fh:
            n = 0
            with open(os.path.join(args.outdir, "library.smi"), "w") as o:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("smiles"):
                        o.write("%s\t%s\n" % (r["smiles"], r.get("ligand_id", "")))
                        n += 1
        print("  library.smi    %5d rows (labels stripped, for inference)" % n)
    else:
        print("  library.tsv absent -> write id<TAB>smiles there to generate library.smi")

    print("\nlabel source: %s | folds: %d | train/val/test = %d/%d/%d"
          % (source.upper(), args.folds, len(train), len(val), len(test)))


if __name__ == "__main__":
    main()
