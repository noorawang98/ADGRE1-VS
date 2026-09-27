#!/usr/bin/env python3
"""Turn the downloaded ZINC15 SMILES into ONE compact property/alerts table, then join to
the docking queue. No 3D files are generated (that is the storage-cheap path).

Steps
  1. stream filters/../zinc15_smiles/by_tranche/*.smi (SMILES<TAB>ZINC_ID)
  2. RDKit descriptors: MW, cLogP, TPSA, HBD, HBA, RotB, Fsp3, HAC, QED
  3. alerts: PAINS_A/B/C + Brenk + 10 toxicophore SMARTS (nitro-aryl, Michael acceptor,
     epoxide, hydrazine, acyl halide, aniline, aldehyde, azide, peroxide, nitroso)
  4. approved-drug gates from filters/approved_drug_thresholds.json (p5-p95 windows)
  5. reconcile with the local PDBQT ids -> which ligands have NO SMILES (degraded mode)
Output: filters/library_props.tsv.gz  (one gzip table, ~40-80 MB for 1.4M rows)

  ~/anaconda3/envs/AmberTools24/bin/python filters/scripts/prepare_library_props.py \
      --smiles-dir /home/workspace/MD/MD/zinc15_smiles/by_tranche
"""
import argparse, glob, gzip, json, os, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ZINC_ROOT = "/home/workspace/MD/MD/zinc"
ALERTS = {
    "nitro_aromatic": "[c;$(c[N+](=O)[O-])]",
    "michael_acceptor": "[CX3]=[CX3][CX3]=[OX1]",
    "epoxide": "C1OC1",
    "hydrazine": "[NX3][NX3]",
    "acyl_halide": "[CX3](=O)[F,Cl,Br,I]",
    "alkyl_halide": "[CX4][F,Cl,Br,I]",
    "aniline": "[NX3;H2,H1;!$(NC=O)]c",
    "aldehyde": "[CX3H1](=O)[#6]",
    "azide": "[NX2]=[NX2+]=[NX1-]",
    "peroxide": "[OX2][OX2]",
    "nitroso": "[NX2]=[OX1]",
    "thiol": "[SX2H]",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles-dir", required=True)
    ap.add_argument("--thresholds", default=os.path.join(ROOT, "filters", "approved_drug_thresholds.json"))
    ap.add_argument("--gate", default="approved_like_p5_p95")
    ap.add_argument("--hac-cap", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(ROOT, "filters", "library_props.tsv.gz"))
    args = ap.parse_args()

    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED, FilterCatalog
    RDLogger.DisableLog("rdApp.*")

    win = {}
    if os.path.exists(args.thresholds):
        thr = json.load(open(args.thresholds))
        win = thr["gates"].get(args.gate, {})
        print("gate=%s windows=%s" % (args.gate, win))
    else:
        print("WARN: %s missing -> no approved-drug window applied" % args.thresholds)

    params = FilterCatalog.FilterCatalogParams()
    for c in (FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A,
              FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B,
              FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C,
              FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK):
        params.AddCatalog(c)
    catalog = FilterCatalog.FilterCatalog(params)
    alert_patts = {k: Chem.MolFromSmarts(v) for k, v in ALERTS.items()}

    files = sorted(glob.glob(os.path.join(args.smiles_dir, "*.smi"))
                   + glob.glob(os.path.join(args.smiles_dir, "*.tsv")))
    if not files:
        sys.exit("no .smi/.tsv under %s" % args.smiles_dir)
    print("smiles files: %d" % len(files))

    cols = ["zinc_id", "smiles", "mw", "clogp", "clogp_zinc", "mw_zinc", "tpsa", "hbd",
            "hba", "rotb", "fsp3", "hac", "qed", "pains_brenk", "alerts", "gate"]
    n = kept = parse_fail = 0
    t0 = time.time()
    with gzip.open(args.out, "wt") as out:
        out.write("\t".join(cols) + "\n")
        for f in files:
            with open(f, errors="ignore") as fh:
                for line in fh:
                    n += 1
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 2 or p[0] == "smiles" or p[1] == "zinc_id":
                        continue                      # header row
                    smi, zid = p[0], p[1]
                    m = Chem.MolFromSmiles(smi)
                    if m is None:
                        parse_fail += 1
                        continue
                    hac = m.GetNumHeavyAtoms()
                    qed = QED.qed(m)
                    mw_zinc = p[3] if len(p) > 3 else ""
                    clogp_zinc = p[4] if len(p) > 4 else ""
                    pains = sorted({e.GetDescription() for e in catalog.GetMatches(m)})
                    fired = sorted(k for k, patt in alert_patts.items()
                                   if patt and m.HasSubstructMatch(patt))
                    d = dict(mw=Descriptors.MolWt(m), clogp=Crippen.MolLogP(m),
                             tpsa=rdMolDescriptors.CalcTPSA(m), hbd=rdMolDescriptors.CalcNumHBD(m),
                             hba=rdMolDescriptors.CalcNumHBA(m),
                             rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
                             fsp3=rdMolDescriptors.CalcFractionCSP3(m))
                    gate = "pass"
                    if win:
                        for k, (lo, hi) in win.items():
                            if k in d and not (lo <= d[k] <= hi):
                                gate = "fail"; break
                    if gate == "pass" and pains:
                        gate = "flag"
                    if gate == "pass" and fired:
                        gate = "flag"
                    if args.hac_cap and hac > args.hac_cap:
                        gate = "fail"
                    out.write("\t".join([zid, smi, "%.1f" % d["mw"], "%.2f" % d["clogp"],
                                         clogp_zinc, mw_zinc,
                                         "%.1f" % d["tpsa"], str(d["hbd"]), str(d["hba"]),
                                         str(d["rotb"]), "%.2f" % d["fsp3"], str(hac),
                                         "%.3f" % qed, ",".join(pains), ",".join(fired), gate]) + "\n")
                    kept += 1
            print("  %s -> cumulative %d rows (%.0f s)" % (os.path.basename(f), kept, time.time() - t0))
    print("rows in=%d written=%d parse_fail=%d" % (n, kept, parse_fail))
    print("table -> %s (%.1f MB)" % (os.path.relpath(args.out, ROOT), os.path.getsize(args.out) / 1e6))

    # reconcile against the local PDBQT library
    if os.path.isdir(ZINC_ROOT):
        have = set()
        for d in glob.glob(os.path.join(ZINC_ROOT, "ZINC15_*")):
            have.update(os.path.basename(f)[:-6] for f in glob.glob(os.path.join(d, "*.pdbqt")))
        import csv
        got = set()
        with gzip.open(args.out, "rt") as fh:
            next(fh)
            for line in fh:
                got.add(line.split("\t", 1)[0])
        miss = have - got
        print("pdbqt ids=%d  with SMILES=%d  without=%d" % (len(have), len(have & got), len(miss)))
        with open(os.path.join(ROOT, "filters", "ids_without_smiles.txt"), "w") as fh:
            for z in sorted(miss):
                fh.write(z + "\n")
        print("ligands without SMILES (must use the bond-order-free QC path) -> "
              "filters/ids_without_smiles.txt")


if __name__ == "__main__":
    main()
