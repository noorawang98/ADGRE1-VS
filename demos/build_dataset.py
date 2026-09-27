#!/usr/bin/env python3
"""Stage 1: build the demo dataset = 811 docking scores + trustworthy SMILES + descriptors.

Joins three things that already exist in the workspace:
  * 811/affinity-811.out            : per-ligand best Vina score for pocket 811 (FDA subset)
  * demo/ligand_smiles.tsv          : SMILES recovered by skeleton-matching drugbank.sdf
  * /home/workspace/MD/MD/drugbank_mk_pdbqt : heavy-atom counts (bond-order independent)

  ~/anaconda3/envs/AmberTools24/bin/python demo/build_dataset.py
-> demo/dataset_811.tsv
"""
import argparse, csv, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ALERTS = {
    "nitro_aromatic": "[c;$(c[N+](=O)[O-])]", "michael_acceptor": "[CX3]=[CX3][CX3]=[OX1]",
    "epoxide": "C1OC1", "hydrazine": "[NX3][NX3]", "acyl_halide": "[CX3](=O)[F,Cl,Br,I]",
    "aniline": "[NX3;H2,H1;!$(NC=O)]c", "aldehyde": "[CX3H1](=O)[#6]",
    "azide": "[NX2]=[NX2+]=[NX1-]", "peroxide": "[OX2][OX2]", "nitroso": "[NX2]=[OX1]",
}


def pdbqt_hac(path):
    n = 0
    for l in open(path, errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            t = l[76:79].strip()
            if t and t[0] != "H":
                n += 1
        elif l.startswith("ENDMDL"):
            break
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--affinity", default=os.path.join(ROOT, "811", "affinity-811.out"),
                    help="legacy '<name>.log <score>' list; superseded by --scan-dir")
    ap.add_argument("--scan-dir", default=os.path.join(ROOT, "cpu_scan", "811_b24"),
                    help="preferred score source: CPU re-scan with the corrected 24 A box")
    ap.add_argument("--smiles", default=os.path.join(ROOT, "demo", "ligand_smiles.tsv"))
    ap.add_argument("--ligdir", default="/home/workspace/MD/MD/drugbank_mk_pdbqt")
    ap.add_argument("--fda-pdbqt", default=os.path.expanduser("~/scripts/vinalc_test/data/fda"))
    ap.add_argument("--out", default=os.path.join(ROOT, "demo", "dataset_811.tsv"))
    ap.add_argument("--ligand-table", default=os.path.join(ROOT, "demo", "ligand_smiles.tsv"),
                    help="ligand_id/smiles table (from recover_smiles.py)")
    args = ap.parse_args()

    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED
    RDLogger.DisableLog("rdApp.*")
    alert_patts = {k: Chem.MolFromSmarts(v) for k, v in ALERTS.items()}

    # trustworthy SMILES only: descriptor-based gates (MW/cLogP/QED) are invalid without
    # bond orders, and mixing in degraded rows would corrupt the approved-drug window test.
    # no_match ligands are handled separately by ligand_qc.py (heavy atoms + alerts).
    smap = {}
    degraded = 0
    with open(args.smiles) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("smiles") and r.get("source") == "sdf_match":
                smap[r["ligand_id"]] = r["smiles"]
            elif r.get("source") == "no_match":
                degraded += 1
    print("trustworthy SMILES: %d | degraded (no_match): %d" % (len(smap), degraded))

    # Score source. The original 811 run used size_z=5, which produced clash artefacts
    # (median -0.05 kcal/mol); the corrected 24 A rescan is the trustworthy source, so it is
    # preferred and the legacy list is only a fallback.
    scores = {}
    if os.path.isdir(args.scan_dir):
        for f in os.listdir(args.scan_dir):
            if not f.endswith("_out.pdbqt"):
                continue
            aff = None
            for l in open(os.path.join(args.scan_dir, f), errors="ignore"):
                if l.startswith("REMARK VINA RESULT"):
                    try: aff = float(l.split()[3])
                    except (IndexError, ValueError): pass
                    break
            if aff is not None:
                scores[f.replace("_out.pdbqt", "")] = aff
        print("scores from %s: %d" % (os.path.relpath(args.scan_dir, ROOT), len(scores)))
    if not scores:
        for line in open(args.affinity):
            p = line.split()
            if len(p) >= 2:
                try:
                    scores[p[0].replace(".log", "")] = float(p[1])
                except ValueError:
                    pass
        print("scores from legacy list: %d (WARNING: size_z=5 artefacts)" % len(scores))

    rows, skipped = [], 0
    for lid in sorted(set(smap) & set(scores)):
        aff = scores[lid]
        if aff > 0 or aff < -50:          # clash artefacts (+1e6) and impossible scores
            skipped += 1
            continue
        m = Chem.MolFromSmiles(smap[lid])
        if m is None:
            skipped += 1
            continue
        hac_pdbqt = 0
        for d in (args.ligdir, args.fda_pdbqt):
            p = os.path.join(d, lid + ".pdbqt")
            if os.path.exists(p):
                hac_pdbqt = pdbqt_hac(p)
                break
        hac = hac_pdbqt or m.GetNumHeavyAtoms()
        fired = sorted(k for k, patt in alert_patts.items() if patt and m.HasSubstructMatch(patt))
        rows.append(dict(
            ligand_id=lid, smiles=smap[lid], affinity=round(aff, 3), hac=hac,
            le=round(-aff / hac, 4),
            mw=round(Descriptors.MolWt(m), 1), clogp=round(Crippen.MolLogP(m), 2),
            tpsa=round(rdMolDescriptors.CalcTPSA(m), 1),
            hbd=rdMolDescriptors.CalcNumHBD(m), hba=rdMolDescriptors.CalcNumHBA(m),
            rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
            fsp3=round(rdMolDescriptors.CalcFractionCSP3(m), 2),
            rings=rdMolDescriptors.CalcNumRings(m), qed=round(QED.qed(m), 3),
            alerts=",".join(fired)))
    if not rows:
        sys.exit("no joined rows; check that the matcher finished")
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["affinity"]):
            w.writerow(r)
    import statistics as st
    print("dataset rows: %d (skipped %d)" % (len(rows), skipped))
    print("affinity: min %.2f median %.2f | MW median %.0f | LE median %.3f"
          % (rows[0]["affinity"] if rows else 0,
             st.median([r["affinity"] for r in rows]), st.median([r["mw"] for r in rows]),
             st.median([r["le"] for r in rows])))
    print("out -> %s" % os.path.relpath(args.out, ROOT))


if __name__ == "__main__":
    main()
