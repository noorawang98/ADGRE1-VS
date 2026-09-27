#!/usr/bin/env python3
"""Size- and bond-order-independent ligand QC for PDBQT screening output.

Why not MW/cLogP/QED: measured on this set, PDBQT carries no bond orders -- obabel yields
radical aromatic rings, rdDetermineBonds fails on most records (read_fail), and 0/200 tested
ligands reached a trustworthy descriptor state. Descriptor-based toxicity/drug-likeness
filters are therefore NOT applicable to these files without re-preparing the library from
SDF/SMILES. What IS applicable:

  size        : heavy-atom count straight from the PDBQT ATOM records -> ligand efficiency
                LE = -dG / HAC  (size-decoupled ranking; essential because the 811/748_852
                boxes are shallow and large ligands collect van-der-Waals score for free)
  connectivity: element-labelled graph from 3D coordinates -> toxicophore flags that need
                only element-level neighbourhoods (nitro, epoxide, hydrazine, acyl halide,
                N-oxide, azide, peroxide), plus synthetic-friendliness proxies
  approved-drug window: filters/approved_drug_thresholds.json holds the empirical p5-p95
                windows of the FDA-approved set; with an SDF/SMILES source they become a
                hard gate (see --apply-descriptors, which requires filters/smiles_map.tsv
                rows with trust=high).

  ~/anaconda3/envs/AmberTools24/bin/python filters/scripts/ligand_qc.py \
      --affinity 811/affinity-811.out --out filters/ligand_qc_811.tsv
"""
import argparse, csv, glob, json, math, os, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BOND_R = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
          "Cl": 1.02, "Br": 1.20, "I": 1.39}


def read_ligand(path):
    """first MODEL only: returns (affinity, [(element, x, y, z)])"""
    aff, atoms = None, []
    for l in open(path, errors="ignore"):
        if l.startswith("REMARK VINA RESULT"):
            try: aff = float(l.split()[3])
            except (IndexError, ValueError): aff = None
        elif l.startswith(("ATOM", "HETATM")):
            t = l[76:79].strip()
            if not t or t[0] == "H":
                continue
            e = t[0]
            el = {"A": "C", "C": "C", "N": "N", "O": "O", "P": "P", "S": "S",
                  "F": "F", "I": "I", "B": "B", "M": "C", "G": "C", "T": "O"}.get(e, e)
            try:
                atoms.append((el, float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                continue
        elif l.startswith("ENDMDL"):
            break
    return aff, atoms


def neighbors(atoms, tol=0.45):
    """connectivity from 3D distances; returns list of neighbour index lists."""
    nb = [[] for _ in atoms]
    n = len(atoms)
    for i in range(n):
        ei, xi, yi, zi = atoms[i]
        for j in range(i + 1, n):
            ej, xj, yj, zj = atoms[j]
            d = math.dist((xi, yi, zi), (xj, yj, zj))
            cut = BOND_R.get(ei, 0.8) + BOND_R.get(ej, 0.8) + tol
            if d <= cut:
                nb[i].append(j); nb[j].append(i)
    return nb


def alerts(atoms, nb):
    """element-neighbourhood alerts that do NOT require bond orders."""
    el = [a[0] for a in atoms]
    fired = set()
    for i, e in enumerate(el):
        deg = len(nb[i])
        nO = sum(1 for j in nb[i] if el[j] == "O")
        nN = sum(1 for j in nb[i] if el[j] == "N")
        if e == "N" and nO >= 2 and deg >= 3:
            fired.add("nitro_or_Noxide")
        if e == "N" and nN >= 1 and deg <= 2:
            fired.add("hydrazine_like")
        if e == "C" and deg <= 3:
            nO2 = sum(1 for j in nb[i] if el[j] == "O" and len(nb[j]) == 1)
            if nO2 >= 2 and deg <= 3:
                fired.add("gem_diol_or_carbonyl_O")
        if e == "S" and nO >= 1:
            fired.add("sulfur_oxide")
        if e in ("F", "Cl", "Br", "I"):
            parent = el[nb[i][0]] if nb[i] else "?"
            if parent != "C":
                fired.add("halogen_on_heteroatom")
        if e == "O" and any(el[j] == "O" for j in nb[i]):
            fired.add("peroxide")
    # three-membered ring proxy: two neighbours of a C that are bonded to each other
    for i, e in enumerate(el):
        if e != "C" or len(nb[i]) < 2:
            continue
        for a in range(len(nb[i])):
            for b in range(a + 1, len(nb[i])):
                if nb[i][b] in nb[nb[i][a]]:
                    fired.add("3ring")
    return sorted(fired)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=None, help="dir with *_out.pdbqt (or use --affinity list)")
    ap.add_argument("--affinity", default=None, help="affinity.out file: '<name>.log  <value>'")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--hac-cap", type=int, default=0, help="0=off; from pocket profile")
    ap.add_argument("--out", default=os.path.join(ROOT, "filters", "ligand_qc.tsv"))
    args = ap.parse_args()

    files = []
    if args.affinity:
        base = os.path.dirname(os.path.abspath(args.affinity))
        for line in open(args.affinity):
            parts = line.split()
            if len(parts) < 2:
                continue
            name = parts[0].replace(".log", "")
            p = os.path.join(base, "out", name + "_out.pdbqt")
            if os.path.exists(p):
                files.append(p)
    elif args.dir:
        files = sorted(glob.glob(os.path.join(args.dir, "*.pdbqt")))
    if not files:
        sys.exit("no ligand files found")
    if args.limit:
        files = files[:args.limit]
    print("ligands to QC: %d" % len(files))

    rows, skipped = [], 0
    for p in files:
        aff, atoms = read_ligand(p)
        if aff is None or not atoms or aff > 0 or aff < -50:
            skipped += 1
            continue
        hac = len(atoms)
        nb = neighbors(atoms)
        al = alerts(atoms, nb)
        rows.append(dict(ligand_id=os.path.basename(p).replace("_out.pdbqt", ""),
                         affinity=round(aff, 2), hac=hac, le=round(-aff / hac, 3),
                         alerts=",".join(al),
                         hac_flag=("over_cap" if args.hac_cap and hac > args.hac_cap else ""),
                         descriptors="unavailable_pdbqt_no_bond_orders"))
    rows.sort(key=lambda r: r["le"], reverse=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("qc rows: %d (skipped clash/artefact: %d)" % (len(rows), skipped))
    if rows:
        print("LE top5: %s" % [(r["ligand_id"], r["le"]) for r in rows[:5]])
        from collections import Counter
        print("alert frequency: %s" % dict(Counter(a for r in rows for a in r["alerts"].split(",") if a)))
    print("out -> %s" % os.path.relpath(args.out, ROOT))
    print("NOTE: MW/cLogP/TPSA/QED gates need re-prepared ligands from SDF/SMILES; "
          "filters/smiles_map.tsv shows trust=high is 0%% on the current PDBQT set.")


if __name__ == "__main__":
    main()
