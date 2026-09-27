#!/usr/bin/env python3
"""Pocket-assignment triage for ADGRE1 VS outputs, with receptor-contact discrimination.
CPU only, stdlib only.

Why contact matters: the 80_316 box (80x60x140 A) and the 748_852 box (60x40x40 A) are
larger than the receptor itself, so a pure 'centre inside box' test puts every one of the
2102 ligands into all three pockets (observed false positive). A pose only counts as
pocket-bound if it both sits in the box and actually contacts the receptor surface.
"""
import argparse, csv, glob, math, os, sys

BOXES = {  # name: (cx, cy, cz, sx, sy, sz) from config/config_*.arg
    "811":     (-65.392, -23.277, 69.679, 10.0, 10.0, 5.0),
    "748_852": (-77.917, -29.517, 77.864, 60.0, 40.0, 40.0),
    "80_316":  (101.425, 26.972, -140.473, 80.0, 60.0, 140.0),
}
CONTACT = 4.5      # A, ligand-receptor heavy-atom contact threshold
MAX_CONTACT = 4.5


def read_receptor(path):
    A = []
    with open(path, errors="ignore") as fh:
        for l in fh:
            if l.startswith(("ATOM", "HETATM")):
                try:
                    A.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                except ValueError:
                    continue
    return A


def poses(path, max_models=9):
    """yield (affinity, [ (x,y,z), ... ] atoms ) per MODEL (falls back to one model)."""
    aff, atoms, out = None, [], []
    with open(path, errors="ignore") as fh:
        for l in fh:
            if l.startswith("REMARK VINA RESULT"):
                try: aff = float(l.split()[3])
                except (IndexError, ValueError): aff = None
            elif l.startswith(("ATOM", "HETATM")):
                try: atoms.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                except ValueError: continue
            elif l.startswith("ENDMDL"):
                if atoms and aff is not None: out.append((aff, atoms))
                atoms, aff = [], None
                if len(out) >= max_models: break
    if atoms and aff is not None: out.append((aff, atoms))
    return out


def centroid(atoms):
    n = len(atoms)
    return (sum(a[0] for a in atoms)/n, sum(a[1] for a in atoms)/n, sum(a[2] for a in atoms)/n)


def in_box(c, box, tol=0.0):
    cx, cy, cz, sx, sy, sz = box
    return (abs(c[0]-cx) <= sx/2+tol and abs(c[1]-cy) <= sy/2+tol and abs(c[2]-cz) <= sz/2+tol)


def min_contact(atoms, rec, cell=6.0, budget=400):
    """min heavy-atom distance ligand<->receptor, using a coarse grid over receptor atoms."""
    best = 99.0
    for a in atoms[:budget]:
        for r in rec:
            d = math.dist(a, r)
            if d < best:
                best = d
                if best < 1.5: return best
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--dirs", nargs="*", default=None)
    ap.add_argument("--receptor", default=None, help="pdbqt of the receptor in the same frame")
    ap.add_argument("--max-affinity", type=float, default=-7.0)
    ap.add_argument("--outdir", default=None, help="output directory (default: <repo>/triage)")
    args = ap.parse_args()
    root = args.root
    out = args.outdir or os.path.join(root, "triage"); os.makedirs(out, exist_ok=True)

    rec_path = args.receptor or os.path.join(root, "811", "receptor.pdbqt")
    rec = read_receptor(rec_path)
    print("receptor atoms: %d (%s)" % (len(rec), os.path.relpath(rec_path, root)))

    srcs = []
    if args.dirs:
        for d in args.dirs:
            srcs.extend(glob.glob(os.path.join(root, d, "*.pdbqt")))
    else:
        for pat in ("*/out/*.pdbqt", "*/*_out/*.pdbqt"):
            srcs.extend(glob.glob(os.path.join(root, pat)))
    srcs = sorted(set(srcs))
    if not srcs: sys.exit("no *_out.pdbqt under %s" % root)

    best = {}   # (pocket, ligand) -> record, keeping the best contacting pose
    stats = dict(files=0, poses=0, clash=0, no_contact=0, kept=0)
    for path in srcs:
        stats["files"] += 1
        lid = os.path.basename(path).replace("_out.pdbqt", "")
        for aff, atoms in poses(path):
            stats["poses"] += 1
            if aff < -50 or aff > 0:          # clash artefact (+1e6, -107)
                stats["clash"] += 1; continue
            c = centroid(atoms)
            hit = [k for k, b in BOXES.items() if in_box(c, b)]
            if not hit: continue
            d = min_contact(atoms, rec)
            if d > MAX_CONTACT:
                stats["no_contact"] += 1; continue
            for k in hit:
                key = (k, lid)
                if key not in best or aff < best[key][1]:
                    best[key] = (c, aff, d, os.path.relpath(path, root))
                    stats["kept"] += 1
    print("files=%d poses=%d clash_dropped=%d no_contact_dropped=%d kept=%d"
          % (stats["files"], stats["poses"], stats["clash"], stats["no_contact"], stats["kept"]))

    keys = ["ligand_id", "best_aff", "min_contact_A", "pockets", "n_pockets", "source_file"]
    per = {k: [] for k in BOXES}
    lig_pockets = {}
    for (k, lid), (c, aff, d, src) in best.items():
        lig_pockets.setdefault(lid, {})[k] = aff
    for lid, d_ in lig_pockets.items():
        for k, aff in d_.items():
            per[k].append({"ligand_id": lid, "best_aff": round(aff, 3),
                           "min_contact_A": round(best[(k, lid)][2], 2),
                           "pockets": ";".join(sorted(d_)), "n_pockets": len(d_),
                           "source_file": best[(k, lid)][3]})
    for k in BOXES:
        sel = sorted(per[k], key=lambda r: r["best_aff"])
        p = os.path.join(out, "candidates_%s.tsv" % k)
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t"); w.writeheader()
            for r in sel: w.writerow(r)
        good = [r for r in sel if r["best_aff"] <= args.max_affinity]
        print("%-8s contacting=%5d  aff<=%.1f: %4d  best=%.2f -> triage/candidates_%s.tsv"
              % (k, len(sel), args.max_affinity, len(good), sel[0]["best_aff"] if sel else float("nan"), k))

    cons = [{"ligand_id": lid, "best_aff": round(min(d_.values()), 3), "min_contact_A": "",
             "pockets": ";".join(sorted(d_)), "n_pockets": len(d_), "source_file": ""}
            for lid, d_ in lig_pockets.items() if len(d_) >= 2]
    cons.sort(key=lambda r: r["best_aff"])
    with open(os.path.join(out, "consensus_multipocket.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t"); w.writeheader()
        for r in cons: w.writerow(r)
    print("multi-pocket consensus: %d -> triage/consensus_multipocket.tsv" % len(cons))


if __name__ == "__main__":
    main()
