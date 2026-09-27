#!/usr/bin/env python3
"""Step 2: characterise the pocket from step-1 CPU poses and emit a GPU-sized box.

Input : cpu_scan/<site>/*_out.pdbqt  (generous 24 A box, CPU vina)
Output: pocket profile (occupied volume, depth proxy, lining residues, ligand size spread)
        + one GPU-sized .arg with centre/edge fitted to the top-hit ensemble, capped by
        --box-limit because Vina-GPU caps the per-run search space.

Fitted box = centre of the top-N pose centroids, edge = min(box_limit,
round_up(2*(R90 + ligand radius))) so the box hugs the occupied pocket instead of the
original hand-set grid.
  ~/anaconda3/envs/adfr/bin/python fit_gpu_box.py --scan cpu_scan/811_b24 --top 100
"""
import argparse, glob, json, math, os
import numpy as np

P = {}


def pose(path):
    """best affinity + all atoms of the best model + centroid + max radius from centroid."""
    best, atoms, cur, aff = None, [], [], None
    for l in open(path, errors="ignore"):
        if l.startswith("REMARK VINA RESULT"):
            try: aff = float(l.split()[3])
            except (IndexError, ValueError): aff = None
        elif l.startswith(("ATOM", "HETATM")):
            try: cur.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError: continue
        elif l.startswith("ENDMDL"):
            if cur and aff is not None and -50 < aff < 0 and (best is None or aff < best[0]):
                best = (aff, np.asarray(cur))
            cur, aff = [], None
            break
    if best is None and cur and aff is not None and -50 < aff < 0:
        best = (aff, np.asarray(cur))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--scan", default="cpu_scan/811_b24")
    ap.add_argument("--receptor", default="811/receptor.pdbqt")
    ap.add_argument("--top", type=int, default=100, help="top-N poses used to fit the box")
    ap.add_argument("--box-limit", type=float, default=24.0, help="GPU max edge (A)")
    ap.add_argument("--out-tag", default=None)
    ap.add_argument("--outdir", default=None, help="output directory (default: <repo>/gpu_rerun/fitted_boxes)")
    args = ap.parse_args()
    root = args.root
    files = sorted(glob.glob(os.path.join(root, args.scan, "*_out.pdbqt")))
    if not files:
        raise SystemExit("no poses in %s" % args.scan)

    rec = []
    for l in open(os.path.join(root, args.receptor), errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            try: rec.append((float(l[30:38]), float(l[38:46]), float(l[46:54]),
                             "%s%s" % (l[17:20].strip()[:1], l[22:26].strip())))
            except ValueError: continue
    R = np.asarray([r[:3] for r in rec]); Rnames = [r[3] for r in rec]

    rows = []
    for f in files:
        b = pose(f)
        if b is None: continue
        aff, atoms = b
        c = atoms.mean(0)
        rad = float(np.linalg.norm(atoms - c, axis=1).max())
        rows.append((aff, c, rad, os.path.basename(f).replace("_out.pdbqt", ""), len(atoms)))
    rows.sort(key=lambda r: r[0])
    print("poses scored: %d | best %.2f | median %.2f"
          % (len(rows), rows[0][0], rows[len(rows)//2][0]))

    top = rows[:min(args.top, len(rows))]
    C = np.asarray([r[1] for r in top]); RAD = np.asarray([r[2] for r in top])
    centre = C.mean(0)
    R90 = float(np.percentile(np.linalg.norm(C - centre, axis=1) + RAD, 90))
    edge = min(args.box_limit, math.ceil(2 * (R90 + 1.5)))
    print("top-%d ensemble: centroid spread R90(centre+ligand radius)=%.1f A -> fitted edge %.0f A"
          % (len(top), R90, edge))

    centre = np.round(centre, 3)
    d = np.linalg.norm(R - centre, axis=1)
    lining = sorted(set(Rnames[i] for i in np.where(d <= edge / 2)[0]))
    depth_proxy = float(np.percentile(d, 5))
    print("centre=(%.2f, %.2f, %.2f)  lining residues within %.0f A: %d"
          % (centre[0], centre[1], centre[2], edge/2, len(lining)))
    print("lining:", ",".join(lining[:28]))
    print("ligand heavy atoms: median %.0f  range %d-%d  (use for MW gate)"
          % (np.median([r[4] for r in top]), min(r[4] for r in top), max(r[4] for r in top)))

    tag = args.out_tag or os.path.basename(args.scan.rstrip("/")) + "_fitted"
    ddir = args.outdir or os.path.join(root, "gpu_rerun", "fitted_boxes"); os.makedirs(ddir, exist_ok=True)
    arg = os.path.join(ddir, "config_%s.arg" % tag)
    with open(arg, "w") as fh:
        fh.write("center_x = %.3f\ncenter_y = %.3f\ncenter_z = %.3f\n" % tuple(centre))
        fh.write("size_x = %.1f\nsize_y = %.1f\nsize_z = %.1f\n" % ((edge,) * 3))
        fh.write("thread = 8000\nsearch_depth = 20\n")
    prof = {"scan": args.scan, "n_poses": len(rows),
            "best_affinity": round(rows[0][0], 2),
            "median_affinity": round(rows[len(rows)//2][0], 2),
            "top_used": len(top), "fitted_centre": [float(v) for v in centre],
            "fitted_edge_A": edge, "box_limit": args.box_limit,
            "lining_residues": lining, "ligand_heavy_atoms_median": float(np.median([r[4] for r in top])),
            "pocket_depth_proxy_A": round(depth_proxy, 2),
            "consistency_with_original_811_centre_A": round(float(np.linalg.norm(centre - np.array([-65.392, -23.277, 69.679]))), 2)}
    jp = os.path.join(ddir, "profile_%s.json" % tag)
    json.dump(prof, open(jp, "w"), indent=2)
    print("config -> %s" % os.path.relpath(arg, root))
    print("profile -> %s" % os.path.relpath(jp, root))
    print("note: calibrate CPU vs Vina-GPU scores on 100-200 shared ligands before merging queues")


if __name__ == "__main__":
    main()
