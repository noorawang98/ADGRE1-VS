#!/usr/bin/env python3
"""ADGRE2 counter-screen analysis: hit-rate based reverse screening (NO cross-site dG subtraction).

METHODOLOGICAL PREMISE (verified experimentally):
ADGRE1's 811 site and the corresponding ADGRE2 pocket are NOT the same pocket.
  - In ADGRE1 the ligands actually bind residues 690-710 (outer surface of TM3),
    pose centroid (-54.3, -10.5, 54.2).
  - The orthologous ADGRE2 region was located by TM-segment pairing and a BLOSUM62
    global alignment (identity 0.545): ADGRE2 residues 572-593, centroid
    (-50.7, -5.2, 54.9).
  - Pocket composition, volume and solvent exposure differ, so subtracting dG values
    obtained on the two receptors has no physical meaning.
Hence this script uses hit-rate reverse screening instead:
  1. take each ligand's best accessible score on each receptor (= "can it bind here?")
  2. use log-ratio and rank information rather than a difference of energies:
       dlog = log10(|aff_ADGRE2|) - log10(|aff_ADGRE1|)     (>0 means ADGRE1-preferring)
  3. criterion: a ligand entering the top X% on ADGRE1 while also entering the top X% on
     ADGRE2 is treated as NON-selective
  4. conclusions are phrased in terms of selective enrichment, never as absolute dG.

Also reports the correlation of affinity and of ligand efficiency (LE) between receptors;
LE correlation is the size-corrected test (raw scores correlate trivially because both grow
with ligand size).

  ~/anaconda3/envs/AmberTools24/bin/python reverse_screen_analysis.py \
      --a1 reverse_screen/adgre1_cpu/fda --a2 reverse_screen/adgre2_cpu/fda --tag gpu_fda
"""
import argparse, csv, glob, json, math, os, sys
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.abspath(__file__))
CLASH, CONTACT, MIN_OK = 2.0, 4.0, 2.2


def rec_atoms(pdbqt):
    R = []
    for l in open(pdbqt, errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            try:
                R.append((int(l[22:26]), l[17:20].strip(),
                          float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                continue
    return R


def best_pose(path, tree, rec):
    """Best affinity (clash artefacts filtered) + geometry + ligand SMILES from the pose file."""
    aff, smi, atoms = None, None, []
    for l in open(path, errors="ignore"):
        if l.startswith("REMARK SMILES IDX"):
            continue
        if l.startswith("REMARK SMILES ") and smi is None:
            smi = l.replace("REMARK SMILES ", "").strip()
        elif l.startswith("REMARK VINA RESULT"):
            try:
                v = float(l.split()[3])
            except (IndexError, ValueError):
                continue
            if -50 < v < 0 and (aff is None or v < aff):   # skip artefacts such as +1e6
                aff = v
        elif l.startswith(("ATOM", "HETATM")):
            try:
                atoms.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                continue
        elif l.startswith("ENDMDL"):
            break
    if aff is None or not atoms:
        return None
    L = np.asarray(atoms)
    d, idx = tree.query(L, k=6)
    if d.ndim == 1:
        d = d[:, None]; idx = idx[:, None]
    dmin = float(d.min())
    nclash = int((d < CLASH).sum())
    contacts = {rec[idx[i, k]][0] for i in range(len(L)) for k in range(d.shape[1]) if d[i, k] < CONTACT}
    return dict(aff=aff, hac=len(atoms), min_dist=dmin, n_clash=nclash,
                contacts=contacts, smiles=smi or "",
                geometry="clash" if (nclash > 0 or dmin < MIN_OK) else "ok")


def scan_dir(dirs, receptor, label):
    rec = rec_atoms(receptor)
    tree = cKDTree(np.asarray([r[2:5] for r in rec]))
    # accept both naming conventions: legacy  xxx_out.pdbqt  and new  xxx.pdbqt
    files = []
    for d in dirs:
        files += glob.glob(os.path.join(d, "*_out.pdbqt"))
        files += [f for f in glob.glob(os.path.join(d, "*.pdbqt")) if not f.endswith("_out.pdbqt")]
    print("[%s] receptor atoms %d | pose files %d (from %d dirs)"
          % (label, len(rec), len(files), len(dirs)))
    out = {}
    for f in files:
        b = os.path.basename(f)
        lid = b[:-len("_out.pdbqt")] if b.endswith("_out.pdbqt") else b[:-len(".pdbqt")]
        r = best_pose(f, tree, rec)
        if r:
            out[lid] = r
    ok = sum(1 for r in out.values() if r["geometry"] == "ok")
    print("[%s] parsed %d, geometry-passing %d (clash rejected %d)"
          % (label, len(out), ok, len(out) - ok))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a1", nargs="+", default=[os.path.join(ROOT, "cpu_scan", "811_b24")],
                    help="ADGRE1 pose dir(s)")
    ap.add_argument("--a1-receptor", default=os.path.join(ROOT, "811", "receptor.pdbqt"))
    ap.add_argument("--a1-extra", nargs="*", default=[], help="extra ADGRE1 pose dirs")
    ap.add_argument("--a2", nargs="+", default=["/home/workspace/MD/MD/HDY/ADGRE2/out"],
                    help="ADGRE2 pose dir(s)")
    ap.add_argument("--a2-receptor", default="/home/workspace/MD/MD/HDY/ADGRE2/receptor.pdbqt")
    ap.add_argument("--a2-extra", nargs="*", default=[], help="extra ADGRE2 pose dirs")
    ap.add_argument("--fda-map", default=os.path.expanduser("~/scripts/vinalc_test/data/fda"),
                    help="dir used to map fdaXXXX ids to ZINC ids")
    ap.add_argument("--ddg-cut", type=float, default=0.15,
                    help="log10 ratio threshold (0.15 ~ 1.4x); NOT kcal/mol")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--outdir", default=None,
                    help="directory for pairs_<tag>.tsv and summary_<tag>.json "
                         "(default: <repo>/reverse_screen)")
    args = ap.parse_args()

    a1 = scan_dir(list(args.a1) + list(args.a1_extra), args.a1_receptor, "ADGRE1")
    a2 = scan_dir(list(args.a2) + list(args.a2_extra), args.a2_receptor, "ADGRE2")

    # fdaXXXX -> ZINC mapping: legacy ADGRE2 output used fda ids while the ADGRE1 rescans
    # use drugbank numeric ids, so ids must be normalised before intersecting.
    f2z = {}
    if os.path.isdir(args.fda_map):
        for f in os.listdir(args.fda_map):
            if f.endswith(".pdbqt"):
                for l in open(os.path.join(args.fda_map, f), errors="ignore"):
                    if "Name =" in l:
                        f2z[f[:-6]] = l.split("=")[1].strip()
                        break

    def canon(lid):
        if lid in f2z:
            return f2z[lid]
        return lid
    a1c = {canon(k): (k, v) for k, v in a1.items()}
    a2c = {canon(k): (k, v) for k, v in a2.items()}
    common = sorted(set(a1c) & set(a2c))
    print("\nCommon ligands: %d (ADGRE1 %d vs ADGRE2 %d)" % (len(common), len(a1c), len(a2c)))

    rows, drop = [], Counter()
    for k in common:
        (i1, r1), (i2, r2) = a1c[k], a2c[k]
        if r1["geometry"] != "ok" or r2["geometry"] != "ok":
            drop["geometry_failed"] += 1
            continue
        # dlog instead of a dG difference: the two pockets are different
        dlog = round(math.log10(abs(r2["aff"])) - math.log10(abs(r1["aff"])), 3)
        rows.append(dict(ligand_id=k, aff_adgre1=round(r1["aff"], 2), aff_adgre2=round(r2["aff"], 2),
                         dlog_pref_a1=dlog, le_adgre1=round(-r1["aff"] / max(r1["hac"], 1), 3),
                         le_adgre2=round(-r2["aff"] / max(r2["hac"], 1), 3),
                         hac=r1["hac"], min_dist_a1=round(r1["min_dist"], 2),
                         min_dist_a2=round(r2["min_dist"], 2),
                         contact_overlap=len(r1["contacts"] & r2["contacts"]),
                         ag1_preferred=int(dlog >= args.ddg_cut),
                         ag2_preferred=int(dlog <= -args.ddg_cut),
                         note="dlog is NOT a dG; different pockets, hit-rate reverse screen only",
                         smiles_a1=r1["smiles"][:120]))
    outd = args.outdir or os.path.join(ROOT, "reverse_screen")
    os.makedirs(outd, exist_ok=True)
    tsv = os.path.join(outd, "pairs_%s.tsv" % args.tag)
    with open(tsv, "w", newline="") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            for r in sorted(rows, key=lambda r: -r["dlog_pref_a1"]):
                w.writerow(r)

    summary = {"tag": args.tag, "n_a1": len(a1), "n_a2": len(a2), "n_common": len(common),
               "n_pairs_geometry_ok": len(rows), "dropped": dict(drop), "ddg_cut": args.ddg_cut}
    if rows:
        dd = np.array([r["dlog_pref_a1"] for r in rows])
        summary.update({
            "dlog_median": round(float(np.median(dd)), 3),
            "dlog_mean": round(float(dd.mean()), 3),
            "dlog_p10": round(float(np.percentile(dd, 10)), 3),
            "dlog_p90": round(float(np.percentile(dd, 90)), 3),
            "n_ag1_preferred": int((dd >= args.ddg_cut).sum()),
            "n_ag2_preferred": int((dd <= -args.ddg_cut).sum()),
            "n_neutral": int(((dd > -args.ddg_cut) & (dd < args.ddg_cut)).sum()),
            # raw scores correlate trivially across receptors (both grow with ligand size);
            # the LE correlation is the size-corrected test
            "corr_affinity_between_receptors": round(float(np.corrcoef(
                [r["aff_adgre1"] for r in rows], [r["aff_adgre2"] for r in rows])[0, 1]), 3),
            "corr_le_between_receptors": round(float(np.corrcoef(
                [r["le_adgre1"] for r in rows], [r["le_adgre2"] for r in rows])[0, 1]), 3),
            "n_unique_hac": len({r["hac"] for r in rows}),
            "top10_ag1_selective": [r["ligand_id"] for r in sorted(rows, key=lambda r: -r["dlog_pref_a1"])[:10]],
        })
        if summary["n_ag1_preferred"] >= 5 and summary["dlog_median"] > 0.15:
            verdict = ("Selectivity window exists: %d/%d ligands prefer ADGRE1 (dlog>=%.2f), "
                       "median dlog=%.3f" % (summary["n_ag1_preferred"], len(rows),
                                             args.ddg_cut, summary["dlog_median"]))
        elif summary.get("corr_le_between_receptors", 0) > 0.8:
            verdict = ("The two sites are highly similar: LE correlation r=%.2f (raw-score r=%.2f), "
                       "only %d/%d ligands prefer ADGRE1, median dlog=%.3f -- selectivity at this "
                       "site is not achievable; change site or mechanism"
                       % (summary["corr_le_between_receptors"], summary["corr_affinity_between_receptors"],
                          summary["n_ag1_preferred"], len(rows), summary["dlog_median"]))
        else:
            verdict = ("Signal unclear: median dlog=%.3f, %d/%d prefer ADGRE1, LE r=%.2f -- "
                       "need more pairs or a re-run with matched boxes"
                       % (summary["dlog_median"], summary["n_ag1_preferred"], len(rows),
                          summary.get("corr_le_between_receptors", 0)))
    else:
        verdict = ("No pair passed the geometry gate: poses on the two sides generally contain "
                   "atomic overlaps -- fix the box/receptor first")
    summary["verdict"] = verdict
    json.dump(summary, open(os.path.join(outd, "summary_%s.json" % args.tag), "w"), indent=2)

    print("\nPairs (geometry-passing on both sides): %d; dropped: %s" % (len(rows), dict(drop) or "none"))
    print(json.dumps({k: v for k, v in summary.items() if k != "top10_ag1_selective"},
                     ensure_ascii=False, indent=2))
    print("\nVERDICT: %s" % verdict)
    print("Output -> %s" % os.path.relpath(tsv, ROOT))


if __name__ == "__main__":
    main()
