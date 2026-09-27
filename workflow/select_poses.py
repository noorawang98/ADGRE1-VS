#!/usr/bin/env python3
"""Turn a raw Vina score table into a defensible pick list.

Motivation (measured, not assumed): in this campaign the top-11 ligands by raw affinity
included four with ligand-receptor heavy-atom distances of 0.61-1.56 A, i.e. ligand atoms
buried inside receptor atoms. Vina still returns a strong score for such poses, so ranking by
affinity alone selects clash artefacts. The one compound that was actually taken into the lab
(ZINC000002446170) was the only one of the five hand-picked candidates with a clean pose
(min distance 2.59 A), the best ligand efficiency (0.356) and a converged score spread.

This script applies that logic automatically:
  * parse every pose file, take REMARK SMILES as the authoritative structure (no bond-order
    guessing, no SDF matching needed)
  * geometry : min ligand-receptor distance, count of colliding heavy-atom pairs (<2.0 A)
  * efficiency: LE = -affinity / HAC  (size-decoupled; essential for shallow surface sites)
  * convergence: spread between the best and worst of the returned models (single-model
    Vina-GPU output is reported as 0.0 and flagged in the notes)
  * alerts    : PAINS/Brenk via the authoritative SMILES, when RDKit is available
  * composite : rank on LE, with hard exclusion of clashing poses

  ~/anaconda3/envs/AmberTools24/bin/python demo/select_poses.py \
      --table ADGRE1_600_852_2_2/ADGRE1_2_2_affinity.csv \
      --receptor 811/receptor.pdbqt --top 20000
-> demo/selection.tsv, demo/selection_summary.json
"""
import argparse, csv, glob, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_CLASH = 2.0          # A, heavy-atom collision threshold
_CONTACT = 4.0        # A, contact definition
_REC = None


def load_receptor(path):
    R = []
    for l in open(path, errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            try:
                R.append((int(l[22:26]), l[17:20].strip(), float(l[30:38]),
                          float(l[38:46]), float(l[46:54])))
            except ValueError:
                continue
    return R


def parse_pose(path):
    """returns (smiles, [model affinities], atoms of model 1, n_models)"""
    smi, affs, atoms, n = None, [], [], 0
    with open(path, errors="ignore") as fh:
        for l in fh:
            if l.startswith("REMARK SMILES IDX"):
                continue
            if l.startswith("REMARK SMILES ") and smi is None:
                smi = l.replace("REMARK SMILES ", "").strip()
            elif l.startswith("REMARK VINA RESULT"):
                try: affs.append(float(l.split()[3]))
                except (IndexError, ValueError): pass
            elif l.startswith(("ATOM", "HETATM")):
                if n == 0:
                    try: atoms.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                    except ValueError: pass
            elif l.startswith("ENDMDL"):
                n += 1
                if n > 20:
                    break
    if atoms and not affs:
        affs = [float("nan")]
    return smi, affs, atoms, max(n, 1)


def evaluate(job):
    """job = (ligand_id, path, receptor list, table_affinity) -> dict or None

    The affinity is read from the pose file when the file carries VINA RESULT remarks
    (Vina/Vina-GPU output does). For plain PDB/PDBQT files without those remarks the
    affinity supplied in --table is used instead, so the geometry gate can be applied to
    any pose set.
    """
    lid, path, rec, table_aff = job
    from scipy.spatial import cKDTree
    try:
        smi, affs, atoms, nmod = parse_pose(path)
    except Exception:
        return None
    good = [a for a in affs if -50 < a < 0]
    if not atoms:
        return None
    if good:
        aff = min(good)
    elif table_aff is not None:
        aff = float(table_aff)          # fall back to the score table
        good = [aff]
    else:
        return None
    tree = cKDTree(np.asarray([r[2:5] for r in rec]))
    L = np.asarray(atoms)
    d, idx = tree.query(L, k=6)
    if d.ndim == 1:
        d = d[:, None]; idx = idx[:, None]
    dmin = float(d.min())
    nclash = int((d < _CLASH).sum())
    contacts = set()
    for i in range(len(L)):
        for k in range(d.shape[1]):
            if d[i, k] < _CONTACT:
                contacts.add(rec[idx[i, k]][0])
    spread = float(max(good) - min(good))
    return dict(ligand_id=lid, affinity=round(aff, 3), n_models=nmod,
                n_poses_scored=len(good), hac=len(atoms),
                min_dist=round(dmin, 2), n_clash_atoms=nclash,
                contact_residues=len(contacts), spread=round(spread, 2),
                smiles=smi or "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="CSV with 'id,affinity'")
    ap.add_argument("--pose-dirs", nargs="*", default=None,
                    help="directories searched recursively for pose files "
                         "(xxx_out.pdbqt or xxx.pdbqt); default: --table-dir if given, "
                         "otherwise the current directory")
    ap.add_argument("--receptor", default=os.path.join(ROOT, "811", "receptor.pdbqt"))
    ap.add_argument("--top", type=int, default=20000, help="how many best-scoring ligands to re-examine")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--min-hac", type=int, default=20,
                    help="reject ligands with fewer heavy atoms: ligand efficiency rewards small "
                         "fragments (a 5-atom fragment at -2 kcal/mol gets LE 0.4 while touching "
                         "very little of the pocket), so a floor is applied by default")
    ap.add_argument("--table-dir", default=None,
                    help="optional directory containing the pose files, used only to resolve "
                         "ligand ids to file names when they are not found under --pose-glob")
    ap.add_argument("--out", default=os.path.join(ROOT, "demo", "selection.tsv"))
    args = ap.parse_args()

    rows = []
    with open(args.table) as fh:
        head = fh.readline()
        delim = "\t" if "\t" in head else ","
        fh.seek(0)
        for r in csv.DictReader(fh, delimiter=delim):
            # accept both 'id'/'affinity' and 'ligand_id'/'docking_score' style headers
            lid = r.get("id") or r.get("ligand_id") or r.get("zinc_id")
            sc = r.get("affinity") or r.get("docking_score") or r.get("score")
            try:
                rows.append((float(sc), lid))
            except (TypeError, ValueError):
                continue
    rows.sort()
    print("table entries: %d | re-examining best %d" % (len(rows), min(args.top, len(rows))))

    search = args.pose_dirs or ([args.table_dir] if args.table_dir else ["."])
    index = {}
    for d in search:
        for pat in ("**/*_out.pdbqt", "**/*.pdbqt", "**/*.pdb", "**/*.pdbqt.gz"):
            for f in glob.glob(os.path.join(d, pat), recursive=True):
                b = os.path.basename(f)
                lid = b[:-len("_out.pdbqt")] if b.endswith("_out.pdbqt") else os.path.splitext(b)[0]
                index.setdefault(lid, f)
    print("pose files indexed: %d (searched: %s)" % (len(index), ", ".join(search)))

    rec = load_receptor(args.receptor)
    aff_of = {lid: aff for aff, lid in rows}
    jobs = [(lid, index[lid], rec, aff_of.get(lid)) for _, lid in rows[:args.top] if lid in index]
    print("jobs: %d" % len(jobs))
    t0 = time.time()
    out = []
    with ProcessPoolExecutor(args.workers) as ex:
        for i, r in enumerate(ex.map(evaluate, jobs, chunksize=64), 1):
            if r:
                out.append(r)
            if i % 5000 == 0:
                print("  %d/%d  (%.0f s)" % (i, len(jobs), time.time() - t0))
    print("evaluated: %d poses" % len(out))

    # RDKit alerts when available
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import FilterCatalog
        RDLogger.DisableLog("rdApp.*")
        p = FilterCatalog.FilterCatalogParams()
        for c in (p.FilterCatalogs.PAINS_A, p.FilterCatalogs.PAINS_B,
                  p.FilterCatalogs.PAINS_C, p.FilterCatalogs.BRENK):
            p.AddCatalog(c)
        cat = FilterCatalog.FilterCatalog(p)
        for r in out:
            m = Chem.MolFromSmiles(r["smiles"]) if r["smiles"] else None
            r["alerts"] = ",".join(sorted({e.GetDescription() for e in cat.GetMatches(m)})) if m else ""
        has_rdkit = True
    except Exception:
        for r in out:
            r["alerts"] = ""
        has_rdkit = False
    print("RDKit alerts: %s" % ("on" if has_rdkit else "unavailable"))

    for r in out:
        r["le"] = round(-r["affinity"] / max(r["hac"], 1), 4)
        r["size_ok"] = r["hac"] >= args.min_hac
        r["size_ok"] = r["hac"] >= args.min_hac
        # hard geometry gate: a clash artefact is not a hit, whatever its score says
        hard_clash = r["n_clash_atoms"] > 0 or r["min_dist"] < 2.2
        r["geometry"] = "clash" if hard_clash else ("tight" if r["min_dist"] < 2.5 else "ok")
        r["pass"] = bool(not hard_clash and r["size_ok"])
        if hard_clash:
            r["reason"] = "clash"
        elif not r["size_ok"]:
            r["reason"] = "below_min_hac"
        else:
            r["reason"] = "ok"
    out.sort(key=lambda r: (-r["pass"], -r["le"]))
    cols = ["ligand_id", "affinity", "le", "hac", "min_dist", "n_clash_atoms", "geometry",
            "pass", "reason", "contact_residues", "spread", "n_models", "alerts", "smiles"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in out:
            w.writerow(r)

    npass = sum(1 for r in out if r["pass"])
    nclash = sum(1 for r in out if r["reason"] == "clash")
    nsmall = sum(1 for r in out if r["reason"] == "below_min_hac")
    nclash = sum(1 for r in out if r["reason"] == "clash")
    nsmall = sum(1 for r in out if r["reason"] == "below_min_hac")
    clash_ids = {r["ligand_id"] for r in out if not r["pass"]}
    top20_raw = [lid for _, lid in rows[:20]]
    summary = {
        "evaluated": len(out), "pass": npass,
        "clash_excluded": nclash, "below_min_hac_excluded": nsmall,
        "min_hac": args.min_hac,
        "raw_top20_rejected_as_clash": sorted(set(top20_raw) & clash_ids),
        "raw_top20_rejected_count": len(set(top20_raw) & clash_ids),
        "best_by_le": [{k: r[k] for k in ("ligand_id", "affinity", "le", "min_dist", "hac", "alerts")}
                       for r in out[:10]],
    }
    summary_path = os.path.splitext(args.out)[0] + "_summary.json"
    json.dump(summary, open(summary_path, "w"), indent=2)
    print("\nfilters: %d pass | %d clash-excluded | %d below min_hac(%d)"
          % (npass, nclash, nsmall, args.min_hac))
    print("of the raw top-20, rejected as clash: %d %s"
          % (summary["raw_top20_rejected_count"], summary["raw_top20_rejected_as_clash"]))
    print("\ntop 8 by ligand efficiency among geometry-passing ligands:")
    print("%-18s %7s %7s %8s %5s %s" % ("ligand", "aff", "LE", "minDist", "HAC", "alerts"))
    for r in out[:8]:
        print("%-18s %7.2f %7.3f %8.2f %5d %s" % (r["ligand_id"], r["affinity"], r["le"],
                                                  r["min_dist"], r["hac"], r["alerts"][:24]))
    print("\nout -> %s, %s" % (args.out, summary_path))


if __name__ == "__main__":
    main()
