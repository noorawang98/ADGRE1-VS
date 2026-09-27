#!/usr/bin/env python3
"""Report how much of the local PDBQT library the downloaded ZINC20 SMILES actually covers.

Reads zinc20_smiles/by_tranche/*.tsv (output of fetch_zinc20_smiles.sh) and
<root>/filters/zinc_ids_bare.txt (local library ids), then prints per-directory counts,
overall coverage, and the size/polarity profile of the matched subset versus the whole
library. This is the go/no-go number for "is downloading all 121 tranche dirs worth it".

  ~/anaconda3/envs/adfr/bin/python repo/scripts/coverage_report.py \
      --dir /home/workspace/MD/MD/HDY/ADGRE1/filters/z20_bench \
      --lib /home/workspace/MD/MD/HDY/ADGRE1/filters/zinc_ids_bare.txt
"""
import argparse, glob, os, statistics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="dir containing by_tranche/*.tsv")
    ap.add_argument("--lib", required=True, help="file with one bare library id per line")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    lib = set(l.strip() for l in open(args.lib) if l.strip())
    files = sorted(glob.glob(os.path.join(args.dir, "by_tranche", "*.tsv")))
    if not files:
        raise SystemExit("no by_tranche/*.tsv under %s" % args.dir)
    print("library ids: %d | tranche dirs downloaded: %d" % (len(lib), len(files)))

    matched = set()
    mw, logp = [], []
    per_dir = []
    for f in files:
        n = 0
        for line in open(f, errors="ignore"):
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[1] == "zinc_id":
                continue
            z = p[1].lstrip("ZINC").lstrip("0") or "0"
            if z in lib:
                n += 1
                matched.add(z)
                if len(p) > 3:
                    try: mw.append(float(p[3]))
                    except ValueError: pass
                if len(p) > 4:
                    try: logp.append(float(p[4]))
                    except ValueError: pass
        per_dir.append((os.path.basename(f), n))
        print("  %-12s matched=%6d" % (os.path.basename(f), n))

    cov = 100.0 * len(matched) / max(len(lib), 1)
    print("\nTOTAL matched %d / %d = %.3f%%" % (len(matched), len(lib), cov))
    if mw:
        print("matched subset: MW median %.1f (p10 %.1f, p90 %.1f) | cLogP median %.2f (p10 %.2f, p90 %.2f)"
              % (statistics.median(mw), sorted(mw)[len(mw)//10], sorted(mw)[len(mw)*9//10],
                 statistics.median(logp), sorted(logp)[len(logp)//10], sorted(logp)[len(logp)*9//10]))
    # extrapolation to all 121 dirs, using the per-dir mean (AA/AB/AC/BA/BB are all low-MW dirs,
    # so treat this as an UPPER bound: later dirs hold larger, more polar molecules)
    mean = sum(n for _, n in per_dir) / len(per_dir)
    print("per-dir mean %.0f -> 121 dirs ~ %.0f ids (%.1f%% of library, upper bound)"
          % (mean, mean * 121, 100 * mean * 121 / max(len(lib), 1)))
    print("decision guide: if the extrapolated coverage stays below ~25%%, do NOT download "
          "all 121 dirs; use ligand_qc.py (heavy-atom count + connectivity alerts, no SMILES "
          "needed) for whole-library screening instead.")
    if args.out:
        with open(args.out, "w") as fh:
            fh.write("dir\tmatched\n")
            for d, n in per_dir:
                fh.write("%s\t%d\n" % (d, n))
            fh.write("TOTAL\t%d\n" % len(matched))


if __name__ == "__main__":
    main()
