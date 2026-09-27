#!/usr/bin/env python3
"""Residue-level ortholog mapping between the ADGRE1 and ADGRE2 7TM pockets.

WHY: the two pockets are NOT the same pocket (different composition, volume, solvent
exposure), so receptor-to-receptor dG subtraction is meaningless. To judge whether
selectivity is achievable at all, the 811 site (ligands actually bind ADGRE1 690-710,
outer surface of TM3) has to be mapped onto the ORTHOLOGOUS ADGRE2 residues.

METHOD: rebuild sequences from both receptor pdbqt files -> keep the 7TM region (located
from hydrophobic segments) -> Needleman-Wunsch global alignment (BLOSUM62, gap -10/-0.5)
-> emit a per-residue ortholog table for the 811 site.
NOTE: both receptors here are PREDICTED models whose numbering does not match UniProt
full length, so the mapping is valid only between these two models.

  ~/anaconda3/envs/AmberTools24/bin/python align_pocket_orthologs.py
-> reverse_screen/pocket_alignment.tsv, pocket_orthologs.json
"""
import argparse, json, os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
HYD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,
       'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
# BLOSUM62 (standard matrix)
AA = "ARNDCQEGHILKMFPSTWYV"
BLOSUM = {}
_rows = """A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4"""
for line in _rows.strip().split("\n"):
    p = line.split()
    for j, v in enumerate(p[1:]):
        BLOSUM[(p[0], AA[j])] = int(v)


def read_seq(pdbqt):
    d = {}
    for l in open(pdbqt, errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            try:
                d.setdefault(int(l[22:26]), THREE2ONE.get(l[17:20].strip(), "X"))
            except ValueError:
                continue
    ks = sorted(d)
    return "".join(d[k] for k in ks), ks[0], ks[-1]


def tm_bounds(seq, start, win=19, thr=1.0):
    idx = list(range(start, start + len(seq)))
    sc = {}
    for i in range(len(seq) - win):
        sc[idx[i]] = sum(HYD.get(seq[i + j], 0) for j in range(win)) / win
    tm = [k for k, v in sc.items() if v > thr]
    g = []
    for k in sorted(tm):
        if g and k - g[-1][-1] <= 8:
            g[-1].append(k)
        else:
            g.append([k])
    return [(x[0], x[-1] + win) for x in g]


def nw(s1, s2, gap_open=-10, gap_ext=-0.5):
    n, m = len(s1), len(s2)
    NEG = -1e9
    M = [[NEG]*(m+1) for _ in range(n+1)]
    X = [[NEG]*(m+1) for _ in range(n+1)]
    Y = [[NEG]*(m+1) for _ in range(n+1)]
    M[0][0] = 0
    for i in range(1, n+1):
        X[i][0] = gap_open + (i-1)*gap_ext
    for j in range(1, m+1):
        Y[0][j] = gap_open + (j-1)*gap_ext
    for i in range(1, n+1):
        for j in range(1, m+1):
            s = BLOSUM.get((s1[i-1], s2[j-1]), -4)
            M[i][j] = max(M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1]) + s
            X[i][j] = max(M[i-1][j] + gap_open, X[i-1][j] + gap_ext)
            Y[i][j] = max(M[i][j-1] + gap_open, Y[i][j-1] + gap_ext)
    # traceback
    a1, a2 = [], []
    st = max((M[n][m], 'M'), (X[n][m], 'X'), (Y[n][m], 'Y'))
    i, j, state = n, m, st[1]
    score = st[0]
    while i > 0 or j > 0:
        if state == 'M':
            a1.append(s1[i-1]); a2.append(s2[j-1])
            prev = max((M[i-1][j-1], 'M'), (X[i-1][j-1], 'X'), (Y[i-1][j-1], 'Y'))[1]
            i -= 1; j -= 1; state = prev
        elif state == 'X':
            a1.append(s1[i-1]); a2.append('-')
            state = 'X' if X[i][j] == X[i-1][j] + gap_ext else 'M'
            i -= 1
        else:
            a1.append('-'); a2.append(s2[j-1])
            state = 'Y' if Y[i][j] == Y[i][j-1] + gap_ext else 'M'
            j -= 1
    return "".join(reversed(a1)), "".join(reversed(a2)), score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a1", default=os.path.join(ROOT, "811", "receptor.pdbqt"))
    ap.add_argument("--a2", default="/home/workspace/MD/MD/HDY/ADGRE2/receptor.pdbqt")
    ap.add_argument("--site-lo", type=int, default=690)
    ap.add_argument("--site-hi", type=int, default=710,
                    help="ADGRE1 residues actually contacted by docked ligands (TM3 outer surface)")
    ap.add_argument("--pad", type=int, default=40, help="extra residues taken on each side for the alignment")
    ap.add_argument("--outdir", default=os.path.join(ROOT, "reverse_screen"))
    args = ap.parse_args()

    s1, a1lo, a1hi = read_seq(args.a1)
    s2, a2lo, a2hi = read_seq(args.a2)
    print("ADGRE1 model: %d-%d (%d aa)" % (a1lo, a1hi, len(s1)))
    print("ADGRE2 model: %d-%d (%d aa)" % (a2lo, a2hi, len(s2)))
    tm1 = tm_bounds(s1, a1lo); tm2 = tm_bounds(s2, a2lo)
    print("ADGRE1 hydrophobic segments: %s" % tm1)
    print("ADGRE2 hydrophobic segments: %s" % tm2)

    # take the ADGRE1 site +/- pad and the corresponding ADGRE2 7TM region (C-term pairing)
    lo1, hi1 = args.site_lo - args.pad, args.site_hi + args.pad
    sub1 = "".join(s1[k - a1lo] for k in range(lo1, hi1 + 1))
    # ADGRE2 side: pair ADGRE1 TM segments to ADGRE2 in C-terminal order
    t1 = [t for t in tm1 if t[1] < 860]
    t2 = [t for t in tm2 if t[0] > 300]
    n = min(len(t1), len(t2))
    pair = list(zip(t1[-n:], t2[-n:]))
    print("TM segment pairs (from C-terminus): %s" % [(a, b) for a, b in pair])
    # locate the TM pair covering the ADGRE1 site, use it to bound the ADGRE2 range
    cov = [b for a, b in pair if a[0] - 5 <= args.site_lo <= a[1] + 5 or a[0] - 5 <= args.site_hi <= a[1] + 5]
    if not cov:
        cov = [pair[-1][1]]
    lo2, hi2 = min(b[0] for b in cov) - args.pad, max(b[1] for b in cov) + args.pad
    sub2 = "".join(s2[k - a2lo] for k in range(lo2, hi2 + 1))
    print("Alignment window: ADGRE1 %d-%d vs ADGRE2 %d-%d" % (lo1, hi1, lo2, hi2))

    A1, A2, score = nw(sub1, sub2)
    ident = sum(1 for x, y in zip(A1, A2) if x == y and x != "-") / max(1, sum(1 for x in A1 if x != "-"))
    print("Alignment score %d, identity %.3f (length %d)" % (score, ident, len(A1)))
    print("ADGRE1 %s" % A1[:120])
    print("ADGRE2 %s" % A2[:120])

    # per-residue mapping
    maps, i1, i2 = [], lo1, lo2
    for x, y in zip(A1, A2):
        if x != "-" and y != "-":
            maps.append((i1, x, i2, y))
        if x != "-":
            i1 += 1
        if y != "-":
            i2 += 1
    os.makedirs(args.outdir, exist_ok=True)
    tsv = os.path.join(args.outdir, "pocket_alignment.tsv")
    with open(tsv, "w") as fh:
        fh.write("adgre1_res\tadgre1_aa\tadgre2_res\tadgre2_aa\tidentical\tin_site\n")
        for r1, x, r2, y in maps:
            fh.write("%d\t%s\t%d\t%s\t%d\t%d\n" % (r1, x, r2, y, int(x == y),
                                                   int(args.site_lo <= r1 <= args.site_hi)))
    site_map = [(r1, x, r2, y) for r1, x, r2, y in maps if args.site_lo <= r1 <= args.site_hi]
    print("\n=== Per-residue orthologs for the 811 site (%d-%d) ===" % (args.site_lo, args.site_hi))
    for r1, x, r2, y in site_map:
        print("   ADGRE1 %s%-4d  <->  ADGRE2 %s%-4d  %s" % (x, r1, y, r2, "same" if x == y else "diff"))
    json.dump({"adgre1_model": [a1lo, a1hi], "adgre2_model": [a2lo, a2hi],
               "aln_region": {"adgre1": [lo1, hi1], "adgre2": [lo2, hi2]},
               "identity": round(ident, 3), "score": score,
               "tm_pairs": [[list(a), list(b)] for a, b in pair],
               "site_map": [{"a1": r1, "aa1": x, "a2": r2, "aa2": y, "same": x == y}
                            for r1, x, r2, y in site_map],
               "note": "both receptors are predicted models whose numbering differs from UniProt; this mapping is valid only between these two models"},
              open(os.path.join(args.outdir, "pocket_orthologs.json"), "w"), indent=2)
    print("\nOutput -> %s / pocket_orthologs.json" % os.path.relpath(tsv, ROOT))


if __name__ == "__main__":
    main()
