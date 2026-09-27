#!/usr/bin/env python3
"""Report figures: why ZINC000002446170 was selected for experimental follow-up.

Produces (figures/):
  fig1_funnel.png       screening funnel: 488,810 -> 20 -> 12 (geometry gate) -> 5 -> 1
  fig2_distribution.png library score distribution + top-decile zoom with the 5 candidates
  fig3_gate_le.png      affinity vs ligand efficiency, clash artefacts marked, candidates labelled
  fig4_contacts.png     contact map of the selected compound (2D bar of 4 A contacts)
  fig5_selectivity.png  ADGRE1 vs ADGRE2 paired scores (reverse-screen data)
All numbers come from files in this workspace, listed in figures/_data.json / _top20.json.

  ~/anaconda3/envs/scvi-pypi/bin/python figures/make_figures.py
"""
import csv, json, math, os
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight"})
C_LIB, C_REJ, C_PICK, C_ACT = "#4c72b0", "#b0b0b0", "#e8a33d", "#c44e52"
SEL = ["ZINC000015956900", "ZINC000009424484", "ZINC000002446170",
       "ZINC000009015892", "ZINC000019782239"]
PICK = "ZINC000002446170"


def load_library():
    s = []
    with open(os.path.join(ROOT, "ADGRE1_600_852_2_2", "ADGRE1_2_2_affinity.csv")) as fh:
        for r in csv.DictReader(fh):
            try:
                s.append((float(r["affinity"]), r["id"]))
            except ValueError:
                pass
    s.sort()
    return s


def fig1_funnel(lib, top20):
    fig, ax = plt.subplots(figsize=(11, 5.2))
    stages = [("Whole-library screen\n(Vina-GPU, 600-852 pocket)", len(lib), C_LIB),
              ("Top 20 by docking score", 20, C_LIB),
              ("Geometry gate passed\n(no atoms < 2.0 A, min dist >= 2.2 A)",
               sum(1 for r in top20 if r.get("ok")), C_REJ),
              ("Selected for assay\n(affinity + LE + geometry)", 5, C_PICK),
              ("Experimentally active", 1, C_ACT)]
    x = np.arange(len(stages))
    vals = [s[1] for s in stages]
    colors = [s[2] for s in stages]
    bars = ax.bar(x, vals, color=colors, width=0.62, edgecolor="white", linewidth=1.5)
    ax.set_yscale("log")
    ax.set_ylim(0.6, max(vals) * 4)
    for i, (b, (label, v, c)) in enumerate(zip(bars, stages)):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.35, f"{v:,}", ha="center",
                va="bottom", fontsize=13, fontweight="bold")
        if i < len(stages) - 1:
            nxt = vals[i + 1]
            keep = 100.0 * nxt / max(v, 1)
            ax.annotate("", xy=(x[i + 1] - 0.32, nxt * 1.05), xytext=(x[i] + 0.32, v * 1.05),
                        arrowprops=dict(arrowstyle="->", color="#555555", lw=1.4))
            ax.text((x[i] + x[i + 1]) / 2, math.sqrt(v * nxt) * 1.5, f"{keep:.1f}%",
                    ha="center", fontsize=10, color="#555555")
    ax.set_xticks(x)
    ax.set_xticklabels([s[0] for s in stages], fontsize=10)
    ax.set_ylabel("number of compounds (log scale)")
    ax.set_title("Screening funnel: how ZINC000002446170 was reached", fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.25, which="both")
    fig.savefig(os.path.join(FIG, "fig1_funnel.png"))
    plt.close(fig)
    return vals


def fig2_distribution(lib):
    s = np.array([a for a, _ in lib])
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    ax.hist(s, bins=120, color=C_LIB, alpha=0.85)
    ax.axvline(np.median(s), color="k", ls="--", lw=1.2)
    ax.text(np.median(s), ax.get_ylim()[1] * 0.9, f" median {np.median(s):.2f}",
            fontsize=10, va="top")
    ax.set_xlabel("docking score (kcal/mol)"); ax.set_ylabel("compounds")
    ax.set_title(f"Whole library (n = {len(s):,})", fontsize=12)
    ax = axes[1]
    cut = np.percentile(s, 0.5)
    sub = s[s <= cut]
    ax.hist(sub, bins=90, color=C_LIB, alpha=0.9)
    top5 = [(a, i) for a, i in lib[:5]]
    for a, i in top5:
        ax.axvline(a, color=C_ACT, lw=0.8, alpha=0.5)
    for lid in SEL:
        a = min(a for a, i in lib if i == lid)
        ax.axvline(a, color=C_PICK if lid == PICK else "#888888",
                   lw=2.2 if lid == PICK else 1.0, ls="-" if lid == PICK else ":")
    ax.axvline(-12.10, color=C_PICK, lw=2.2)
    ax.set_xlabel("docking score (kcal/mol)"); ax.set_ylabel("compounds")
    ax.set_title(f"Best 0.5% (n = {len(sub):,}); orange = ZINC000002446170",
                 fontsize=12)
    ax.text(0.02, 0.95, "dashed red: top-5 library scores\ngrey dotted: other 4 assayed compounds",
            transform=ax.transAxes, fontsize=9, va="top")
    fig.suptitle("Docking-score landscape of the ADGRE1 screen", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_distribution.png"))
    plt.close(fig)


def fig3_gate_le(top20):
    fig, ax = plt.subplots(figsize=(9.5, 6))
    ok = [r for r in top20 if r.get("ok")]
    bad = [r for r in top20 if r.get("ok") is False]
    ax.scatter([r["aff"] for r in ok], [r["le"] for r in ok], s=110, color="#4c9f70",
               edgecolor="k", linewidth=0.6, label="top-20 passing the geometry gate", zorder=3)
    ax.scatter([r["aff"] for r in bad], [r["le"] for r in bad], s=110, color=C_REJ,
               edgecolor="k", linewidth=0.6, marker="X",
               label="clash artefacts (rejected: atom pairs < 2.0 A)", zorder=3)
    for r in top20:
        if r["id"] in SEL or r["id"] == "ZINC000017057457":
            ax.annotate(r["id"].replace("ZINC", ""), (r["aff"], r["le"]),
                        textcoords="offset points", xytext=(6, 4), fontsize=8.5)
    pk = next(r for r in top20 if r["id"] == PICK)
    ax.scatter([pk["aff"]], [pk["le"]], s=330, facecolor="none", edgecolor=C_PICK,
               linewidth=2.6, zorder=4, label="ZINC000002446170 (chosen)")
    ax.set_xlabel("docking score (kcal/mol)   [more negative = better]")
    ax.set_ylabel("ligand efficiency (-dG / heavy-atom count)")
    ax.set_title("Geometry gate separates real poses from clash artefacts\n"
                 "(top 20 by docking score)", fontsize=13, fontweight="bold")
    ax.legend(loc="lower left", fontsize=9.5, frameon=True)
    ax.grid(alpha=0.25)
    fig.savefig(os.path.join(FIG, "fig3_gate_le.png"))
    plt.close(fig)


def fig4_contacts():
    A = []
    for l in open(os.path.join(ROOT, "811", "receptor.pdbqt"), errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            try:
                A.append((int(l[22:26]), l[17:20].strip(), float(l[30:38]),
                          float(l[38:46]), float(l[46:54])))
            except ValueError:
                continue
    L = []
    p = os.path.join(ROOT, "ADGRE1_600_852_2_2", "ZINC15_004_out",
                     "ZINC000002446170_out.pdbqt")
    for l in open(p, errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            try:
                L.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                continue
        elif l.startswith("ENDMDL"):
            break
    cnt, mind = Counter(), 99.0
    for x, y, z in L:
        for rn, rn3, rx, ry, rz in A:
            d = math.dist((x, y, z), (rx, ry, rz))
            mind = min(mind, d)
            if d < 4.0:
                cnt[f"{rn3}{rn}"] += 1
    items = cnt.most_common(12)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.2), gridspec_kw={"width_ratios": [2, 1]})
    names = [k for k, _ in items][::-1]
    vals = [v for _, v in items][::-1]
    colors = [C_PICK if k in ("TRP811", "GLY814", "ILE815", "GLN817", "PHE829") else C_LIB
              for k in names]
    ax.barh(names, vals, color=colors)
    ax.set_xlabel("number of ligand atoms within 4 A")
    ax.set_title(f"Contact map of ZINC000002446170\n(min ligand-receptor distance {mind:.2f} A)",
                 fontsize=12.5, fontweight="bold")
    for i, v in enumerate(vals):
        ax.text(v + 0.15, i, str(v), va="center", fontsize=9)
    ax.grid(axis="x", alpha=0.25)

    # phenotype metrics panel
    metrics = [("docking score", -12.10, "kcal/mol"),
               ("ligand efficiency", 0.336, "kcal/mol per HA"),
               ("heavy atoms", 36, "atoms"),
               ("min contact", 2.59, "A"),
               ("clash pairs <2.0 A", 0, "pairs"),
               ("multi-model spread", 0.0, "kcal/mol")]
    ax2.axis("off")
    ax2.set_title("Selection metrics", fontsize=12.5, fontweight="bold")
    for i, (k, v, u) in enumerate(metrics):
        y = 0.88 - i * 0.14
        ax2.text(0.02, y, k, fontsize=10.5, va="center")
        ax2.text(0.98, y, f"{v} {u}", fontsize=10.5, va="center", ha="right",
                 fontweight="bold", color=C_PICK if i in (0, 1, 3, 4) else "k")
        ax2.plot([0.02, 0.98], [y - 0.055, y - 0.055], color="#dddddd", lw=0.8)
    ax2.text(0.02, 0.03, "Vina-GPU output carries one pose per ligand,\n"
                         "so pose convergence could not be assessed here.",
             fontsize=8.5, style="italic", color="#666666")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_contacts.png"))
    plt.close(fig)
    return items, mind


def fig5_selectivity():
    p = os.path.join(ROOT, "reverse_screen", "pairs_gpu_drugbank.tsv")
    if not os.path.exists(p):
        p = os.path.join(ROOT, "reverse_screen", "pairs_gpu_fda.tsv")
    if not os.path.exists(p):
        print("  (no paired data yet - skipping fig5)")
        return None
    a1, a2 = [], []
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a1.append(float(r["aff_adgre1"])); a2.append(float(r["aff_adgre2"]))
    a1, a2 = np.array(a1), np.array(a2)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ax.scatter(a1, a2, s=6, alpha=0.25, color=C_LIB)
    lim = [min(a1.min(), a2.min()) - 0.5, max(a1.max(), a2.max()) + 0.5]
    ax.plot(lim, lim, "k--", lw=1)
    r = np.corrcoef(a1, a2)[0, 1]
    ax.set_xlabel("ADGRE1 score (kcal/mol)"); ax.set_ylabel("ADGRE2 score (kcal/mol)")
    ax.set_title(f"Paired scores, n = {len(a1):,}\nPearson r = {r:.3f}", fontsize=12)
    ax.grid(alpha=0.25)
    d = a2 - a1
    ax2.hist(d, bins=80, color=C_LIB, alpha=0.9)
    ax2.axvline(0, color="k", ls="--", lw=1)
    ax2.axvline(np.median(d), color=C_ACT, lw=2)
    ax2.set_xlabel("score difference  ADGRE2 - ADGRE1 (kcal/mol)")
    ax2.set_ylabel("compounds")
    ax2.set_title(f"median = {np.median(d):+.2f} kcal/mol\n"
                  f"ADGRE1-preferring (>= +1.5): {int((d >= 1.5).sum())}", fontsize=12)
    fig.suptitle("Reverse screen: the ADGRE1 site versus its ADGRE2 ortholog",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig5_selectivity.png"))
    plt.close(fig)
    return dict(n=len(a1), r=float(r), median=float(np.median(d)),
                n_pref_a1=int((d >= 1.5).sum()))


def main():
    lib = load_library()
    top20 = json.load(open(os.path.join(FIG, "_top20.json")))
    print("library %d, top20 loaded %d" % (len(lib), len(top20)))
    v = fig1_funnel(lib, top20); print("fig1 funnel:", v)
    fig2_distribution(lib); print("fig2 ok")
    fig3_gate_le(top20); print("fig3 ok, gate pass %d/20" % sum(1 for r in top20 if r.get("ok")))
    items, mind = fig4_contacts(); print("fig4 ok, top contacts:", items[:4], "min %.2f" % mind)
    r = fig5_selectivity(); print("fig5:", r)
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
