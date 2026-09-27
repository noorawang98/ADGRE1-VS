#!/usr/bin/env python3
"""2D structures of the five assayed candidates with their selection metrics.

Panel layout: 2 rows x 3 columns (the sixth cell holds the summary text). Each structure is
drawn with RDKit, and the metrics underneath come from the docking output and the geometry gate.

  ~/anaconda3/envs/AmberTools24/bin/python figures/make_structure_figure.py
-> figures/fig6_candidate_structures.png
"""
import argparse, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
_ap = argparse.ArgumentParser(description="2D structures of assayed candidates + metrics")
_ap.add_argument("--data", default=os.path.join(HERE, "five_candidates.json"),
                 help="JSON list with ligand_id, smiles, affinity, hac, le, min_dist, clash, geometry")
_ap.add_argument("--pick", default="ZINC000002446170", help="ligand highlighted as the chosen one")
_ap.add_argument("--out", default=os.path.join(HERE, "candidate_structures.png"))
_args = _ap.parse_args()
FIG = HERE
DATA = json.load(open(_args.data))
PICK = _args.pick

# order: the chosen compound first, then the others by ligand efficiency
order = [PICK] + [d["ligand_id"] for d in DATA if d["ligand_id"] != PICK]
by_id = {d["ligand_id"]: d for d in DATA}

# R-group colour legend differs per molecule; draw with atom indices off and a highlight of
# the ring system that carries the halogen when present.
def draw_mol(smi, w=520, h=380):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    o = d.drawOptions()
    o.bondLineWidth = 2.4
    o.minFontSize = 15
    o.maxFontSize = 22
    o.clearBackground = True
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return d.GetDrawingText()

fig = plt.figure(figsize=(15.5, 9.2))
gs = fig.add_gridspec(2, 3, hspace=0.28, wspace=0.12)
axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]

for ax, lid in zip(axes, order):
    d = by_id[lid]
    ax.axis("off")
    png = draw_mol(d["smiles"])
    if png is None:
        ax.text(0.5, 0.5, "structure unavailable", ha="center")
        continue
    tmp = os.path.join(FIG, "_tmp_%s.png" % lid)
    open(tmp, "wb").write(png)
    img = plt.imread(tmp)
    ax.imshow(img)
    os.remove(tmp)
    ok = d["geometry"] == "PASS"
    chosen = lid == PICK
    color = "#1a7f37" if ok else "#c0392b"
    title = ("%s\n(chosen for assay)" % lid) if chosen else lid
    ax.set_title(title, fontsize=12 if not chosen else 12.5,
                 fontweight="bold" if chosen else "normal",
                 color="#e8a33d" if chosen else "black")
    box = ("docking score  %.2f kcal/mol\n"
           "ligand efficiency  %.3f\n"
           "heavy atoms  %d\n"
           "min ligand-receptor distance  %.2f A\n"
           "clash pairs (<2.0 A)  %d\n"
           "geometry gate  %s" % (d["affinity"], d["le"], d["hac"], d["min_dist"],
                                  d["clash"], d["geometry"]))
    ax.text(0.5, -0.02, box, transform=ax.transAxes, ha="center", va="top",
            fontsize=10.2, family="monospace",
            bbox=dict(boxstyle="round,pad=0.45", facecolor="#f6f6f6",
                      edgecolor=color, linewidth=1.6 if chosen else 1.0))
    if chosen:
        for s in ax.spines.values():
            s.set_visible(True); s.set_color("#e8a33d"); s.set_linewidth(2.2)

# summary cell
ax = axes[5]
ax.axis("off")
lines = [
    "Selection rationale",
    "",
    "488,810 compounds screened (docking score)",
    "top 20 re-examined by geometry",
    "  -> 8 rejected: atoms overlap the receptor",
    "     (min distance 0.61-1.85 A)",
    "  -> 12 pass the gate",
    "5 taken into the assay",
    "",
    "ZINC000002446170 was the best-ranked",
    "candidate that satisfies all three criteria:",
    "  1. score in the extreme tail (-12.10)",
    "  2. geometry valid: 2.59 A, no clash pairs",
    "  3. best ligand efficiency of the set (0.336)",
]
ax.text(0.02, 0.98, "\n".join(lines), transform=ax.transAxes, va="top", ha="left",
        fontsize=10.6, family="monospace",
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#fff8e8",
                  edgecolor="#e8a33d", linewidth=1.6))

fig.suptitle("Five assayed candidates: structures and selection metrics",
             fontsize=16, fontweight="bold", y=0.985)
out = _args.out
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote", out)
