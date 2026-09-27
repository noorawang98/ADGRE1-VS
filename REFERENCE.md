# Reference values used in this toolkit

These constants were measured in the campaign and are the starting points the scripts assume.
They are given here so that the thresholds can be reviewed rather than taken on faith.

## 1. Geometry gate

| Constant | Value | Where used |
|---|---|---|
| heavy-atom collision threshold | **2.0 A** | `workflow/select_poses.py`, `workflow/pocket_assign.py`, `workflow/reverse_screen_analysis.py` |
| minimum acceptable ligand-receptor distance | **2.2 A** | same |
| contact definition | **4.0 A** | contact sets and residue lining |
| buried-atom criterion | >= 5 receptor atoms within 5.0 A | surface-site triage |

## 2. Site definition (ADGRE1, 7TM region)

Docked ligands contact residues **690-710** (outer surface of TM3). Residue names below are the
correct `resName` values (an earlier note in the campaign mis-parsed the PDB columns).

| Item | Value |
|---|---|
| occupied-site centroid | **(-54.3, -10.5, 54.2)** |
| contact residues (4 A, selected compound) | GLY814, PHE652, ILE815, PHE829, ARG741, GLN817, TRP811, TRP743, PHE676, LEU755, THR746 |
| transmembrane segments (19-residue hydrophobicity window > 1.0) | centres 609, 647, 683, 723, 762, 823 |
| pocket character | solvent-exposed surface groove; no enclosed cavity larger than 40 A^3 |

## 3. Ortholog mapping to the closest paralogue site

| Item | Value |
|---|---|
| alignment | full-length BLOSUM62 global alignment, gap -10/-0.5 |
| identity over the 7TM region | 0.545 (full length 0.356) |
| site-residue conservation | 13/21 identical (0.62) for the measured patch; 12/16 (0.75) for the compound site |
| orthologous residue range | patch site 572-593; compound site 488-715 |
| ortholog box centre | patch site (-50.7, -5.2, 54.9); compound site (-27.5, 0.1, 31.3) |

Rule enforced by `workflow/reverse_screen_analysis.py`: results are expressed as
`dlog = log10|aff_B| - log10|aff_A|` and rank statistics, never as a difference of dG values,
because the two pockets are not the same pocket.

## 4. Ranking key

`LE = -dG / heavy-atom count` (heavy atoms counted from the pose file, so no bond-order
information is required). Raw docking scores grow with molecular size; ligand efficiency removes
that bias and is therefore the primary key, with score used only as a tie-breaker.

## 5. GPU screen settings that were actually used

| Setting | Value | Rationale |
|---|---|---|
| box | 24 A cube, centred on the occupied site | a blind 80x60x140 A box cannot run on Vina-GPU and wastes search space |
| `thread` | 3000 | 8000 is oversaturated on a 10 GB device |
| `search_depth` | 10 | halves sampling with no loss for reverse screening |
| concurrency | one Vina-GPU instance at a time | two instances on one card measured 5x slower |
| pilot | 200 ligands before any full run | converts throughput from a guess into a measurement |

## 6. Environment findings worth reusing

* `repo.anaconda.com` and `conda.anaconda.org` were unreachable from the working host while the
  Tsinghua mirror served both; `defaults` in `~/.condarc` therefore had to be dropped.
* `~/.cache`, site-packages and `$HOME/anaconda3` were read-only in the container, so conda
  needed `--prefix` into a writable path plus a custom `CONDA_PKGS_DIRS`.
* `uv venv` does not seed pip and uv-python builds often omit `ensurepip`; install with
  `uv pip install --python <env>/bin/python ...` instead of relying on pip inside the env.
