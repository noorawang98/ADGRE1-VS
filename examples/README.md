# Runnable example

Everything here is small enough to commit, so the geometry gate can be exercised without any
docking run of your own.

```
examples/
  receptor.pdb               ADGRE1 600-852 model (6818 heavy atoms, same frame as the poses)
  pose_files/*.pdb           six docked poses for six example ligands, one per file
  scores_6.tsv               id <TAB> affinity table for those six ligands
  selection_example.tsv      produced by the command below
```

## Run it

```bash
python workflow/select_poses.py \
  --table      examples/scores_6.tsv \
  --pose-dirs  examples/pose_files \
  --receptor   examples/receptor.pdb \
  --min-hac    20 \
  --out        examples/selection_example.tsv
```

Expected output (this is what the committed `selection_example.tsv` contains):

| ligand_id | affinity | LE | heavy atoms | min distance | clash pairs | verdict |
|---|---|---|---|---|---|---|
| 102697 | -10.50 | 0.420 | 25 | 3.30 A | 0 | pass |
| 10178705 | -8.61 | 0.297 | 29 | 2.86 A | 0 | pass |
| 10156 | -8.83 | 0.294 | 30 | 3.22 A | 0 | pass |
| 10026128 | -9.45 | 0.278 | 34 | 2.85 A | 0 | pass |
| 10024590 | -8.18 | 0.273 | 30 | 2.82 A | 0 | pass |
| 10298385 | -8.87 | 0.246 | 36 | 3.17 A | 0 | pass |

Note the ordering: the compound with the best docking score (-10.50) also has the best ligand
efficiency here, but on the full campaign that is not the general rule, which is why both
columns are reported.

## Larger example set

`data/example_200_ligands/` holds the 200 ligands used for the GPU pilot, with three tables:

| File | Columns |
|---|---|
| `ligands.tsv` | `ligand_id`, `smiles`, `heavy_atoms` |
| `scores.tsv` | `id`, `affinity` |
| `scored_with_LE.tsv` | `ligand_id`, `smiles`, `docking_score`, `heavy_atoms`, `ligand_efficiency` |
| `scored_with_geometry.tsv` | the above plus `min_distance`, `clash_pairs`, `geometry` |

On those 200 ligands the geometry gate rejects **38 (19 %)** as clash artefacts, and 162 pass;
the best passing compound scores -10.50 kcal/mol at 3.30 A minimum distance. The SMILES in these
tables were taken from the `REMARK SMILES` records written by Vina-GPU itself, so they carry
correct bond orders (see `data/README.md` for why that matters).
