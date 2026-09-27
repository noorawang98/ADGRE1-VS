# Data

No ligand structures, docking outputs or experimental measurements are shipped in this
repository: the docking sets used in the campaign are several gigabytes and the activity data
belong to the experimental collaborators.

## Included example data

| Path | Content |
|---|---|
| `example_200_ligands/ligands.tsv` | the 200 ligands used for the GPU pilot: `ligand_id`, `smiles`, `heavy_atoms` |
| `example_200_ligands/scores.tsv` | `id`, `affinity` (what `workflow/select_poses.py --table` reads) |
| `example_200_ligands/scored_with_LE.tsv` | adds docking score and ligand efficiency, sorted by LE |
| `example_200_ligands/scored_with_geometry.tsv` | adds `min_distance`, `clash_pairs`, `geometry` |

Result on this set: **38 of 200 (19 %) rejected as clash artefacts**, 162 pass the geometry gate.
SMILES were extracted from the `REMARK SMILES` records of the docking output, i.e. with correct
bond orders.

Expected inputs for a full run, all of which the runbooks can generate:

| Path | Content | Produced by |
|---|---|---|
| `<screen>/out/` | AutoDock Vina or Vina-GPU pose files (`xxx_out.pdbqt` or `xxx.pdbqt`) | `runbooks/cpu_pocket_scan.sh`, `runbooks/gpu_run_serial.sh` |
| `affinity.csv` | `id,affinity` table of best score per ligand | see the small helper inside `runbooks/`, or any summary of the pose directory |
| `<receptor>.pdbqt` | receptor in the same coordinate frame as the poses | AutoDockTools / Meeko |
| `ligand_smiles.tsv` | `ligand_id, smiles, ...` tables with trustworthy SMILES | `demos/recover_smiles.py` |
| `library_props.tsv.gz` | descriptor/alert table for a large library | `ligand_qc/prepare_library_props.py` |

**PDBQT carries no bond orders.** Structures recovered from PDBQT alone give wrong descriptors
(a 138 Da drug is reported as 132 Da). `demos/recover_smiles.py` recovers correct SMILES by
matching a neutralised skeleton against an SDF that does contain bond orders
(DrugBank SDF in that campaign); if no such SDF is available, restrict property-based filtering
and use `ligand_qc/ligand_qc.py`, which needs only heavy-atom counts and element-level
connectivity.
