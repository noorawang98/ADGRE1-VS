# adgre1-vs -- docking triage and counter-screen toolkit for a 7TM site of ADGRE1

Toolkit distilled from a virtual-screening campaign against the human adhesion GPCR
**ADGRE1 (UniProt Q14246)**, built around three questions that come up in almost every
structure-based screen:

1. **Which docking poses are physically real?** Energy-based ranking alone selects poses whose
   atoms overlap the receptor. A geometry gate removes them; on this target it removed **24.8 %
   of the top 20,000 poses** and **8 of the top 20** compounds.
2. **Where does the ligand actually bind?** A grid box inherited from an earlier run can be tens
   of angstroms away from the site the ligands occupy; two independent checks are provided.
3. **Is the site selective against the closest paralogue?** A hit-rate counter-screen (never a
   cross-site dG subtraction, which is meaningless for two different pockets) quantifies this.

Nothing here is a black box: every decision is a documented threshold or a simple statistical
model, and negative results are reported as such (`RESULTS.md`).

## Contents

- [Layout](#layout)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Usage: workflow scripts](#usage-workflow-scripts)
- [Usage: runbooks](#usage-runbooks)
- [Usage: library preparation](#usage-library-preparation)
- [Usage: environment and deep models](#usage-environment-and-deep-models)
- [Conventions the code relies on](#conventions-the-code-relies-on)
- [Results and caveats](#results-and-caveats)
- [License](#license)

## Layout

```
examples/         runnable example: receptor, six poses, score table (no data needed)
data/             200-ligand example set + input format specification
workflow/         analysis code (CPU; no GPU, no external services)
runbooks/         shell drivers for CPU and GPU screens (host-side, resumable)
ligand_qc/        library-side filtering: properties, structural alerts, ZINC20 retrieval
environment/      environment inventory, mirror preflight, deep-learning env setup
report_figures/   publication figures and pose animations
demos/            self-contained demos of the surrogate and multi-task models
```

## Requirements

| Purpose | Requirement |
|---|---|
| analysis (`workflow/`, `ligand_qc/`) | Python 3.9+, `numpy`, `scipy`; `rdkit` for descriptor and property work |
| docking runs (`runbooks/`) | AutoDock Vina 1.2.3+ (CPU) and/or Vina-GPU 2.1 (GPU); `obabel` for format conversion |
| deep models | Chemprop 2.x (needs Python >= 3.11); installed by `environment/setup_deep_env.sh` |
| figures / animations | `matplotlib`; VMD + `tachyon` for molecular renders, `ffmpeg` for video |

Environments are never assumed. `environment/check_env.sh` prints a capability-by-capability
matrix (GPU visibility, RDKit/scikit-learn availability, data sources) and
`environment/check_mirrors.sh` reports which package indexes and mirrors are reachable.

## Quick start

```bash
# 0. run the committed example (no data of your own needed)
python workflow/select_poses.py \
  --table examples/scores_6.tsv --pose-dirs examples/pose_files \
  --receptor examples/receptor.pdb --min-hac 20 \
  --out examples/selection_example.tsv

# 1. what can this machine do?
bash environment/check_env.sh

# 2. geometry gate on an existing Vina/Vina-GPU output directory
python workflow/select_poses.py \
  --table affinity.tsv --pose-dirs screen/out/ \
  --receptor receptor.pdbqt --top 20000 --out selection.tsv

# 3. paired counter-screen between two receptors
python workflow/reverse_screen_analysis.py \
  --a1 screen_A/ --a1-receptor recA.pdbqt \
  --a2 screen_B/ --a2-receptor recB.pdbqt --tag pair1 --outdir results/

# 4. residue-level ortholog mapping between the two receptors
python workflow/align_pocket_orthologs.py --a1 recA.pdbqt --a2 recB.pdbqt \
  --site-lo 690 --site-hi 710 --outdir results/
```

## Usage: workflow scripts

Common options: `--out` / `--outdir` control every output path, so nothing is written outside the
directory you choose. `--help` lists the full set for each script.

### `workflow/select_poses.py` -- geometry gate and ligand-efficiency ranking

Turns a raw score table plus a directory of poses into a defensible pick list. Poses with any
ligand-receptor heavy-atom pair below 2.0 A (or minimum distance below 2.2 A) are rejected;
ligands below `--min-hac` heavy atoms are rejected because ligand efficiency rewards fragments;
survivors are ranked by ligand efficiency.

| Option | Default | Meaning |
|---|---|---|
| `--table` | required | `id,affinity` or `ligand_id,docking_score` table; comma or tab separated, auto-detected |
| `--pose-dirs` | current dir | directories searched recursively for `*_out.pdbqt`, `*.pdbqt`, `*.pdb` |
| `--receptor` | -- | receptor in the same coordinate frame as the poses |
| `--top` | 20000 | how many of the best-scoring ligands to re-examine |
| `--min-hac` | 20 | heavy-atom floor (0 disables) |
| `--workers` | 16 | process pool size |
| `--out` | `<repo>/selection.tsv` | output table; a `<out>_summary.json` is written beside it |

Output columns: `ligand_id, affinity, le, hac, min_dist, n_clash_atoms, geometry, pass, reason,
contact_residues, spread, n_models, alerts, smiles`. If a pose file carries no `REMARK VINA
RESULT` record, the affinity from `--table` is used, so plain PDB pose sets work too.

### `workflow/pocket_assign.py` -- assign poses to candidate pockets

Tests each pose centroid against candidate boxes **and** requires real receptor contact, so a
large box that contains most of the protein does not label every ligand as "in pocket".

```bash
python workflow/pocket_assign.py --root . --dirs screen/out/ --receptor receptor.pdbqt \
  --outdir triage/
```
Writes one `candidates_<pocket>.tsv` per box plus a multi-pocket consensus table.

### `workflow/fit_gpu_box.py` -- measure the pocket from poses and fit a GPU-sized box

Vina-GPU caps the search space per run, so a blind box cannot be used. This script characterises
the occupied volume from the best poses, reports lining residues and the ligand size
distribution, and writes a fitted box (centre, edge capped by `--box-limit`) plus a profile JSON.

```bash
python workflow/fit_gpu_box.py --scan scan_dir/ --receptor receptor.pdbqt \
  --top 100 --box-limit 24 --outdir fitted/
```

### `workflow/reverse_screen_analysis.py` -- paired counter-screen across two receptors

Computes, per ligand, `dlog = log10|aff_B| - log10|aff_A|` and rank statistics, using only poses
that pass the geometry gate **on both sides**. No dG is subtracted across sites, because the two
pockets are different pockets. Also reports the affinity correlation and the size-corrected
ligand-efficiency correlation between receptors.

```bash
python workflow/reverse_screen_analysis.py \
  --a1 adgre1/out/ --a1-receptor adgre1.pdbqt \
  --a2 adgre2/out/ --a2-receptor adgre2.pdbqt --tag pair1 --outdir results/
```
Accepts several directories per side (`--a1-extra`, `--a2-extra`), a `fdaXXXX -> ZINC` id map
(`--fda-map`) and a `--ddg-cut` on the log ratio (default 0.15, i.e. about 1.4x). Verdicts are
phrased as selective enrichment, never as an absolute energy difference.

### `workflow/align_pocket_orthologs.py` -- residue-level ortholog mapping

Aligns the receptor sequences (BLOSUM62 global alignment), locates the 7TM region from
hydrophobicity, and returns a per-residue correspondence for the site of interest plus the
centroid of the orthologous residues, which is what a cross-receptor box centre needs.

```bash
python workflow/align_pocket_orthologs.py --a1 recA.pdbqt --a2 recB.pdbqt \
  --site-lo 690 --site-hi 710 --pad 40 --outdir results/
```

### `workflow/surrogate_model.py` -- descriptor surrogate for docking scores

Trains gradient-boosted models on Morgan-fingerprint SVD features plus physico-chemical
descriptors, with a property-only control model so that "specificity" (pocket-specific signal
beyond molecular size) can be separated from trivial size effects.

```bash
python workflow/surrogate_model.py --pocket 811 --out triage/ --top 500
```

### `workflow/prepare_deep_inputs.py` + `workflow/evaluate_deep_scores.py`

Scaffold-grouped train/val/test splits for Chemprop (prefers measured activity labels in
`activity.tsv` over docking labels when present), and metric computation with a measured
GO/NO-GO bar.

```bash
python workflow/prepare_deep_inputs.py --dataset dataset.tsv --folds 5 --outdir input/
python workflow/evaluate_deep_scores.py --pred test_preds.csv --test input/test.csv \
  --out metrics.json
```

## Usage: runbooks

All runbooks are resumable (existing outputs are skipped) and take a pilot measurement before any
full run. Run them on the machine that owns the GPU; the analysis scripts stay on CPU.

| Runbook | Purpose |
|---|---|
| `runbooks/cpu_pocket_scan.sh` | CPU scan of a small library to define a pocket before committing GPU time |
| `runbooks/gpu_reverse_prepare.sh` | stage ligand directories and write both Vina-GPU configs; `prepare` needs no GPU |
| `runbooks/gpu_sweep.sh` | measure 4 (thread, search_depth) combinations on 200 ligands and report s/ligand |
| `runbooks/gpu_run_serial.sh` | run batches **one instance at a time** (two instances on one card measured 5x slower) |
| `runbooks/gpu_background.sh` | one-shot background pipeline: stop old jobs -> sweep -> serial runs -> analysis |
| `runbooks/stop_gpu.sh` | stop all Vina-GPU jobs and verify output growth has stopped |
| `runbooks/run_deep_train.sh` | Chemprop 5-fold training, ensembling, metrics and library ranking |
| `runbooks/run_eval_deep.sh` | wait for folds, produce ensemble predictions and the GO/NO-GO verdict |

```bash
# GPU reverse screen, pilot first
bash runbooks/gpu_reverse_prepare.sh prepare
bash runbooks/gpu_sweep.sh                      # or the 200-ligand pilot in the docs
GPU_BATCH=fda bash runbooks/gpu_reverse_prepare.sh run-adgre2
GPU_BATCH=fda bash runbooks/gpu_reverse_prepare.sh run-adgre1

# deep model: training, ensemble, metrics, ranking
bash runbooks/run_deep_train.sh
bash runbooks/run_eval_deep.sh
```

## Usage: library preparation

| Script | Purpose |
|---|---|
| `ligand_qc/ligand_qc.py` | bond-order-free QC: heavy-atom counts and element-level toxicophore flags, usable when only PDBQT exists |
| `ligand_qc/prepare_library_props.py` | full descriptor table (MW, cLogP, TPSA, HBD/HBA, RotB, QED), PAINS/Brenk and toxicophore alerts, approved-drug property windows |
| `ligand_qc/fetch_zinc20_smiles.sh` | download ZINC20 2D tranches and keep only ids present in a local library |
| `ligand_qc/coverage_report.py` | report how much of a local library a downloaded SMILES set covers |

```bash
python ligand_qc/ligand_qc.py --affinity affinity.out --out qc.tsv
python ligand_qc/prepare_library_props.py --smiles-dir smiles_by_tranche/ --out library_props.tsv.gz
bash ligand_qc/fetch_zinc20_smiles.sh --probe --dirs AA,AB,AC
python ligand_qc/coverage_report.py --dir zinc20_out/ --lib library_ids.txt
```

## Usage: environment and deep models

```bash
bash environment/check_env.sh          # capability inventory (safe, read-only)
bash environment/check_mirrors.sh      # which package indexes/mirrors are reachable
bash environment/setup_deep_env.sh     # Chemprop + GPU torch + PyG; falls back to uv if conda fails
SKIP_GPU_CHECK=1 bash environment/setup_deep_env.sh    # CPU-only install
```
`setup_deep_env.sh` handles the failure modes measured in this project: unreachable conda
`defaults`, read-only `$HOME`, and `uv venv` shipping without pip (installs go through
`uv pip --python <env>`).

## Conventions the code relies on

* **Score orientation**: lower is better (Vina convention) everywhere.
* **Geometry gate**: reject when any ligand-receptor heavy-atom pair is below **2.0 A**, or the
  minimum distance is below **2.2 A**. Clash artefacts often carry *better* scores than valid
  poses, so the gate is applied before ranking.
* **Ranking key**: ligand efficiency `LE = -dG / heavy-atom count`, with a heavy-atom floor,
  because raw scores grow with molecular size and LE alone rewards fragments.
* **Cross-receptor comparisons**: `dlog = log10|aff_B| - log10|aff_A|` and rank/percentile
  statistics, never a difference of dG values.
* **Engine consistency**: both sides of a pair must be docked with the same engine and the same
  search settings, otherwise engine differences leak into the selectivity estimate.
* **Contacts**: a residue is "in contact" when a ligand heavy atom is within **4.0 A**.

## Results and caveats

`RESULTS.md` lists every measured number, including two negative results that should not be
hidden: a machine-learning surrogate did **not** beat ligand-efficiency ranking on a
scaffold-grouped test set, and the two paralogue sites proved functionally equivalent, which rules
out selectivity at that particular site. `REFERENCE.md` documents the thresholds and the residue
tables so they can be reviewed rather than trusted. `examples/` runs without any external data.

## License

MIT -- see `LICENSE`. No experimental data are included; see `data/README.md`.
