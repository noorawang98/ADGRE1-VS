# Measured results

All numbers below come from the campaign that produced this toolkit and are reproducible from
the workspace artefacts listed with each item. They are reported together with the negative
results, because those are the most transferable part of the exercise.

## 1. Library and screening scale

| Item | Value |
|---|---|
| Library screened (Vina-GPU) | **488,810** compounds, one pocket, 30 A box |
| Second library screened (CPU Vina, corrected box) | 9,714 DrugBank compounds |
| Receptors docked in the counter-screen | 2 (ADGRE1 site and its paralogue ortholog), 23,604 poses total |
| Wall time, GPU counter-screen | ~5.1 h on one RTX 3080 (thread 3000, search_depth 10) |

## 2. Clash artefacts dominate the raw ranking

Energy-only ranking is unreliable at this site:

| Measure | Value |
|---|---|
| Top-20,000 poses rejected by the geometry gate | **4,851 (24.8 %)** |
| Of the top-20 compounds by score | **8 (35 %)** have ligand-receptor heavy-atom distances of 0.61-1.85 A |
| Effect on selection | the compound taken into the assay was the best-scoring candidate **whose pose passes the gate** (min distance 2.59 A, zero atom pairs < 2.0 A) |

Reproduce with `workflow/select_poses.py` on any Vina/Vina-GPU output directory.

## 2b. The five assayed candidates, one by one

| ligand | score | LE | heavy atoms | min distance | clash pairs < 2.0 A | geometry |
|---|---|---|---|---|---|---|
| **ZINC000002446170** (chosen) | **-12.10** | **0.336** | 36 | **2.59 A** | **0** | **pass** |
| ZINC000009015892 | -12.10 | 0.327 | 37 | 2.85 A | 0 | pass |
| ZINC000019782239 | -12.10 | 0.346 | 35 | 1.56 A | 5 | reject |
| ZINC000015956900 | -12.30 | 0.293 | 42 | 0.71 A | 5 | reject |
| ZINC000009424484 | -12.20 | 0.290 | 42 | 0.77 A | 7 | reject |

Three of the five hand-picked compounds have poses whose atoms overlap the receptor; the chosen
compound is the best-scoring candidate that also passes the geometry gate. Reproduce with
`workflow/select_poses.py` and plot with `report_figures/make_structure_figure.py`.

## 3. The grid box was 23 A from the occupied site

| Measure | Value |
|---|---|
| Box centre of the earlier pocket definition | (-65.4, -23.3, 69.7) |
| Centroid of the residues actually contacted by docked ligands | (-54.3, -10.5, 54.2) |
| Distance | **23.0 A** |
| Consequence | that box returned a median score around -0.05 kcal/mol and produced clash artefacts; after recentring, the same library returned a normal distribution with a best score of -8.5 kcal/mol for a known active compound |

Two independent checks are provided: contact analysis (`workflow/pocket_assign.py`) and cavity
detection (voxel flood fill; no enclosed cavity larger than 40 A^3 was found on this receptor,
i.e. the site is a surface groove rather than a buried pocket).

## 4. Machine learning did not improve enrichment (negative result)

Same data, scaffold-grouped 5-fold cross-validation, 4,414 ligands with descriptors:

| Task / model | R2 | MAE | Spearman |
|---|---|---|---|
| affinity, ridge | 0.612 | 0.347 | 0.761 |
| affinity, kNN-10 (Tanimoto) | **0.690** | 0.314 | **0.786** |
| paralogue-pair affinity, kNN-10 | 0.749 | 0.415 | 0.839 |
| site selectivity (dlog), kNN-10 | 0.095 | 0.030 | 0.333 |

Enrichment against a size-matched decoy set: ML **EF1% 3.18** versus plain ligand efficiency
**5.00**. On the scaffold-grouped test set the ligand-efficiency baseline itself is weak
(precision@10% = 0.045, Spearman -0.259 against the docking labels), so the honest statement is:
**the surrogate reproduces docking scores but does not rank better than the simple physical
descriptor, and it must not be presented as an enrichment improvement.**

Reproduce with `demos/ai_screen_surrogate.py` and `demos/ai_screen_multitask.py`.

## 5. The paralogue sites are functionally equivalent (negative result with consequences)

Counter-screen of the ADGRE1 site against its orthologous site in the closest paralogue, with
identical engine and box settings on both sides:

| Metric | FDA batch | DrugBank batch |
|---|---|---|
| pairs passing the geometry gate on both sides | 978 | **5,752** |
| affinity correlation | 0.854 | 0.891 |
| **ligand-efficiency correlation (size-corrected)** | **0.965** | **0.967** |
| median score difference | -0.070 | -0.060 |
| ligands preferring the ADGRE1 site | 0 / 978 | **5 / 5,752 (0.09 %)** |

Removing the size term *increases* the correlation, so this is not a trivial size effect. The
site pairs are functionally equivalent and no small-molecule selectivity window exists at that
position; a programme targeting it selectively should change site or modality. Reproduce with
`workflow/reverse_screen_analysis.py`; the ortholog mapping that defines the second box comes
from `workflow/align_pocket_orthologs.py` (identity 0.545 over the 7TM region).

## 6. Deep model: pipeline prepared, baseline fixed before use

* Chemprop 2.3.1 / torch 2.6.0+cu124 environment builds and trains (`environment/setup_deep_env.sh`,
  5-epoch smoke run completes).
* Training labels are docking scores, so any model trained on them is a surrogate of that engine.
* The GO/NO-GO bar is measured, not assumed: **Spearman > 0.6 and precision@10% > 0.045**
  (the ligand-efficiency baseline on the same scaffold-grouped test set).
* `runbooks/run_deep_train.sh` performs training, ensembling, metric computation and library
  ranking in one command.

## 7. Throughput measurements that shaped the runbooks

| Measurement | Value |
|---|---|
| CPU cost per ligand (exhaustiveness 16) | ~20 CPU-seconds, independent of `--cpu` |
| CPU throughput, 64 cores saturated | 6.7 ligands/min |
| GPU, 200-ligand pilot | 1.53 s/ligand |
| GPU, two instances on one card | 8.3 s/ligand (**5x slower**) |
| GPU, serial with thread 3000 / depth 10 | 23,604 poses in 5.1 h |

Consequences implemented in the runbooks: run one Vina-GPU instance at a time, prefer a
site-centred 24 A box over a blind box, and measure with a 200-ligand pilot before committing.
Grid reuse (`--maps` in this Vina build) gave only a 9 % speed-up and changed the grid spacing,
so it is not used.
