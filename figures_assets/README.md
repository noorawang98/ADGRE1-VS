# Figure assets

The campaign produced slide-ready figures and pose animations:

| File | Content |
|---|---|
| `fig1_funnel.png` | screening funnel: 488,810 -> top 20 -> 12 passing the geometry gate -> 5 assayed -> 1 active |
| `fig2_distribution.png` | library score distribution plus a zoom on the best 0.5 %, with the assayed compounds marked |
| `fig3_gate_le.png` | affinity versus ligand efficiency, clash artefacts flagged, candidates labelled |
| `fig4_contacts.png` | contact map of the selected compound and its selection metrics |
| `fig5_selectivity.png` | paired counter-screen, ADGRE1 site versus paralogue ortholog |
| `fig6_candidate_structures.png` | 2D structures of the five assayed candidates with score, ligand efficiency, heavy atoms, minimum distance, clash pairs and geometry verdict |
| `anim_binding_site_spin.mp4` | pose animation of the binding site (VMD + tachyon, 36 frames) |
| `five_candidates_pose.png` | the five assayed compounds overlaid in the same pocket |

`fig6_candidate_structures.png` **is** committed (compressed, 1.2 MB) together with the input
table `report_figures/five_candidates.json`, because it is the slide that explains the compound
choice. The other figures and the animation are large binaries; regenerate them with
`report_figures/make_figures.py` (statistics) and the VMD script in
`report_figures/README.md` (animation).
