# Report figures and pose animation

## Statistical figures (CPU only)

`make_figures.py` reads the numbers produced by the workflow and writes 300 dpi PNGs:

```bash
python report_figures/make_figures.py --help
python report_figures/make_figures.py
```

Inputs it expects (paths configurable at the top of the script): a library score table
(`id,affinity`), the top-20 geometry table (from `workflow/select_poses.py`) and the paired
counter-screen output (`workflow/reverse_screen_analysis.py`).

## Structure panel for the assayed compounds (RDKit)

`make_structure_figure.py` draws the 2D structures of the assayed candidates and annotates each
with its docking score, ligand efficiency, heavy-atom count, minimum ligand-receptor distance,
clash-pair count and geometry verdict; the last cell carries the selection rationale.

```bash
python report_figures/make_structure_figure.py \
  --data report_figures/five_candidates.json --out candidate_structures.png
```

Input is a JSON list with `ligand_id, smiles, affinity, hac, le, min_dist, clash, geometry`;
`five_candidates.json` in this directory is the set used in the campaign. SMILES must carry
correct bond orders (they were taken from the `REMARK SMILES` records of the docking output).

## Pose animation (VMD + tachyon)

1. Prepare PDB files in a common coordinate frame: the receptor, the ligand of interest, and (for
   the comparison frame) each candidate ligand. `demos/recover_smiles.py` is unrelated to this
   step; use any PDBQT -> PDB conversion that keeps the columns aligned, e.g.

   ```bash
   obabel receptor.pdbqt -O receptor.pdb
   obabel ligand_out.pdbqt -O ligand.pdb -r   # -r keeps only the first model
   ```

2. Render the rotating scene. VMD's internal renderer locks output to the display size, so a scene
   file is written and traced externally, which is what gives full resolution:

   ```tcl
   # spin.tcl
   mol new receptor.pdb type pdb waitfor all
   mol new ligand.pdb  type pdb waitfor all
   mol modstyle 0 0 NewCartoon
   mol modstyle 1 0 Licorice
   color Display Background white
   display resetview
   scale by 1.5
   molinfo top set center_matrix [molinfo top get center_matrix]
   for {set i 0} {$i < 36} {incr i} {
     rotate y by 10
     render Tachyon [format "scene_%04d.dat" $i]
   }
   ```

   ```bash
   vmd -dispdev text -e spin.tcl
   for f in scene_*.dat; do
     tachyon_LINUXAMD64 -aasamples 32 -res 1280 960 -format TARGA -o "${f%.dat}.tga" "$f"
   done
   ffmpeg -framerate 12 -i scene_%04d.tga -c:v libx264 -crf 18 -pix_fmt yuv420p spin.mp4
   ```

3. `tga2gif.py` converts TGA frame directories into GIFs if a video container is not wanted.

Notes: one tachyon frame at 1280x960 with 32x anti-aliasing takes about 1.1 s; VMD's
`display resize` and `render ... width/height` are ignored by the internal renderer in the build
used here, which is the reason for the two-step route above.
