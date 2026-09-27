#!/bin/bash
# Step 1: CPU small-library scan to DEFINE the pocket (fda/drugbank), box deliberately
# generous (24 A) so the pocket shape is measured rather than assumed.
# Next: fit_gpu_box.py turns the top-hit poses into a GPU-sized box.
set -euo pipefail
ROOT=/home/workspace/MD/MD/HDY/ADGRE1
REC=$ROOT/811/receptor.pdbqt
LIG=/home/workspace/MD/MD/drugbank_mk_pdbqt
OUT=$ROOT/cpu_scan/811_b24
LOG=$ROOT/cpu_scan/811_b24_log
mkdir -p "$OUT" "$LOG"
JOBS=8; CPU=8; EXH=16; BOX=24
list=$ROOT/cpu_scan/ligand_list.txt
[ -s "$list" ] || ls "$LIG"/*.pdbqt | sed 's#.*/##; s#\.pdbqt$##' > "$list"
echo "ligands: $(wc -l < "$list")  out=$OUT"
time xargs -a "$list" -P "$JOBS" -I{} \
  vina --receptor "$REC" --ligand "$LIG/{}.pdbqt" \
       --center_x -65.392 --center_y -23.277 --center_z 69.679 \
       --size_x $BOX --size_y $BOX --size_z $BOX \
       --exhaustiveness $EXH --cpu $CPU --num_modes 5 --seed 42 \
       --out "$OUT/{}_out.pdbqt" > "$LOG/{}.log" 2>&1
echo "done: $(ls "$OUT"/*_out.pdbqt 2>/dev/null | wc -l) outputs"
