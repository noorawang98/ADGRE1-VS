#!/bin/bash
# =============================================================================
# Parameter sweep: measure 4 (thread, search_depth) combinations with 200 ligands,
# report seconds per ligand and a projection for the full workload.
# RUN ON HOST:  bash gpu_sweep.sh
# Finished combinations are skipped, so the script can be re-run to fill gaps.
#
# NOTE on timing: seconds-per-ligand is derived from the model runtime reported by
# Vina-GPU ("AutoDockVina-GPU3 total runtime = ... s"), NOT from file mtimes -- Vina-GPU
# flushes outputs at the end of a batch, so mtime spans give nonsense (measured 0.07 s).
# =============================================================================
set -uo pipefail
ROOT=/home/workspace/MD/MD/HDY/ADGRE1
OUT=$ROOT/reverse_screen
VINA_GPU=${VINA_GPU:-$HOME/softwares/Vina-GPU-2.1-main/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1}
REC=/home/workspace/MD/MD/HDY/ADGRE2/receptor.pdbqt
SWEEP=$OUT/sweep
mkdir -p "$SWEEP"
LOG=$ROOT/gpu_sweep.log
exec > >(tee -a "$LOG") 2>&1
echo "### gpu_sweep.sh $(date)"
nvidia-smi -L || { echo "GPU not visible, aborting"; exit 1; }
bash "$ROOT/gpu_reverse_prepare.sh" prepare >/dev/null 2>&1 || true

LIG=$SWEEP/ligs
if [ "$(ls $LIG/*.pdbqt 2>/dev/null | wc -l)" -lt 200 ]; then
  mkdir -p "$LIG"; ls $OUT/lig_fda/*.pdbqt | head -200 | xargs -I{} cp -n {} "$LIG/" 2>/dev/null
fi
N=$(ls $LIG/*.pdbqt | wc -l); echo "benchmark ligands: $N; receptor: ADGRE2 corrected box (-50.7,-5.2,54.9), 24 A"

write_cfg () {  # $1=file $2=thread $3=depth
cat > "$1" <<ARG
center_x = -50.7
center_y = -5.2
center_z = 54.9
size_x = 24
size_y = 24
size_z = 24
thread = $2
search_depth = $3
ARG
}

for cfg in "8000 20" "8000 10" "3000 20" "3000 10"; do
  set -- $cfg; TH=$1; DP=$2
  DST=$SWEEP/t${TH}_d${DP}
  mkdir -p "$DST"
  if [ "$(ls $DST/*.pdbqt 2>/dev/null | wc -l)" -ge $((N-5)) ]; then
    echo "=== thread=$TH depth=$DP : results already present, skipping run ==="
  else
    write_cfg "$SWEEP/cfg.arg" "$TH" "$DP"
    echo "=== thread=$TH depth=$DP : running ==="
    "$VINA_GPU" --config "$SWEEP/cfg.arg" --receptor "$REC" \
      --ligand_directory "$LIG" --output_directory "$DST" \
      --opencl_binary_path "$(dirname "$VINA_GPU")" 2>&1 | tee "$SWEEP/log_t${TH}_d${DP}.txt" | tail -3
    echo "  runtime line: $(grep -ao 'total runtime = [0-9.]*' "$SWEEP/log_t${TH}_d${DP}.txt" | tail -1)"
  fi
done

echo; echo "=== SUMMARY (seconds per ligand -> projection) ==="
python3 - <<'PY'
import glob, os, re
base='/home/workspace/MD/MD/HDY/ADGRE1/reverse_screen/sweep'
rows=[]
for lg in sorted(glob.glob(base+'/log_t*_d*.txt')):
    tag=os.path.basename(lg)[4:-4]
    m=re.search(r'total runtime = ([0-9.]+)', open(lg, errors='ignore').read())
    n=len(glob.glob(lg.replace('log_','').replace('.txt','')+'/*.pdbqt'))
    if m and n>50:
        per=float(m.group(1))/n
        rows.append((tag,n,per))
rows.sort(key=lambda r:r[2])
print("  %-12s %6s %11s %13s %12s"%("config","ligands","s/ligand","full 23632","fda only 4204"))
for tag,n,per in rows:
    print("  %-12s %6d %11.2f %10.1f h %10.1f h"%(tag,n,per,23632*per/3600,4204*per/3600))
if rows:
    print("\n  FASTEST: %s (%.2f s/ligand)"%(rows[0][0],rows[0][2]))
    if len(rows)>1:
        print("  speedup vs slowest: %.1fx"%(rows[-1][2]/rows[0][2]))
print("\n  NOTE: both receptors must be re-run with the SAME config, otherwise parameter")
print("        differences leak into the cross-receptor comparison.")
PY
echo "### done $(date)   log: gpu_sweep.log"
