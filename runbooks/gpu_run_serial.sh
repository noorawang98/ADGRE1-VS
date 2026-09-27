#!/bin/bash
# =============================================================================
# SERIAL runner: one Vina-GPU instance at a time, so the two sides do not compete for the
# same GPU (which measured 5x slower). Sides/batches run in sequence and the analysis for
# each batch runs right after its two sides finish.
#
# RUN ON HOST:  bash gpu_run_serial.sh [fda|drugbank|all]
# Override parameters:  THREAD=3000 DEPTH=10 bash gpu_run_serial.sh all
# Quick self-test (first N ligands only):  LIMIT=50 SKIP_SWEEP=1 bash gpu_run_serial.sh fda
# =============================================================================
set -uo pipefail
ROOT=/home/workspace/MD/MD/HDY/ADGRE1
MODE=${1:-all}
LOG=$ROOT/gpu_run_serial.log
LIMIT=${LIMIT:-0}
exec > >(tee -a "$LOG") 2>&1
echo "### gpu_run_serial.sh $(date)  MODE=$MODE  THREAD=${THREAD:-3000} DEPTH=${DEPTH:-10}  LIMIT=$LIMIT"
bash "$ROOT/gpu_reverse_prepare.sh" prepare || exit 1
grep -H "thread\|search_depth" "$ROOT"/reverse_screen/cfg_*.arg

run () {  # $1=batch $2=side
  echo; echo "=============== $2 / $1  $(date +%H:%M) ==============="
  local t0=$(date +%s)
  if [ "$LIMIT" -gt 0 ]; then
    local tmp=$ROOT/reverse_screen/lig_$1_smoke
    rm -rf "$tmp"; mkdir -p "$tmp"
    ls "$ROOT/reverse_screen/lig_$1"/*.pdbqt | head -"$LIMIT" | xargs -I{} cp -n {} "$tmp/" 2>/dev/null
    local cfg=$ROOT/reverse_screen/cfg_$2_site.arg
    local rec=/home/workspace/MD/MD/HDY/ADGRE2/receptor.pdbqt
    [ "$2" = adgre1 ] && rec=$ROOT/811/receptor.pdbqt
    if [ "$2" = adgre2 ]; then local dst=$ROOT/reverse_screen/adgre2_cpu/$1; else local dst=$ROOT/reverse_screen/adgre1_cpu/$1; fi
    mkdir -p "$dst"
    echo "  [self-test] $LIMIT ligands -> $dst"
    "$HOME/softwares/Vina-GPU-2.1-main/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1" \
      --config "$cfg" --receptor "$rec" --ligand_directory "$tmp" --output_directory "$dst" \
      --opencl_binary_path "$HOME/softwares/Vina-GPU-2.1-main/AutoDock-Vina-GPU-2.1" 2>&1 | tail -5
  else
    GPU_BATCH="$1" bash "$ROOT/gpu_reverse_prepare.sh" "run-$2" | tail -6
  fi
  local t1=$(date +%s) d n
  if [ "$2" = "adgre2" ]; then d=$ROOT/reverse_screen/adgre2_cpu/$1; else d=$ROOT/reverse_screen/adgre1_cpu/$1; fi
  n=$(ls $d/*.pdbqt 2>/dev/null | wc -l)
  echo "  $2/$1: wall $((t1-t0)) s, cumulative outputs $n"
}

analyse () {  # $1=batch
  echo; echo "--------------- analysis: $1 ----------------"
  ~/anaconda3/envs/AmberTools24/bin/python "$ROOT/reverse_screen_analysis.py" \
    --a1 "$ROOT/reverse_screen/adgre1_cpu/$1" --a2 "$ROOT/reverse_screen/adgre2_cpu/$1" \
    --tag "gpu_$1" 2>&1 | tail -14
}

case "$MODE" in
  fda)      run fda adgre2; run fda adgre1; [ "$LIMIT" -gt 0 ] || analyse fda ;;
  drugbank) run drugbank adgre2; run drugbank adgre1; [ "$LIMIT" -gt 0 ] || analyse drugbank ;;
  all)      run fda adgre2; run fda adgre1; [ "$LIMIT" -gt 0 ] || analyse fda
            run drugbank adgre2; run drugbank adgre1; [ "$LIMIT" -gt 0 ] || analyse drugbank ;;
  *) echo "usage: bash gpu_run_serial.sh [fda|drugbank|all]"; exit 2 ;;
esac
echo; echo "### done $(date)  summaries: reverse_screen/summary_gpu_*.json"
