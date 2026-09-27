#!/bin/bash
# =============================================================================
# One-shot background GPU pipeline:
#   stop old jobs -> parameter sweep (pick fastest) -> serial fda -> serial drugbank
#   -> automatic analysis
#
# RUN ON HOST:
#   nohup bash gpu_background.sh > /dev/null 2>&1 &
#   tail -f status.txt gpu_background.log          # progress
#
# Options:
#   SKIP_SWEEP=1        skip the sweep, use current cfg thread/depth
#   ONLY=fda            run only the fda batch
#   THREAD=3000 DEPTH=10   override parameters
# =============================================================================
set -uo pipefail
ROOT=/home/workspace/MD/MD/HDY/ADGRE1
cd "$ROOT"
LOG=$ROOT/gpu_background.log
ST=$ROOT/status.txt
exec > >(tee -a "$LOG") 2>&1
stage () { echo "$(date +%F' '%H:%M:%S)  $1" | tee -a "$ST"; echo "--- $1 ---"; }

stage "START GPU background pipeline  THREAD=${THREAD:-3000} DEPTH=${DEPTH:-10}  ONLY=${ONLY:-all}"

stage "1/5 stop any running Vina-GPU jobs"
bash "$ROOT/stop_gpu.sh" 2>&1 | tail -8 || true

if [ "${SKIP_SWEEP:-0}" != 1 ]; then
  stage "2/5 parameter sweep (200 ligands x 4 combos, finished combos are skipped)"
  bash "$ROOT/gpu_sweep.sh" 2>&1 | tail -25
  BEST=$(grep -oE "fastest: t[0-9]+_d[0-9]+|FASTEST: t[0-9]+_d[0-9]+" "$ROOT/gpu_sweep.log" 2>/dev/null | tail -1 | grep -oE "t[0-9]+_d[0-9]+")
  if [ -n "${BEST:-}" ]; then
    export THREAD=$(echo "$BEST" | sed 's/^t//; s/_d.*//')
    export DEPTH=$(echo "$BEST" | sed 's/.*_d//')
    stage "sweep picked thread=$THREAD depth=$DEPTH (cfg files are rewritten below)"
  else
    stage "sweep inconclusive -> keep THREAD=${THREAD:-3000} DEPTH=${DEPTH:-10}"
  fi
else
  stage "2/5 sweep skipped (SKIP_SWEEP=1)"
fi

export THREAD=${THREAD:-3000} DEPTH=${DEPTH:-10}
bash "$ROOT/gpu_reverse_prepare.sh" prepare >/dev/null 2>&1
grep -H "thread\|search_depth" "$ROOT"/reverse_screen/cfg_*.arg | tee -a "$ST"

case "${ONLY:-all}" in
  fda)      stage "3/5 serial run: fda both sides";      bash "$ROOT/gpu_run_serial.sh" fda ;;
  drugbank) stage "3/5 serial run: drugbank both sides"; bash "$ROOT/gpu_run_serial.sh" drugbank ;;
  *)        stage "3/5 serial run: fda both sides";      bash "$ROOT/gpu_run_serial.sh" fda
            stage "4/5 serial run: drugbank both sides"; bash "$ROOT/gpu_run_serial.sh" drugbank ;;
esac

stage "5/5 done; summary files:"
ls -l "$ROOT"/reverse_screen/summary_gpu_*.json 2>/dev/null | tee -a "$ST" || echo "  (no summary, see log)" | tee -a "$ST"
echo "$(date +%F' '%H:%M:%S)  ALL DONE" | tee -a "$ST"
