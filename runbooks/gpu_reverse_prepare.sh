#!/bin/bash
# =============================================================================
# GPU reverse-screen preparation: stage ligand directories, write both Vina-GPU
# configs, and align output paths with the CPU side so the analysis scripts and the
# monitor can read GPU results without modification.
#
#   bash gpu_reverse_prepare.sh                 # prepare only (no GPU needed)
#   GPU_BATCH=fda      bash gpu_reverse_prepare.sh run-adgre2
#   GPU_BATCH=fda      bash gpu_reverse_prepare.sh run-adgre1
#   GPU_BATCH=drugbank bash gpu_reverse_prepare.sh run-adgre2
#   GPU_BATCH=drugbank bash gpu_reverse_prepare.sh run-adgre1
#
# Output paths (identical to the CPU side):
#   reverse_screen/adgre2_cpu/<batch>/   reverse_screen/adgre1_cpu/<batch>/
#
# Speed parameters (from measurement): the pilot ran 200 ligands in 306 s = 1.53 s/ligand,
# while two concurrent instances dropped to 7.2 ligands/min = 8.3 s/ligand (5x slower) --
# caused by both instances competing on one GPU plus an oversized thread count.
#   2) search_depth 20 -> 10   (halves sampling; sufficient for reverse screening; MUST match on both sides)
#   3) thread       8000 -> 3000 (8000 is oversaturated on a 10 GB RTX 3080)
# =============================================================================
set -uo pipefail
ROOT=/home/workspace/MD/MD/HDY/ADGRE1
OUT=$ROOT/reverse_screen
VINA_GPU=${VINA_GPU:-$HOME/softwares/Vina-GPU-2.1-main/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1}
[ -x "$VINA_GPU" ] || VINA_GPU=$HOME/softwares/Vina-GPU-2.1-main/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1
THREAD=${THREAD:-3000}; DEPTH=${DEPTH:-10}; BOX=${BOX:-24}
A2_REC=/home/workspace/MD/MD/HDY/ADGRE2/receptor.pdbqt
A1_REC=$ROOT/811/receptor.pdbqt
# Corrected box centers (see SOLUTION.md)
CX2=-50.7; CY2=-5.2; CZ2=54.9      # ADGRE2 orthologs 572-593 centroid
CX1=-54.3; CY1=-10.5; CZ1=54.2     # ADGRE1 contacted residues 690-710 centroid
BATCH=${GPU_BATCH:-fda}

prepare () {
  mkdir -p "$OUT/logs" "$OUT/lig_fda" "$OUT/lig_drugbank"
  if [ "$(ls $OUT/lig_fda 2>/dev/null | wc -l)" -lt 2000 ]; then
    cp -n /home/wr/scripts/vinalc_test/data/fda/*.pdbqt "$OUT/lig_fda/" 2>/dev/null || true
  fi
  if [ "$(ls $OUT/lig_drugbank 2>/dev/null | wc -l)" -lt 9000 ]; then
    for f in /home/workspace/MD/MD/drugbank_mk_pdbqt/*.pdbqt; do
      ln -sf "$f" "$OUT/lig_drugbank/$(basename "$f")"
    done
  fi
  write_cfg () {  # $1=file $2..4=center
    cat > "$1" <<ARG
center_x = $2
center_y = $3
center_z = $4
size_x = $BOX
size_y = $BOX
size_z = $BOX
thread = $THREAD
search_depth = $DEPTH
ARG
  }
  write_cfg "$OUT/cfg_adgre2_site.arg" $CX2 $CY2 $CZ2
  write_cfg "$OUT/cfg_adgre1_site.arg" $CX1 $CY1 $CZ1
  echo "ligands: fda $(ls $OUT/lig_fda/*.pdbqt 2>/dev/null|wc -l) | drugbank $(ls $OUT/lig_drugbank/*.pdbqt 2>/dev/null|wc -l)"
  echo "configs: cfg_adgre2_site.arg (box $BOX A, center $CX2,$CY2,$CZ2)"
  echo "         cfg_adgre1_site.arg (box $BOX A, center $CX1,$CY1,$CZ1)"
  echo "thread=$THREAD  search_depth=$DEPTH"
}

run_one () {   # $1=adgre2|adgre1  $2=batch
  local side=$1 batch=$2 cfg rec dest lig
  if [ "$side" = adgre2 ]; then cfg=$OUT/cfg_adgre2_site.arg; rec=$A2_REC; dest=$OUT/adgre2_cpu/$batch; else
    cfg=$OUT/cfg_adgre1_site.arg; rec=$A1_REC; dest=$OUT/adgre1_cpu/$batch; fi
  lig=$OUT/lig_$batch
  mkdir -p "$dest"
  [ "$(ls $lig/*.pdbqt 2>/dev/null | wc -l)" -gt 0 ] || { echo "!! missing ligand dir $lig"; exit 1; }
  echo "=== GPU $side / $batch -> $dest ($(ls $lig/*.pdbqt|wc -l) ligands) ==="
  time "$VINA_GPU" --config "$cfg" --receptor "$rec" \
      --ligand_directory "$lig" --output_directory "$dest" \
      --opencl_binary_path "$(dirname "$VINA_GPU")" 2>&1 | tail -20
  echo "outputs: $(ls $dest/*_out.pdbqt 2>/dev/null | wc -l)"
}

case "${1:-prepare}" in
  prepare) prepare; echo; echo "next (host): bash reverse_screen_runbook.sh pilot" ;;
  run-adgre2) run_one adgre2 "$BATCH" ;;
  run-adgre1) run_one adgre1 "$BATCH" ;;
  *) echo "usage: bash gpu_reverse_prepare.sh [prepare|run-adgre2|run-adgre1]  (GPU_BATCH=fda|drugbank)"; exit 2 ;;
esac
