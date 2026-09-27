#!/bin/bash
# Environment inventory for the ADGRE1 7TM screening pipeline.
# Reports, per capability, whether it is present in THIS container and where.
echo "=== SYSTEM ==="
echo "cpu: $(nproc) cores | mem: $(free -g 2>/dev/null | awk 'NR==2{print $2" GB"}') | disk: $(df -h /home/workspace | awk 'NR==2{print $4" free"}')"
echo "gpu: $(nvidia-smi -L 2>/dev/null || echo 'UNAVAILABLE (no /dev/nvidia*)') | opencl platforms: $(clinfo 2>/dev/null | awk '/Number of platforms/{print $4}')"
echo
echo "=== CLI TOOLS ==="
for t in vina obabel gmx pmemd sander MMPBSA.py python3 conda git curl wget; do
  printf "  %-12s %s\n" "$t" "$(command -v $t 2>/dev/null || echo MISSING)"
done
echo
echo "=== CONDA ENVS AND KEY PACKAGES ==="
printf "  %-16s %-10s %-10s %-10s %-10s %-10s\n" ENV rdkit numpy scipy sklearn torch
for e in adfr AmberTools24 scvi-pypi cellpose-gpu meeko hovernet; do
  p=~/anaconda3/envs/$e/bin/python
  [ -x "$p" ] || { printf "  %-16s %s\n" "$e" "(absent)"; continue; }
  line=$(for m in rdkit numpy scipy sklearn torch; do
      $p -c "import $m,sys;print(getattr($m,'__version__','?'))" 2>/dev/null || echo "-"
    done | paste -sd' ' -)
  printf "  %-16s %s\n" "$e" "$line"
done
echo
echo "=== DEEP-LEARNING / SCREENING PACKAGES (expected missing) ==="
for m in chemprop torch_geometric unimol deepchem boltz diffdock; do
  found=""
  for e in adfr AmberTools24 scvi-pypi cellpose-gpu; do
    p=~/anaconda3/envs/$e/bin/python
    [ -x "$p" ] && $p -c "import $m" 2>/dev/null && found="$found $e"
  done
  printf "  %-18s %s\n" "$m" "${found:-MISSING (in every env)}"
done
echo
echo "=== GPU TOOLCHAIN ==="
for t in nvcc "$HOME/softwares/Vina-GPU-2.1-main/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1" "/usr/local/gromacs/bin/gmx" "$HOME/softwares/vmd-1.9.3"; do
  [ -e "$t" ] && printf "  present: %s\n" "$t" || printf "  MISSING: %s\n" "$t"
done
echo
echo "=== DATA SOURCES ==="
for d in /home/workspace/MD/MD/drugbank_mk_pdbqt /home/workspace/MD/MD/drugbank.sdf \
         /home/workspace/MD/MD/zinc /home/workspace/MD/MD/zinc15_smiles /home/workspace/MD/MD/HDY/ADGRE2/out; do
  [ -e "$d" ] && printf "  present: %-52s %s files\n" "$d" "$(find "$d" -maxdepth 1 | wc -l)" || printf "  MISSING: %s\n" "$d"
done
