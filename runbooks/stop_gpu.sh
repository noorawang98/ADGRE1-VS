#!/bin/bash
# =============================================================================
# Stop all Vina-GPU jobs. RUN ON THE HOST (the container cannot signal host processes).
#   bash stop_gpu.sh
# Verifies that output counts stop growing afterwards.
# =============================================================================
echo "=== before: matching processes ==="
pgrep -af "Vina-GPU" || echo "  (none)"
echo "=== stopping ==="
pkill -f "AutoDock-Vina-GPU-2-1"
pkill -f "AutoDock-Vina-GPU"
pkill -f "Vina-GPU"
sleep 2
pkill -9 -f "AutoDock-Vina-GPU" 2>/dev/null
sleep 2
echo "=== after ==="
pgrep -af "Vina-GPU" && echo "  !! processes still alive" || echo "  all stopped OK"
echo
echo "=== verify output growth has stopped (3 samples, 10 s apart) ==="
for i in 1 2 3; do
  echo "  $(date +%H:%M:%S)  total $(ls /home/workspace/MD/MD/HDY/ADGRE1/reverse_screen/*/fda/*.pdbqt 2>/dev/null | wc -l)"
  sleep 10
done
