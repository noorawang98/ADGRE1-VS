#!/bin/bash
# Connectivity and writability preflight: decides whether to use the pip route or conda.
# Measured context: repo.anaconda.com and conda.anaconda.org are unreachable from this host
# (HTTP 000) while conda's default channels are bioconda/conda-forge/defaults -> conda fails.
# This script only READS state and writes a log; it changes no configuration.
probe () { printf "  %-4s %s\n" "$(curl -s -o /dev/null -w '%{http_code}' --max-time 25 -r 0-0 "$1" 2>/dev/null)" "$1"; }
writable () { touch "$1/.wtest" 2>/dev/null && { echo "  writable OK  $1"; rm -f "$1/.wtest"; } || echo "  READ-ONLY    $1"; }

OUT=${OUT:-/home/workspace/MD/MD/HDY/ADGRE1/check_mirrors.log}
exec > >(tee -a "$OUT") 2>&1
echo "### check_mirrors.sh  $(date)  (read-only preflight; no configuration is modified)"

echo "=== 1. NETWORK CHANNELS ==="
probe https://pypi.org/simple/
probe https://pypi.tuna.tsinghua.edu.cn/simple/
probe https://download.pytorch.org/whl/cu124/
probe https://data.pyg.org/whl/torch-2.6.0+cu124.html
probe https://repo.anaconda.com/pkgs/main/linux-64/repodata.json
probe https://conda.anaconda.org/conda-forge/linux-64/repodata.json
probe https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/linux-64/repodata.json
probe https://zinc.docking.org/

echo "=== 2. WRITABLE DIRECTORIES ==="
writable "$HOME/anaconda3/envs"
writable "$HOME/anaconda3/pkgs"
writable /home/workspace/MD/MD/HDY/ADGRE1

echo "=== 3. GPU ==="
nvidia-smi -L 2>/dev/null || echo "  GPU not visible (missing /dev/nvidia*)"
ls /dev/nvidia* 2>/dev/null || echo "  /dev/nvidia* absent"
clinfo 2>/dev/null | grep -m1 "Number of platforms" || echo "  clinfo unavailable (or OpenCL platforms = 0)"

echo "=== 4. RECOMMENDATION (no configuration changed) ==="
code () { curl -s -o /dev/null -w "%{http_code}" --max-time 25 -r 0-0 "$1" 2>/dev/null; }
echo "  --- reachability (30 s timeout; 000 usually means a transient mirror hiccup) ---"
printf "    %-22s %s\n" "PyPI"                "$(code https://pypi.org/simple/)"
printf "    %-22s %s\n" "PyPI (Tsinghua)"     "$(code https://pypi.tuna.tsinghua.edu.cn/simple/)"
printf "    %-22s %s\n" "PyTorch cu124"       "$(code https://download.pytorch.org/whl/cu124/)"
printf "    %-22s %s\n" "PyG wheels"          "$(code https://data.pyg.org/whl/torch-2.6.0+cu124.html)"
printf "    %-22s %s\n" "anaconda defaults"   "$(code https://repo.anaconda.com/pkgs/main/linux-64/repodata.json)"
printf "    %-22s %s\n" "conda-forge official" "$(code https://conda.anaconda.org/conda-forge/linux-64/repodata.json)"
printf "    %-22s %s\n" "conda-forge Tsinghua" "$(code https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/linux-64/repodata.json)"
echo "  --- recommended route ---"
echo "    1) pip first: python -m venv / uv venv + pip (PyPI direct; torch from download.pytorch.org)"
echo "       This reaches chemprop 2.x (requires python >= 3.11)."
echo "    2) conda fallback: add ONLY the Tsinghua conda-forge mirror, drop defaults/bioconda:"
echo "       conda create -y --prefix <workspace>/envs/deepscreen \\"
echo "         -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge python=3.11"
echo "    3) last resort: setup_deep_env.sh falls back to 'uv venv --python 3.11' (no conda at all)"
echo "    NOTE: $HOME/anaconda3/envs and $HOME/anaconda3/pkgs are read-only here, so --prefix"
echo "          into the workspace plus a custom CONDA_PKGS_DIRS are required."
