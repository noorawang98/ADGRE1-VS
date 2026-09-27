#!/bin/bash
# =============================================================================
# Install the deep-screening stack: Chemprop + GPU torch + PyG (+ optional Boltz/DiffDock)
# RUN ON THE HOST (needs network and a writable home). In this container site-packages and
# ~/anaconda3/envs are read-only, so a --prefix environment inside the workspace is used.
#
# Version choices (measured):
#   * driver 535.54.03 supports CUDA 12.x -> use cu124 wheels (torch 2.6.0+cu124 verified)
#   * the existing scvi-pypi env has torch 2.7.1+cu118; mixing cu118 wheels with the system
#     CUDA 12.4 headers breaks PyG extension builds ("CUDA version mismatch"), so do not reuse it
#   * torch-geometric 2.6.x matches torch 2.6; PyG extension wheels come from data.pyg.org
#   * chemprop 2.x requires python >= 3.11 (1.x requires <3.9 and has a different CLI)
# =============================================================================
set -euo pipefail
ENV_NAME=${ENV_NAME:-deepscreen}
WORKSPACE=/home/workspace/MD/MD/HDY/ADGRE1
# Environment location: prefer ~/anaconda3/envs/<name> so conda can activate it by short name;
# fall back to the workspace when that directory is not writable (measured: $HOME is read-only
# in the sandbox, while on the host it is normally writable).
if [ -z "${PREFIX:-}" ]; then
  if mkdir -p "$HOME/anaconda3/envs" 2>/dev/null && [ -w "$HOME/anaconda3/envs" ]; then
    PREFIX="$HOME/anaconda3/envs/$ENV_NAME"
  else
    PREFIX="$WORKSPACE/envs/$ENV_NAME"
  fi
fi
LINK_NAME=${LINK_NAME:-0}   # =1 also creates a symlink under ~/anaconda3/envs (if writable)
CONDA=${CONDA:-$HOME/anaconda3/bin/conda}
TORCH_VER=${TORCH_VER:-2.6.0}
CUDA_TAG=${CUDA_TAG:-cu124}
PYG_TAG="torch-${TORCH_VER}+${CUDA_TAG}"
export PIP_CACHE_DIR=${PIP_CACHE_DIR:-$WORKSPACE/.pipcache}
export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS:-$WORKSPACE/.condapkgs}
mkdir -p "$PIP_CACHE_DIR" "$CONDA_PKGS_DIRS"

SKIP_GPU_CHECK=${SKIP_GPU_CHECK:-0}    # =1 allows continuing without a GPU (CPU torch)
SKIP_SMOKE=${SKIP_SMOKE:-0}            # =1 skips the smoke training run
LOG=${LOG:-$WORKSPACE/setup_deep_env.log}
exec > >(tee -a "$LOG") 2>&1
echo "### setup_deep_env.sh  $(date)  PREFIX=$PREFIX  torch=$TORCH_VER+$CUDA_TAG"

echo "=== [0/6] PREFLIGHT ==="
if [ "$SKIP_GPU_CHECK" = 1 ]; then
  echo "SKIP_GPU_CHECK=1 -> skipping GPU check (CPU torch will be installed)"
else
  nvidia-smi -L || { echo "!! nvidia-smi failed: fix GPU device passthrough first, then re-run.";
                     echo "   For a CPU-only install: SKIP_GPU_CHECK=1 bash $0"; exit 1; }
fi
# A read-only $HOME means conda/uv/pip cannot install anything.
if touch "$HOME/.wtest_$$" 2>/dev/null; then rm -f "$HOME/.wtest_$$"; else
  echo "!! \$HOME ($HOME) is read-only: conda/uv will fail. Run from a writable shell or set"
  echo "   PREFIX= CONDA_PKGS_DIRS= PIP_CACHE_DIR= to writable paths."
fi

# pkgs cache may also be read-only: prefer ~/anaconda3/pkgs, else use the workspace
if [ -z "${CONDA_PKGS_DIRS_OVERRIDE:-}" ]; then
  if mkdir -p "$HOME/anaconda3/pkgs" 2>/dev/null && [ -w "$HOME/anaconda3/pkgs" ]; then
    CONDA_PKGS_DIRS="$HOME/anaconda3/pkgs"
  fi
fi
echo "PREFIX           = $PREFIX"
echo "CONDA_PKGS_DIRS  = $CONDA_PKGS_DIRS"
echo "PIP_CACHE_DIR    = $PIP_CACHE_DIR"
case "$PREFIX" in
  "$HOME"/anaconda3/envs/*) echo "-> standard conda location: conda activate $ENV_NAME" ;;
  *) echo "-> workspace location (fallback). Activate with:"
     echo "     conda activate $PREFIX        # full path"
     echo "     source $PREFIX/bin/activate    # no conda registration needed" ;;
esac

echo "=== [1a] CONDA CHANNEL CHECK (measured: repo.anaconda.com unreachable -> defaults fails) ==="
if curl -s -o /dev/null -w "%{http_code}" --max-time 25 -r 0-0 \
     https://repo.anaconda.com/pkgs/main/linux-64/repodata.json 2>/dev/null | grep -qE "^(200|206)$"; then
  echo "  defaults reachable"
  CH=()
else
  echo "  defaults unreachable -> drop it, use conda-forge only (Tsinghua mirror first)"
  CH=(-c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
      -c https://conda.anaconda.org/conda-forge)
fi

echo "=== [1b] CREATE ENV (--prefix works around read-only envs/ and pkgs/) ==="
if [ -x "$PREFIX/bin/python" ] && "$PREFIX/bin/python" -c "import sys" >/dev/null 2>&1; then
  echo "already present and usable: $PREFIX (skipping creation)"
elif [ -e "$PREFIX" ]; then
  echo "directory exists but is not usable (half-built) -> rebuilding"
  rm -rf "$PREFIX"
  "$CONDA" create -y --prefix "$PREFIX" "${CH[@]}" python=3.11 || true
fi
if [ ! -x "$PREFIX/bin/python" ]; then
  if "$CONDA" create -y --prefix "$PREFIX" "${CH[@]}" python=3.11; then
    echo "conda env created"
  else
    echo "conda create failed -> falling back to uv standalone python (no conda involved)"
    BOOT_PY=$("$CONDA" run -n base python -c "import sys;print(sys.executable)" 2>/dev/null || echo /usr/bin/python3)
    "$BOOT_PY" -m pip install -q -U --user uv 2>/dev/null || \
      "$BOOT_PY" -m pip install -q -U --target "$WORKSPACE/.uvlib" uv
    UV=$(command -v uv || echo "$HOME/.local/bin/uv")
    [ -x "$UV" ] || UV="$BOOT_PY -m uv"
    if [ "$UV" = "$BOOT_PY -m uv" ]; then
      PYTHONPATH="$WORKSPACE/.uvlib" "$BOOT_PY" -m uv venv --seed --python 3.11 "$PREFIX"
    else
      "$UV" venv --seed --python 3.11 "$PREFIX"     # --seed also installs pip/setuptools
    fi
  fi
fi

PY="$PREFIX/bin/python"; PIP="$PREFIX/bin/pip"
"$PY" -V

echo "=== [1c] INSTALLER SETUP (uv preferred; works without pip inside the env) ==="
find_uv () {
  command -v uv 2>/dev/null && return 0
  [ -x "$HOME/.local/bin/uv" ] && { echo "$HOME/.local/bin/uv"; return 0; }
  local BP
  BP=$("$CONDA" run -n base python -c "import sys;print(sys.executable)" 2>/dev/null || echo /usr/bin/python3)
  if [ -n "$BP" ]; then
    "$BP" -m pip install -q --target "$WORKSPACE/.uvlib" uv >/dev/null 2>&1 || true
    if [ -x "$WORKSPACE/.uvlib/bin/uv" ]; then echo "$WORKSPACE/.uvlib/bin/uv"; return 0; fi
    if PYTHONPATH="$WORKSPACE/.uvlib" "$BP" -c "import uv" >/dev/null 2>&1; then
      echo "MODULE:$BP"; return 0
    fi
  fi
  return 1
}
UVSPEC=$(find_uv || true)
if [ -z "$UVSPEC" ]; then
  echo "  !! uv not found and could not be installed. Run: pip install -U uv"; exit 1
fi
if [ "${UVSPEC#MODULE:}" != "$UVSPEC" ]; then
  UV_BOOT=${UVSPEC#MODULE:}; UV_MODE=module
else
  UV_BOOT=$UVSPEC; UV_MODE=bin
fi
echo "  uv source: $UV_BOOT (mode=$UV_MODE)"

# Unified installer: prefer uv pip (no pip needed inside the env), else fall back to $PIP
uvpip () {
  if [ "$UV_MODE" = module ]; then
    PYTHONPATH="$WORKSPACE/.uvlib" "$UV_BOOT" -m uv pip install --python "$PY" "$@"
  else
    "$UV_BOOT" pip install --python "$PY" "$@"
  fi
}

# uv venv does not seed pip, and uv-python builds often omit ensurepip, so bootstrap with uv
if [ ! -x "$PIP" ]; then
  echo "  no pip inside the env; bootstrapping with uv:"
  if uvpip pip setuptools wheel; then echo "  pip installed"; else
    echo "  (pip bootstrap failed; installs will still work through uv pip)"
  fi
fi
[ -x "$PIP" ] && "$PIP" --version || echo "  pip unavailable -> all installs go through uv pip"

echo "=== [2/6] INSTALL CUDA TORCH (first, so chemprop cannot pull a CPU build) ==="
if [ "$SKIP_GPU_CHECK" = 1 ]; then
  uvpip "torch==${TORCH_VER}" torchvision torchaudio
else
  uvpip --index-url "https://download.pytorch.org/whl/${CUDA_TAG}" \
    "torch==${TORCH_VER}" torchvision torchaudio
fi

echo "=== [3/6] VERIFY TORCH SEES THE GPU ==="
"$PY" - <<'PY'
import torch, os
print("torch", torch.__version__, "| cuda build", torch.version.cuda)
ok = torch.cuda.is_available()
print("cuda.is_available:", ok, "| device_count:", torch.cuda.device_count())
if ok:
    print("device:", torch.cuda.get_device_name(0),
          "| cap:", torch.cuda.get_device_capability(0),
          "| mem: %.1f GB" % (torch.cuda.get_device_properties(0).total_memory / 1e9))
    a = torch.randn(1000, 1000, device="cuda"); print("matmul ran:", bool((a @ a).sum() == (a @ a).sum()))
elif os.environ.get("SKIP_GPU_CHECK") == "1":
    print("SKIP_GPU_CHECK=1 -> CPU training accepted (1-2 orders of magnitude slower)")
else:
    print("!! GPU not visible: check /dev/nvidia* and the driver; for CPU only use SKIP_GPU_CHECK=1")
    raise SystemExit(3)
PY

echo "=== [4/6] PyG EXTENSION WHEELS (prebuilt; avoids source builds) ==="
uvpip torch-geometric
uvpip torch-scatter torch-sparse -f "https://data.pyg.org/whl/${PYG_TAG}.html" || \
  echo "(PyG extension wheels unavailable; optional for chemprop v2)"

echo "=== [5/6] CHEMPROP AND SCREENING PACKAGES ==="
uvpip "chemprop>=2.0"
uvpip rdkit numpy scipy pandas scikit-learn pyyaml tensorboard
"$PY" - <<'PY'
import chemprop, torch, rdkit, sys
v = chemprop.__version__
print("chemprop", v, "| torch", torch.__version__, "| rdkit", rdkit.__version__)
if not v.startswith("2."):
    print("!! got chemprop", v, "(1.x CLI differs from the runbook); confirm python >= 3.11")
    sys.exit(4)
PY

echo "=== [6/6] SMOKE TRAINING (2-epoch equivalent; expect under a minute) ==="
if [ "$SKIP_SMOKE" = 1 ]; then echo "SKIP_SMOKE=1 -> skipped"; exit 0; fi
ROOT=$WORKSPACE
if [ -s "$ROOT/deep_demo/input/train.csv" ]; then
  # chemprop 2.x requires epochs > warmup epochs (default warmup=2), so --epochs 2 is rejected
  if [ "$SKIP_GPU_CHECK" = 1 ]; then ACC=cpu; else ACC=gpu; fi
  "$PREFIX/bin/chemprop" train --data-path "$ROOT/deep_demo/input/train.csv" \
    --task-type regression --output-dir "$ROOT/deep_demo/output/smoke" \
    --epochs 5 --warmup-epochs 1 --batch-size 32 --message-hidden-dim 64 --depth 2 \
    --accelerator "$ACC" --devices 1 2>&1 | tail -15
  echo "smoke training finished -> deep_demo/output/smoke"
else
  echo "missing $ROOT/deep_demo/input/train.csv (run deep_demo/scripts/prepare_deep_inputs.py first)"
fi

cat <<'NOTE'

=== OPTIONAL COMPONENTS (uncomment to install) ===
# Boltz-2 co-folding (weights are several GB, tens of seconds per molecule)
#   uvpip boltz
# DiffDock (diffusion pose correction; install inside the cloned repo)
#   git clone https://github.com/gcorso/DiffDock.git && cd DiffDock && uvpip -e .
# Uni-Mol (3D pretrained, weights mol_pre_no_h_220816.pt ~1.5 GB)
#   git clone https://github.com/deepmodeling/Uni-Mol.git && cd Uni-Mol/unimol && uvpip -e .

=== AMBER24 ENGINE (currently missing: only AmberTools modelling tools are present) ===
# ~/softwares/amber24_src exists but is not built; pmemd.cuda and MMPBSA.py are absent.
# With an Amber licence:  export AMBERHOME=$HOME/amber24
#   cd $HOME/softwares/amber24_src/build && ./configure -cuda gnu && make install
#   source $AMBERHOME/amber.sh && which pmemd.cuda MMPBSA.py
# Without a licence: conda install -c conda-forge ambertools (gives MMPBSA.py/cpptraj);
#   use the existing GROMACS 2021.3 (/usr/local/gromacs/bin/gmx) as the MD engine.

=== USAGE ===
export PATH=<PREFIX>/bin:$PATH
chemprop train --data-path deep_demo/input/train.csv --task-type regression \
  --output-dir deep_demo/output/fold1 --epochs 100 --batch-size 32 \
  --message-hidden-dim 300 --depth 3 --split-sizes 0.9 0.1 0.0 --accelerator gpu --devices 1
NOTE
