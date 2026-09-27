#!/bin/bash
# Fetch SMILES from ZINC20 2D tranches (the ZINC15 path is dead: zinc15.docking.org now
# 302s to /captcha/, and the ZINC20 *API* also 302s, so only the static export works).
#
# Layout (verified with curl):
#   index : https://files.docking.org/2D/                     -> 121 two-letter dirs (AA..LK)
#   files : https://files.docking.org/2D/<XX>/<XXXX>.txt      -> header-less TSV:
#           smiles  zinc_id  inchikey  mwt  logp  reactive  purchasable  tranche_name  features
#           (the .smi variant only has smiles+zinc_id; .txt is richer and same price)
#           zinc_id has NO 'ZINC' prefix (e.g. 4371221)
#   The 3rd/4th letters of <XXXX> are reactivity and purchasability CLASSES, i.e. properties
#   of the molecule and NOT derivable from the id, so id->filename cannot be computed.
#   Hence: walk each tranche directory, download the .smi files, keep only ids we own.
#
# Usage (LOCAL, needs network):
#   bash fetch_zinc20_smiles.sh --probe            # list size before downloading
#   bash fetch_zinc20_smiles.sh --dirs AA,AB,AC    # restricted, resumable
set -euo pipefail
OUT=${OUT:-/home/workspace/MD/MD/zinc20_smiles}
ZINC_ROOT=${ZINC_ROOT:-/home/workspace/MD/MD/zinc}
JOBS=${JOBS:-4}   # >4 concurrent connections get reset by files.docking.org
DIRS=${DIRS:-}
PROBE=0
while [ $# -gt 0 ]; do
  case "$1" in
    --probe) PROBE=1; shift;;
    --dirs) DIRS="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --jobs) JOBS="$2"; shift 2;;
    *) echo "unknown arg $1" >&2; exit 2;;
  esac
done
mkdir -p "$OUT/raw" "$OUT/by_tranche" "$OUT/logs"
BASE=https://files.docking.org/2D

# library ids WITHOUT the ZINC prefix, for matching ZINC20's bare numeric ids
IDS=$OUT/../zinc_ids_bare.txt
if [ ! -s "$IDS" ]; then
  python3 - "$ZINC_ROOT" "$IDS" <<'PY'
import os, sys
root, dst = sys.argv[1], sys.argv[2]
ids = set()
for base, dirs, files in os.walk(root):
    if os.path.basename(base).startswith("ZINC15_"):
        for f in files:
            if f.endswith(".pdbqt"):
                ids.add(f[4:-6].lstrip("0") or "0")   # ZINC000000307570 -> 307570
open(dst, "w").write("\n".join(sorted(ids)) + "\n")
print("library bare ids: %d" % len(ids))
PY
fi

# discover tranche dirs and their data files.
# NOTE: do not pattern-match the filename against the directory name -- CA/ has 73 entries
# that such a regex misses entirely. Take every href ending in .txt (or .smi) instead, and
# skip directory links. The 3rd/4th letters are reactivity/purchasability classes, so the
# filename cannot be predicted from an id; listing is the only correct method.
if [ -z "$DIRS" ]; then
  DIRS=$(curl -s --max-time 60 "$BASE/" | grep -oE 'href="[A-Z]{2}/"' | tr -d 'href="/' | paste -sd, -)
fi
echo "tranche dirs: $(echo "$DIRS" | tr ',' '\n' | wc -l)"
: > "$OUT/logs/filelist.txt"
for d in ${DIRS//,/ }; do
  body=$(curl -s --retry 5 --retry-all-errors --max-time 90 "$BASE/$d/" || true)
  n=$(printf '%s' "$body" | grep -oE 'href="[A-Za-z0-9_.-]+\.txt"' | sed 's/href="//; s/"//' | sort -u | wc -l)
  [ "$n" = 0 ] && echo "WARN empty/no-listing: $d" >> "$OUT/logs/empty_dirs.txt"
  printf '%s' "$body" | grep -oE 'href="[A-Za-z0-9_.-]+\.txt"' | sed 's/href="//; s/"//' \
    | sort -u | sed "s#^#$d #" >> "$OUT/logs/filelist.txt"
done
echo "smi files discovered: $(wc -l < "$OUT/logs/filelist.txt")"
if [ "$PROBE" = 1 ]; then
  echo "probe sample (first 10, size via curl -I):"
  head -10 "$OUT/logs/filelist.txt" | while read -r d f; do
    sz=$(curl -sI --max-time 20 "$BASE/$d/$f" | awk 'BEGIN{IGNORECASE=1}/[Cc]ontent-[Ll]ength/{print $2}' | tr -d '\r')
    echo "  $d/$f ${sz:-unknown} bytes"
  done
  echo "estimate: $(wc -l < "$OUT/logs/filelist.txt") files x ~0.1-3 MB"
  exit 0
fi

# download + filter to our ids, one output per tranche directory
fetch() {
  d=$1; f=$2; out=$3
  tmp="$out/raw/$f"
  [ -s "$tmp" ] && return 0
  ok=0
  for attempt in 1 2 3 4; do
    if curl -fsSL --retry 5 --retry-all-errors --retry-delay 2 --max-time 300 \
         -o "$tmp.part.$$" "$BASE/$d/$f"; then ok=1; break; fi
    sleep $((attempt * 3))
  done
  if [ "$ok" = 1 ]; then
    mv -f "$tmp.part.$$" "$tmp"
  else
    rm -f "$tmp.part.$$"
    echo "MISS $d/$f" >> "$out/logs/missing.txt"
  fi
}
export -f fetch; export BASE
cat "$OUT/logs/filelist.txt" | xargs -P "$JOBS" -n 2 bash -c 'fetch "$1" "$2" "$0"' "$OUT"

python3 - "$IDS" "$OUT" <<'PY'
import glob, os, sys
want=set(l.strip() for l in open(sys.argv[1]))
out=sys.argv[2]; kept=set()
for f in sorted(glob.glob(os.path.join(out,'raw','*.txt')) + glob.glob(os.path.join(out,'raw','*.smi'))):
    d=os.path.basename(f)[:2]
    dst=os.path.join(out,'by_tranche','%s.tsv'%d)
    new = not os.path.exists(dst)
    with open(f,errors='ignore') as fh, open(dst,'a') as o:
        if new:
            o.write('smiles\tzinc_id\tinchikey\tmwt\tclogp_zinc\n')
        for line in fh:
            p=line.split()
            if len(p)>=2:
                z=p[1].lstrip('ZINC').lstrip('0') or '0'
                if z in want:
                    # .txt rows carry extra columns; keep smiles, id, inchikey, mwt, logp
                    extra = p[2:5] if len(p) >= 5 else ['', '', '']
                    o.write('%s\tZINC%s\t%s\n' % (p[0], z.zfill(12), '\t'.join(extra)))
                    kept.add(z)
print('matched library ids: %d / %d (%.1f%%)' % (len(kept), len(want), 100*len(kept)/max(len(want),1)))
print('NOTE: ZINC20 and the local ZINC15 library are DIFFERENT molecule sets; the measured '
      'overlap on a 13k-row sample was 0.67%, i.e. full coverage of this library is not '
      'expected from ZINC20. Ligands without SMILES fall back to ligand_qc.py (heavy-atom '
      'count + connectivity alerts) or to the ZINC15 API per-id.')
PY
echo "done -> $OUT/by_tranche/ ; next: prepare_library_props.py --smiles-dir $OUT/by_tranche"
