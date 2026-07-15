#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- reproduce the leakage-audit benchmark end to end (Unit 12).
#
# Runs, in order:
#   0. hash-verify the input artifacts against ARTIFACT_HASHES.txt  (ABORT on mismatch)
#   1. graph_stats.py            -> data/processed/graph_stats.json
#   2. graph_stats_hetionet.py   -> data/processed/graph_stats_hetionet.json
#   3. build_deleaked_splits.py  -> data/processed/splits/*.csv + split_manifest.json
#   4. build_degree_null.py      -> data/processed/results/null/degree_null.json (+ R2 split)
#   5. run_baselines.py R0-R3    -> data/processed/results/baselines/baselines_R*.json
#   6. run_kge.py (GPU)          -> data/processed/results/kge/*.json + kge_summary.json
#   7. hetionet_audit.py R0-R2   -> data/processed/results/hetionet/baselines_R*.json
#   8. make_figures.py           -> figures/fig1-3.{png,pdf}
#   9. make_tables.py            -> tables/table1-4.{md,tex}
#  10. hash-verify reproduced splits (informational)
#
# Acceptance: a clean clone + the data channel reproduces Table 2 and Fig 2.
# Runtimes below are approximate on a modern workstation; the KGE step needs a
# CUDA GPU and dominates the wall-clock. Paths are repo-relative; run from anywhere.
#
# Usage:
#   ./run_all.sh                 # full pipeline (KGE included -- needs a GPU)
#   ./run_all.sh --skip-kge      # everything except KGE (reuses committed kge_summary.json)
#   ./run_all.sh --verify-only   # only run the hash checks, then exit
#   ./run_all.sh --skip-hash-check   # skip integrity check (NOT recommended)
#   PYTHON=venv/bin/python ./run_all.sh   # choose the interpreter
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

PYTHON="${PYTHON:-python}"
SKIP_KGE=0
SKIP_HASH=0
VERIFY_ONLY=0

usage() { sed -n '2,29p' "$0"; }

while [ $# -gt 0 ]; do
  case "$1" in
    --skip-kge)        SKIP_KGE=1 ;;
    --skip-hash-check) SKIP_HASH=1 ;;
    --verify-only)     VERIFY_ONLY=1 ;;
    -h|--help)         usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage; exit 2 ;;
  esac
  shift
done

STEP=0
banner() {
  echo
  echo "=============================================================================="
  echo "  $1"
  echo "=============================================================================="
}

# run_py "<title>" "<expected runtime>" "<output location>" <command...>
run_py() {
  local title="$1" est="$2" out="$3"; shift 3
  STEP=$((STEP + 1))
  banner "STEP ${STEP}: ${title}"
  echo "  expected runtime : ${est}"
  echo "  writes           : ${out}"
  echo "  command          : $*"
  echo
  "$@"
  echo
  echo "  [ok] step ${STEP} done -> ${out}"
}

echo "run_all.sh  |  repo: ${ROOT}"
echo "interpreter : ${PYTHON}  ($("$PYTHON" --version 2>&1))"

# ---------------------------------------------------------------- 0. hash-verify
banner "PRE-FLIGHT: verify input artifacts against ARTIFACT_HASHES.txt"
if [ "$SKIP_HASH" -eq 1 ]; then
  echo "  --skip-hash-check set: SKIPPING integrity check (NOT recommended)."
else
  "$PYTHON" scripts/verify_hashes.py --require \
      data/processed/edges_clean_integrated.csv \
      data/external/hetionet/hetionet-v1.0-edges.sif \
      data/external/hetionet/hetionet-v1.0-nodes.tsv
fi

if [ "$VERIFY_ONLY" -eq 1 ]; then
  banner "VERIFY-ONLY: checking every ARTIFACT_HASHES.txt entry present on disk"
  "$PYTHON" scripts/verify_hashes.py --check-existing || true
  exit 0
fi

# --------------------------------------------------------------------- pipeline
run_py "Graph statistics (Monarch)" "~1-2 min" \
  "data/processed/graph_stats.json" \
  "$PYTHON" scripts/graph_stats.py

run_py "Graph statistics (Hetionet robustness graph)" "<1 min" \
  "data/processed/graph_stats_hetionet.json" \
  "$PYTHON" scripts/graph_stats_hetionet.py

run_py "Build de-leaked splits R0/R1/R3 (seed 42)" "~5-15 min" \
  "data/processed/splits/*.csv + split_manifest.json" \
  "$PYTHON" scripts/build_deleaked_splits.py \
      --edges data/processed/edges_clean_integrated.csv \
      --out data/processed/splits --seed 42

run_py "Build R2 degree-preserving null (10 replicates, seed 42)" "~5-15 min" \
  "data/processed/results/null/degree_null.json + splits/train_R2_degree_null_seed42.csv" \
  "$PYTHON" scripts/build_degree_null.py \
      --splits-dir data/processed/splits --replicates 10 --seed 42

run_py "Topological baselines across regimes R0,R1,R2,R3 (seeds 42,1,7)" "~3-8 min" \
  "data/processed/results/baselines/baselines_R{0,1,2,3}.json" \
  "$PYTHON" scripts/run_baselines.py \
      --regimes R0,R1,R2,R3 --seeds 42,1,7 --splits-dir data/processed/splits

# --------------------------------------------------------------------- 6. KGE
if [ "$SKIP_KGE" -eq 1 ]; then
  STEP=$((STEP + 1))
  banner "STEP ${STEP}: KGE (TransE + RotatE)  --  SKIPPED (--skip-kge)"
  echo "  Reusing the committed data/processed/results/kge/kge_summary.json."
  echo "  To (re)compute on a CUDA GPU (~40 min/cell, 24 cells -> many hours):"
  echo "    $PYTHON scripts/run_kge.py --models TransE RotatE --regimes R0 R1 R2 R3 \\"
  echo "        --seeds 42 1 7 --dim 64 --epochs 300"
else
  run_py "KGE full benchmark: TransE + RotatE x R0,R1,R2,R3 x 3 seeds" \
    "GPU REQUIRED; ~40 min/cell, 24 cells -> many hours" \
    "data/processed/results/kge/*.json + kge_summary.json" \
    "$PYTHON" scripts/run_kge.py --models TransE RotatE --regimes R0 R1 R2 R3 \
        --seeds 42 1 7 --dim 64 --epochs 300
fi

run_py "Hetionet cross-graph robustness (DaG; R0,R1,R2)" "~2-5 min" \
  "data/processed/results/hetionet/baselines_R{0,1,2}.json" \
  "$PYTHON" scripts/hetionet_audit.py --regimes R0,R1,R2 --seeds 42,1,7

# ---------------------------------------------------------------- 8-9. deliverables
run_py "Manuscript figures (Fig 1-3, png + pdf)" "<1 min" \
  "figures/fig1-3.{png,pdf}" \
  "$PYTHON" scripts/make_figures.py

run_py "Manuscript tables (Tables 1-4, md + tex)" "<1 min" \
  "tables/table1-4.{md,tex}" \
  "$PYTHON" scripts/make_tables.py

# ------------------------------------------------------- 10. verify reproduction
banner "POST-RUN: verify reproduced artifacts against ARTIFACT_HASHES.txt"
if "$PYTHON" scripts/verify_hashes.py --check-existing; then
  echo "  [ok] reproduced splits are byte-identical to the frozen manifest."
else
  echo "  [warn] some reproduced artifacts did not match the manifest -- see above."
fi

banner "DONE. Outputs:"
cat <<'EOF'
  graph stats   : data/processed/graph_stats.json, graph_stats_hetionet.json
  splits        : data/processed/splits/
  degree null   : data/processed/results/null/degree_null.json
  baselines     : data/processed/results/baselines/
  KGE           : data/processed/results/kge/
  Hetionet      : data/processed/results/hetionet/
  figures       : figures/fig1-3.{png,pdf}
  tables        : tables/table1-4.{md,tex}
EOF
