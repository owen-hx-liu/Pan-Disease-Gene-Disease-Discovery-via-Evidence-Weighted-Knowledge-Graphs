#!/usr/bin/env bash
# run_fullrank_kgearm.sh -- complete the filtered full-ranking arm of the KGE benchmark.
#
# The existing data/processed/results/kge_fullrank/ holds TransE x {R0,R1,R2} x {42,1,7},
# produced by `run_kge.py --full-rank` (dim 64, epochs 300, no subgraph, batch 16384).
# This script fills the two gaps through the IDENTICAL code path:
#   * RotatE x {R0,R1,R2,R3} x {42,1,7}   (12 cells)
#   * TransE x {R3}         x {42,1,7}   ( 3 cells)
#
# One invocation PER REGIME (not per cell): run_kge.load_regimes() holds every requested
# regime's 5.8M-edge training array in RAM at once, and main() builds the ~3-minute
# category-pool cache per regime. Per-regime invocations keep peak RAM at one graph and
# pay the pool build once, while run_kge's own "skip (exists)" resume logic checkpoints
# each cell to its own JSON, so an interrupted run restarts at the first missing cell.
#
# --validate-fullrank 150 is attached to the FIRST regime of each model, which
# brute-forces the vectorized full-pool rank against a per-candidate scorer loop
# (run_kge._validate_fullrank) and aborts on any disagreement.
#
# Budget on the RTX 5070, from the train_seconds of the runs already on disk:
#   RotatE ~4460 s/cell (~3905 s on R3), TransE ~2350 s/cell (~2050 s on R3),
#   plus ~14 min of full-ranking evaluation per cell -> ~20 GPU-hours total.
#
# Usage:  bash scripts/run_fullrank_kgearm.sh
# Tail :  tail -f data/processed/results/kge_fullrank/run_fullrank_kgearm.log
set -u

cd "$(dirname "$0")/.." || exit 1
PY=./venv/Scripts/python.exe
LOG=data/processed/results/kge_fullrank/run_fullrank_kgearm.log
mkdir -p "$(dirname "$LOG")"

run() {
    echo "" >> "$LOG"
    echo "=== $(date '+%Y-%m-%d %H:%M:%S') :: run_kge.py --full-rank $* ===" >> "$LOG"
    "$PY" -u scripts/run_kge.py --full-rank "$@" >> "$LOG" 2>&1
    echo "=== exit $? :: $(date '+%Y-%m-%d %H:%M:%S') ===" >> "$LOG"
}

echo "=== KGE full-rank arm started $(date) ===" >> "$LOG"

# ---- RotatE: all four regimes x three seeds -------------------------------- #
run --models RotatE --regimes R0 --seeds 42 1 7 --validate-fullrank 150
run --models RotatE --regimes R1 --seeds 42 1 7
run --models RotatE --regimes R2 --seeds 42 1 7
run --models RotatE --regimes R3 --seeds 42 1 7

# ---- TransE: the one regime missing from the frozen full-rank set ---------- #
run --models TransE --regimes R3 --seeds 42 1 7 --validate-fullrank 150

# ---- final merge over every frozen per-cell JSON --------------------------- #
run --merge-only

echo "=== KGE full-rank arm finished $(date) ===" >> "$LOG"
