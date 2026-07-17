<#
=============================================================================
 run_all.ps1 -- reproduce the leakage-audit benchmark end to end (Unit 12).

 Runs, in order:
   0. hash-verify the input artifacts against ARTIFACT_HASHES.txt (ABORT on mismatch)
   1. graph_stats.py            -> data/processed/graph_stats.json
   2. graph_stats_hetionet.py   -> data/processed/graph_stats_hetionet.json
   3. build_deleaked_splits.py  -> data/processed/splits/*.csv + split_manifest.json
   4. build_degree_null.py      -> data/processed/results/null/degree_null.json (+ R2 split)
   5. run_baselines.py R0-R3    -> data/processed/results/baselines/baselines_R*.json
   6. run_kge.py (GPU)          -> data/processed/results/kge/*.json + kge_summary.json
   7. hetionet_audit.py R0-R2   -> data/processed/results/hetionet/baselines_R*.json
   8. degree_stratified.py      -> data/processed/results/degree_stratified/degree_stratified.json
   9. case_study.py             -> data/processed/results/case_study/case_study.json
  10. make_figures.py           -> figures/fig1-4.{png,pdf}
  11. make_tables.py            -> tables/table1-6.{md,tex}
  12. hash-verify reproduced splits (informational)

 Acceptance: a clean clone + the data channel reproduces Table 2 and Fig 2.
 Runtimes below are approximate on a modern workstation; the KGE step needs a
 CUDA GPU and dominates the wall-clock. Paths are repo-relative.

 Usage:
   ./run_all.ps1                  # full pipeline (KGE included -- needs a GPU)
   ./run_all.ps1 -SkipKge         # everything except KGE (reuses committed kge_summary.json)
   ./run_all.ps1 -VerifyOnly      # only run the hash checks, then exit
   ./run_all.ps1 -SkipHashCheck   # skip integrity check (NOT recommended)
   ./run_all.ps1 -Python venv\Scripts\python.exe   # choose the interpreter
=============================================================================
#>
[CmdletBinding()]
param(
    [switch]$SkipKge,
    [switch]$SkipHashCheck,
    [switch]$VerifyOnly,
    [string]$Python
)

$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot

# Default interpreter: the repo venv if present, else python on PATH.
if (-not $Python) {
    $venvPy = Join-Path $PSScriptRoot "venv\Scripts\python.exe"
    if (Test-Path $venvPy) { $Python = $venvPy } else { $Python = "python" }
}

$script:Step = 0

function Banner([string]$Text) {
    Write-Host ""
    Write-Host "=============================================================================="
    Write-Host "  $Text"
    Write-Host "=============================================================================="
}

function Invoke-Step {
    param(
        [string]$Title,
        [string]$Est,
        [string]$Out,
        [Parameter(Mandatory)][string[]]$PyArgs
    )
    $script:Step++
    Banner "STEP $($script:Step): $Title"
    Write-Host "  expected runtime : $Est"
    Write-Host "  writes           : $Out"
    Write-Host "  command          : $Python $($PyArgs -join ' ')"
    Write-Host ""
    & $Python @PyArgs
    if ($LASTEXITCODE -ne 0) {
        throw "step $($script:Step) failed (exit $LASTEXITCODE): $Title"
    }
    Write-Host ""
    Write-Host "  [ok] step $($script:Step) done -> $Out"
}

try {
    Write-Host "run_all.ps1  |  repo: $PSScriptRoot"
    Write-Host "interpreter : $Python  ($(& $Python --version 2>&1))"

    # ------------------------------------------------------------ 0. hash-verify
    Banner "PRE-FLIGHT: verify input artifacts against ARTIFACT_HASHES.txt"
    if ($SkipHashCheck) {
        Write-Host "  -SkipHashCheck set: SKIPPING integrity check (NOT recommended)."
    } else {
        & $Python scripts/verify_hashes.py --require `
            data/processed/edges_clean_integrated.csv `
            data/external/hetionet/hetionet-v1.0-edges.sif `
            data/external/hetionet/hetionet-v1.0-nodes.tsv
        if ($LASTEXITCODE -ne 0) { throw "input hash verification failed -- ABORT." }
    }

    if ($VerifyOnly) {
        Banner "VERIFY-ONLY: checking every ARTIFACT_HASHES.txt entry present on disk"
        & $Python scripts/verify_hashes.py --check-existing
        exit 0
    }

    # ------------------------------------------------------------------ pipeline
    Invoke-Step -Title "Graph statistics (Monarch)" -Est "~1-2 min" `
        -Out "data/processed/graph_stats.json" `
        -PyArgs @("scripts/graph_stats.py")

    Invoke-Step -Title "Graph statistics (Hetionet robustness graph)" -Est "<1 min" `
        -Out "data/processed/graph_stats_hetionet.json" `
        -PyArgs @("scripts/graph_stats_hetionet.py")

    Invoke-Step -Title "Build de-leaked splits R0/R1/R3 (seed 42)" -Est "~5-15 min" `
        -Out "data/processed/splits/*.csv + split_manifest.json" `
        -PyArgs @("scripts/build_deleaked_splits.py",
                  "--edges", "data/processed/edges_clean_integrated.csv",
                  "--out", "data/processed/splits", "--seed", "42")

    Invoke-Step -Title "Build R2 degree-preserving null (10 replicates, seed 42)" -Est "~5-15 min" `
        -Out "data/processed/results/null/degree_null.json + splits/train_R2_degree_null_seed42.csv" `
        -PyArgs @("scripts/build_degree_null.py",
                  "--splits-dir", "data/processed/splits", "--replicates", "10", "--seed", "42")

    Invoke-Step -Title "Topological baselines across regimes R0,R1,R2,R3 (seeds 42,1,7)" -Est "~3-8 min" `
        -Out "data/processed/results/baselines/baselines_R{0,1,2,3}.json" `
        -PyArgs @("scripts/run_baselines.py",
                  "--regimes", "R0,R1,R2,R3", "--seeds", "42,1,7",
                  "--splits-dir", "data/processed/splits")

    # -------------------------------------------------------------------- 6. KGE
    if ($SkipKge) {
        $script:Step++
        Banner "STEP $($script:Step): KGE (TransE + RotatE)  --  SKIPPED (-SkipKge)"
        Write-Host "  Reusing the committed data/processed/results/kge/kge_summary.json."
        Write-Host "  To (re)compute on a CUDA GPU (~40 min/cell, 24 cells -> many hours):"
        Write-Host "    $Python scripts/run_kge.py --models TransE RotatE --regimes R0 R1 R2 R3 ``"
        Write-Host "        --seeds 42 1 7 --dim 64 --epochs 300"
    } else {
        Invoke-Step -Title "KGE full benchmark: TransE + RotatE x R0,R1,R2,R3 x 3 seeds" `
            -Est "GPU REQUIRED; ~40 min/cell, 24 cells -> many hours" `
            -Out "data/processed/results/kge/*.json + kge_summary.json" `
            -PyArgs @("scripts/run_kge.py",
                      "--models", "TransE", "RotatE",
                      "--regimes", "R0", "R1", "R2", "R3",
                      "--seeds", "42", "1", "7", "--dim", "64", "--epochs", "300")
    }

    Invoke-Step -Title "Hetionet cross-graph robustness (DaG; R0,R1,R2)" -Est "~2-5 min" `
        -Out "data/processed/results/hetionet/baselines_R{0,1,2}.json" `
        -PyArgs @("scripts/hetionet_audit.py", "--regimes", "R0,R1,R2", "--seeds", "42,1,7")

    # ---------------------------------------------- 8. degree-stratified re-analysis
    Invoke-Step -Title "Degree-stratified re-analysis (gene vs disease degree; no retraining)" `
        -Est "<1 min" `
        -Out "data/processed/results/degree_stratified/degree_stratified.json" `
        -PyArgs @("scripts/degree_stratified.py")

    # ---------------------------------------------------- 9. worked case study
    Invoke-Step -Title "Worked case study (degree false-positives vs missed rare edges)" `
        -Est "<1 min" `
        -Out "data/processed/results/case_study/case_study.json" `
        -PyArgs @("scripts/case_study.py")

    # ----------------------------------------------------------- 10-11. deliverables
    Invoke-Step -Title "Manuscript figures (Fig 1-4, png + pdf)" -Est "<1 min" `
        -Out "figures/fig1-4.{png,pdf}" `
        -PyArgs @("scripts/make_figures.py")

    Invoke-Step -Title "Manuscript tables (Tables 1-6 + S1, md + tex)" -Est "<1 min" `
        -Out "tables/table1-6.{md,tex}" `
        -PyArgs @("scripts/make_tables.py")

    # ------------------------------------------------- 10. verify reproduction
    Banner "POST-RUN: verify reproduced artifacts against ARTIFACT_HASHES.txt"
    & $Python scripts/verify_hashes.py --check-existing
    if ($LASTEXITCODE -eq 0) {
        Write-Host "  [ok] reproduced splits are byte-identical to the frozen manifest."
    } else {
        Write-Host "  [warn] some reproduced artifacts did not match the manifest -- see above."
    }

    Banner "DONE. Outputs:"
    Write-Host "  graph stats   : data/processed/graph_stats.json, graph_stats_hetionet.json"
    Write-Host "  splits        : data/processed/splits/"
    Write-Host "  degree null   : data/processed/results/null/degree_null.json"
    Write-Host "  baselines     : data/processed/results/baselines/"
    Write-Host "  KGE           : data/processed/results/kge/"
    Write-Host "  Hetionet      : data/processed/results/hetionet/"
    Write-Host "  figures       : figures/fig1-3.{png,pdf}"
    Write-Host "  tables        : tables/table1-4.{md,tex}"
}
catch {
    Write-Host ""
    Write-Host "run_all.ps1 ABORTED: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
}
