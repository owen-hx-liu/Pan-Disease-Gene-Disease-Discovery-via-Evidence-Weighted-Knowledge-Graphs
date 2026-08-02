<#
=============================================================================
 run_r2_sweep.ps1 -- evidence against the "R2 collapse is an optimization
 artifact" objection.

 THE OBJECTION
 -------------
 R2 is a degree-preserving, type-preserving rewired null of the R0 training graph
 (Maslov-Sneppen style; scripts/build_degree_null.py). Every KGE model loses a large
 fraction of its MRR there, and the paper reads that as "the topological signal was
 degree". A reviewer can object that this is confounded: the hyperparameters were
 selected on R0, so the R2 drop might be the recipe FAILING TO OPTIMIZE on a rewired
 graph rather than the rewired graph containing no learnable structure.

 The five non-KGE baselines are parameter-free heuristics, so their R2 collapse cannot
 be a tuning artifact at all -- but that argument does not reach the KGE arm, and
 PreferentialAttachment being flat at R2 only proves the rewiring preserved degree
 (it is untrained, so it has no optimization that could fail). This script produces the
 evidence that does reach the KGE arm.

 PHASE A -- train-fit diagnostic (the decisive single measurement)
 ----------------------------------------------------------------
 Score ~4000 sampled TRAINING edges of the target relation with the SAME evaluator that
 produces the test numbers (run_kge.py --train-fit-n). Train fit upper-bounds what the
 optimizer achieved on the graph it was actually given:

   * R2 train fit ~= R0 train fit, while R2 TEST collapses
         -> optimization succeeded; the collapse is a GENERALIZATION result. Objection
            refuted.
   * R2 train fit ALSO collapses
         -> optimization genuinely failed on the rewired graph. The reviewer is right,
            and that must be reported as the finding.

 The full per-epoch training-loss trajectory is recorded for every cell so the R0 and
 R2 curves can be shown side by side.

 Cells: {TransE, DistMult} x {R0, R2}, seed 42, dim 64, 300 epochs -- byte-identical to
 the headline recipe (no recipe flags are passed; DistMult's softplus + regularizer=none
 are already its defaults in run_kge.model_config). The test MRRs are therefore expected
 to REPRODUCE the headline seed-42 runs, which is a free correctness check on the rerun.

 PHASE B -- matched hyperparameter sweep
 ---------------------------------------
 Sweep TransE on R2 AND run the identical grid on R0 as a matched control. A sweep on R2
 alone could not show the gap is config-independent: without the R0 arm there is no
 reference for what the same grid buys on an un-rewired graph.

 Grid: lr in {1e-3, 3e-3, 1e-2} x negatives in {16, 64}, dim 64, 100 epochs, seed 42
 = 6 configs x 2 regimes = 12 cells.

 100 epochs, NOT the headline 300, is deliberate: this sweep is asked to establish the
 config ORDERING (does any configuration close the R0-R2 gap?), not headline numbers.
 Its absolute MRRs are undertrained and must not be quoted as results.

 The dim axis {64,128} from the original plan was CUT. Measured on this GPU (RTX 5070
 Laptop), adding dim 128 at 16 negatives costs a further ~3.2 h, taking A+B to ~12 h
 against a ~6-9 h budget. Consequence, stated plainly: this sweep varies OPTIMIZATION
 hyperparameters (step size, negative count), not model CAPACITY, so it cannot answer
 "would a wider model fit R2?". Phase A does reach that question -- a capacity shortfall
 would show up as depressed R2 TRAIN fit, which is exactly what Phase A measures.

 OUTPUT
 ------
   data/processed/results/kge_r2_sweep/trainfit/   Phase A run JSONs (+ kge_summary.json)
   data/processed/results/kge_r2_sweep/sweep/      Phase B cell JSONs (+ sweep_summary.json)
   data/processed/results/kge_r2_sweep/r2_sweep.log / .status
 Then:  python scripts/r2_sweep_record.py     -> FINDINGS_r2_sweep.md

 Resumable: one JSON per cell; existing cells are skipped unless -Force.

 Usage (detached, survives this shell):
   Start-Process powershell -WindowStyle Hidden `
     -ArgumentList '-NoProfile','-ExecutionPolicy','Bypass','-File','scripts/run_r2_sweep.ps1'
   Get-Content data/processed/results/kge_r2_sweep/r2_sweep.log -Wait
=============================================================================
#>
[CmdletBinding()]
param(
    [string]$Python,
    [ValidateSet('A', 'B', 'AB')]
    [string]$Phase        = 'AB',
    [string[]]$FitModels  = @('TransE', 'DistMult'),
    [string[]]$Regimes    = @('R0', 'R2'),
    [int]$Seed            = 42,
    [int]$TrainFitN       = 4000,
    [int]$FitEpochs       = 300,
    [string]$SweepModel   = 'TransE',
    [int[]]$SweepDims     = @(64),
    [int[]]$SweepEpochs   = @(100),
    [double[]]$SweepLrs   = @(1e-3, 3e-3, 1e-2),
    [int[]]$SweepNegs     = @(16, 64),
    [switch]$Force
)

$RepoRoot = Split-Path -Parent $PSScriptRoot
Set-Location -Path $RepoRoot

if (-not $Python) {
    $venvPy = Join-Path $RepoRoot 'venv\Scripts\python.exe'
    if (Test-Path $venvPy) { $Python = $venvPy } else { $Python = 'python' }
}

$OutDir  = Join-Path $RepoRoot 'data\processed\results\kge_r2_sweep'
$FitDir  = Join-Path $OutDir 'trainfit'
New-Item -ItemType Directory -Force -Path $FitDir | Out-Null
$Log    = Join-Path $OutDir 'r2_sweep.log'
$Status = Join-Path $OutDir 'r2_sweep.status'

$env:PYTORCH_CUDA_ALLOC_CONF = 'expandable_segments:True'

"STARTED $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Out-File -FilePath $Status -Encoding utf8
"  python : $Python"                | Out-File -FilePath $Status -Append -Encoding utf8
"  phase  : $Phase"                 | Out-File -FilePath $Status -Append -Encoding utf8
"  regimes: $($Regimes -join ',')"  | Out-File -FilePath $Status -Append -Encoding utf8
"  seed   : $Seed"                  | Out-File -FilePath $Status -Append -Encoding utf8

# ---------------------------------------------------------------------------
# PHASE A -- train fit vs test fit at the HEADLINE recipe.
# No --regularizer / --*-loss flags: model_config already carries each model's
# published recipe, so this reruns the headline configuration exactly and its test
# MRRs should reproduce the headline seed-42 JSONs.
# ---------------------------------------------------------------------------
if ($Phase -eq 'A' -or $Phase -eq 'AB') {
    "PHASE_A_START $(Get-Date -Format 'HH:mm:ss')" | Out-File -FilePath $Status -Append -Encoding utf8
    $aArgs = @(
        'scripts/run_kge.py',
        '--models') + $FitModels + @(
        '--regimes') + $Regimes + @(
        '--seeds', "$Seed",
        '--dim', '64',
        '--epochs', "$FitEpochs",
        '--train-fit-n', "$TrainFitN",
        '--results-dir', $FitDir
    )
    if ($Force) { $aArgs += '--force' }
    "  cmd: $Python $($aArgs -join ' ')" | Out-File -FilePath $Status -Append -Encoding utf8
    & $Python @aArgs 2>&1 | Out-File -FilePath $Log -Append -Encoding utf8
    "PHASE_A_DONE $(Get-Date -Format 'HH:mm:ss') exit=$LASTEXITCODE" |
        Out-File -FilePath $Status -Append -Encoding utf8
}

# ---------------------------------------------------------------------------
# PHASE B -- the matched grid. --sweep-only skips the headline dim/epoch cell (the
# 300-epoch runs already exist and must not be recomputed); --sweep-regimes carries
# BOTH regimes so R0 is a matched control rather than an absent baseline.
# ---------------------------------------------------------------------------
if ($Phase -eq 'B' -or $Phase -eq 'AB') {
    "PHASE_B_START $(Get-Date -Format 'HH:mm:ss')" | Out-File -FilePath $Status -Append -Encoding utf8
    $bArgs = @(
        'scripts/run_kge.py',
        '--models', $SweepModel,
        '--regimes') + $Regimes + @(
        '--seeds', "$Seed",
        '--sweep-only',
        '--sweep-regimes') + $Regimes + @(
        '--sweep-dims') + ($SweepDims  | ForEach-Object { "$_" }) + @(
        '--sweep-epochs') + ($SweepEpochs | ForEach-Object { "$_" }) + @(
        '--sweep-lrs') + ($SweepLrs  | ForEach-Object { "$_" }) + @(
        '--sweep-negs') + ($SweepNegs | ForEach-Object { "$_" }) + @(
        '--results-dir', $OutDir
    )
    if ($Force) { $bArgs += '--force' }
    "  cmd: $Python $($bArgs -join ' ')" | Out-File -FilePath $Status -Append -Encoding utf8
    & $Python @bArgs 2>&1 | Out-File -FilePath $Log -Append -Encoding utf8
    "PHASE_B_DONE $(Get-Date -Format 'HH:mm:ss') exit=$LASTEXITCODE" |
        Out-File -FilePath $Status -Append -Encoding utf8
}

"DONE $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Out-File -FilePath $Status -Append -Encoding utf8
