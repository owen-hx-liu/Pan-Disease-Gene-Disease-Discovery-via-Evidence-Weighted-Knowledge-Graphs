<#
=============================================================================
 run_bilinear_grid.ps1 -- headline grid for the bilinear KGE models
 (DistMult, Yang et al. ICLR 2015; ComplEx, Trouillon et al. ICML 2016).

 WHY THIS EXISTS
 ---------------
 The earlier DistMult grid (scripts/run_distmult.ps1) and the ComplEx rescue runs
 collapsed to the trivial all-scores-equal optimum: training loss pinned at exactly
 the loss evaluated at zero score gap (softplus -> ln 2 = 0.6931; NSSA margin 9 ->
 4.5; margin-ranking 1.0 -> 1.0) from the first epoch onward.

 A single-variable sweep over the whole recipe (loss in {softplus, NSSA, margin-
 ranking} x batch in {16384, 4096} x negatives in {16, 64} x lr in {1e-3, 1e-2})
 showed the collapse is invariant to ALL of them. The cause is PyKEEN's per-model
 default regularizer, which TransE and RotatE do not carry:

   TransE    entity_constrainer=normalize          regularizer: NONE
   RotatE    relation_constrainer=complex_normalize regularizer: NONE
   DistMult  entity_constrainer=normalize          relation LpRegularizer w=0.1
   ComplEx   (no constrainer)                      entity+relation Lp w=0.01

 DistMult clamps entity embeddings to unit norm, so |score| <= |r| and the relation
 embeddings are the ONLY source of score scale; an unopposed w=0.1 L2 penalty on them
 drives |r| -> 0, after which the entity gradient (proportional to r (*) t) vanishes
 and the model is frozen. ComplEx shows the same signature more slowly (|r| 10.84 ->
 0.10). Setting regularizer=none therefore does not privilege the bilinear models --
 it brings them to PARITY with TransE and RotatE, which are unregularized. Every other
 hyperparameter (sLCWA, NSSA loss, dim 64, 300 epochs, batch 16384, 16 negatives,
 PyKEEN 1.11.1) is left byte-identical to the translational recipe.

 Grid: {DistMult, ComplEx} x {R0,R1,R2,R3} x seeds {42,1,7} = 24 cells.
 --rank-seed is fixed at 42 (run_kge.py default) so candidate sets stay byte-identical
 to every other method in the benchmark. Resumable: one JSON per cell, existing cells
 are skipped unless -Force.

 Logs   -> data/processed/results/kge/bilinear_grid.log
 Status -> data/processed/results/kge/bilinear_grid.status

 Usage (detached, survives this shell):
   Start-Process powershell -WindowStyle Hidden `
     -ArgumentList '-NoProfile','-ExecutionPolicy','Bypass','-File','scripts/run_bilinear_grid.ps1'
   Get-Content data/processed/results/kge/bilinear_grid.log -Wait
=============================================================================
#>
[CmdletBinding()]
param(
    [string]$Python,
    [string[]]$Models  = @('DistMult', 'ComplEx'),
    [int[]]$Seeds      = @(42, 1, 7),
    [string[]]$Regimes = @('R0', 'R1', 'R2', 'R3'),
    [int]$Dim          = 64,
    [int]$Epochs       = 300,
    [switch]$Force
)

$RepoRoot = Split-Path -Parent $PSScriptRoot
Set-Location -Path $RepoRoot

if (-not $Python) {
    $venvPy = Join-Path $RepoRoot 'venv\Scripts\python.exe'
    if (Test-Path $venvPy) { $Python = $venvPy } else { $Python = 'python' }
}

$KgeDir = Join-Path $RepoRoot 'data\processed\results\kge'
New-Item -ItemType Directory -Force -Path $KgeDir | Out-Null
$Log    = Join-Path $KgeDir 'bilinear_grid.log'
$Status = Join-Path $KgeDir 'bilinear_grid.status'

$env:PYTORCH_CUDA_ALLOC_CONF = 'expandable_segments:True'

"STARTED $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Out-File -FilePath $Status -Encoding utf8
"  python  : $Python"                  | Out-File -FilePath $Status -Append -Encoding utf8
"  models  : $($Models -join ',')"     | Out-File -FilePath $Status -Append -Encoding utf8
"  regimes : $($Regimes -join ',')"    | Out-File -FilePath $Status -Append -Encoding utf8
"  seeds   : $($Seeds -join ',')"      | Out-File -FilePath $Status -Append -Encoding utf8
"  dim     : $Dim   epochs: $Epochs"   | Out-File -FilePath $Status -Append -Encoding utf8
"  recipe  : sLCWA, softplus (canonical bilinear objective), batch 16384, 16 negs, regularizer=none" |
    Out-File -FilePath $Status -Append -Encoding utf8

# One model at a time so a failure in one does not strand the other's cells.
foreach ($m in $Models) {
    "MODEL_START $m $(Get-Date -Format 'HH:mm:ss')" | Out-File -FilePath $Status -Append -Encoding utf8

    $pyArgs = @(
        'scripts/run_kge.py',
        '--models', $m,
        '--regimes') + $Regimes + @(
        '--seeds') + ($Seeds | ForEach-Object { "$_" }) + @(
        '--dim', "$Dim",
        '--epochs', "$Epochs",
        '--regularizer', 'none'
    )
    # Both bilinear models use their own published objective (logistic/softplus). ComplEx
    # REQUIRES it: with the regularizer off it still collapses under NSSA (loss pinned at
    # 4.5004 = the margin-9 zero-gap value from epoch 21 through 60, MRR 0.0797, AUROC
    # 0.500), whereas under softplus it reaches MRR 0.5081 / AUROC 0.860 at 60 epochs.
    # DistMult trains under either loss and takes softplus for consistency with ComplEx.
    if ($m -eq 'DistMult') { $pyArgs += @('--distmult-loss', 'softplus') }
    if ($m -eq 'ComplEx')  { $pyArgs += @('--complex-loss',  'softplus') }
    if ($Force)            { $pyArgs += '--force' }

    "  cmd: $Python $($pyArgs -join ' ')" | Out-File -FilePath $Status -Append -Encoding utf8
    & $Python @pyArgs 2>&1 | Out-File -FilePath $Log -Append -Encoding utf8
    $code = $LASTEXITCODE
    "MODEL_DONE $m $(Get-Date -Format 'HH:mm:ss') exit=$code" | Out-File -FilePath $Status -Append -Encoding utf8
}

"DONE $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Out-File -FilePath $Status -Append -Encoding utf8
