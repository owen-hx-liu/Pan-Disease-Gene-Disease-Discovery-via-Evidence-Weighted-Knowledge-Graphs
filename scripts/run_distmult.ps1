<#
=============================================================================
 run_distmult.ps1 -- launch the DistMult (third KGE baseline) headline grid.

 Adds DistMult (Yang et al., ICLR 2015) as a third KGE row alongside the frozen
 TransE + RotatE results, on the SAME splits and the SAME lib_eval harness, so it
 drops straight into Table 2 / Fig 1-2 when it lands. Writes per-run JSONs into
 data/processed/results/kge/ and re-merges kge_summary.json at the end (the
 existing TransE/RotatE JSONs are NOT touched -- no --force).

 Grid: DistMult x {R0,R1,R2,R3} x seeds {42,1,7}, dim 64, 300 epochs.
   12 cells, ~54 min/cell at ~10.5 s/epoch  ->  ~11 h wall on the RTX 5070 (8 GB).

 Recipe (see scripts/run_kge.py::model_config): SLCWA + softplus, 16 negs, at
 parity with the translational models. DistMult is bilinear, so -- exactly like
 ComplEx -- its unbounded scores saturate NSSA; softplus is the non-saturating
 loss. The 2- and 20-epoch smokes show this recipe COLLAPSES to the trivial
 optimum (loss pinned at ln 2 = 0.6931, MRR ~= chance, AUROC ~= 0.50), and the
 canonical LCWA cross-entropy alternative is intractable (OOM at 452k entities on
 8 GB). This full grid documents that collapse rigorously across seeds/regimes --
 the honest "caveat, not a clean third row" outcome. If any cell escapes chance it
 will show in the log; otherwise the caveat stands on 3 seeds x 4 regimes.

 Logs   -> data/processed/results/kge/distmult_grid.log   (all epochs, live-tailable)
 Status -> data/processed/results/kge/distmult_grid.status (STARTED/DONE + exit code)

 Usage (detached, survives this shell):
   Start-Process powershell -WindowStyle Hidden `
     -ArgumentList '-NoProfile','-ExecutionPolicy','Bypass','-File','scripts/run_distmult.ps1'
   Get-Content data/processed/results/kge/distmult_grid.log -Wait   # follow it
=============================================================================
#>
[CmdletBinding()]
param(
    [string]$Python,
    [int[]]$Seeds   = @(42, 1, 7),
    [string[]]$Regimes = @('R0', 'R1', 'R2', 'R3'),
    [int]$Dim       = 64,
    [int]$Epochs    = 300
)

$RepoRoot = Split-Path -Parent $PSScriptRoot
Set-Location -Path $RepoRoot

if (-not $Python) {
    $venvPy = Join-Path $RepoRoot 'venv\Scripts\python.exe'
    if (Test-Path $venvPy) { $Python = $venvPy } else { $Python = 'python' }
}

$KgeDir = Join-Path $RepoRoot 'data\processed\results\kge'
New-Item -ItemType Directory -Force -Path $KgeDir | Out-Null
$Log    = Join-Path $KgeDir 'distmult_grid.log'
$Status = Join-Path $KgeDir 'distmult_grid.status'

# Reduce fragmentation-driven OOM (harmless if the platform ignores it).
$env:PYTORCH_CUDA_ALLOC_CONF = 'expandable_segments:True'

"STARTED $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Out-File -FilePath $Status -Encoding utf8
"  python  : $Python"                                | Out-File -FilePath $Status -Append -Encoding utf8
"  models  : DistMult"                               | Out-File -FilePath $Status -Append -Encoding utf8
"  regimes : $($Regimes -join ',')"                  | Out-File -FilePath $Status -Append -Encoding utf8
"  seeds   : $($Seeds -join ',')"                    | Out-File -FilePath $Status -Append -Encoding utf8
"  dim     : $Dim   epochs: $Epochs"                 | Out-File -FilePath $Status -Append -Encoding utf8

$pyArgs = @(
    'scripts/run_kge.py',
    '--models', 'DistMult',
    '--regimes') + $Regimes + @(
    '--seeds') + ($Seeds | ForEach-Object { "$_" }) + @(
    '--dim', "$Dim",
    '--epochs', "$Epochs"
)

"GRID_START $(Get-Date -Format 'HH:mm:ss')" | Out-File -FilePath $Status -Append -Encoding utf8
& $Python @pyArgs 2>&1 | Out-File -FilePath $Log -Append -Encoding utf8
$code = $LASTEXITCODE
"DONE $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') exit=$code" | Out-File -FilePath $Status -Append -Encoding utf8
exit $code
