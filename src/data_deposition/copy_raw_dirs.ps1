<#
.SYNOPSIS
    Copy the raw files listed in a TSV into a single flat PRIDE staging folder.

.DESCRIPTION
    PowerShell has no rsync; robocopy is the equivalent. Only *.raw files from
    the top level of each source directory are copied, and they all land
    directly in -Destination - no per-experiment subfolders, no subdirectories
    of any kind. Raw files nested inside subdirectories of a source are skipped,
    and the script warns about each one it finds.

    QC wash injections (Wash*.raw) are excluded by default - they are not part
    of the submission and are the main source of filename collisions. Override
    with -IncludeWash, or set your own patterns with -ExcludeFile.

    Because everything is flattened into one folder, two sources holding a
    same-named .raw file would overwrite each other. The script scans for that
    up front and refuses to copy anything until it is resolved (or -Force is
    passed). Sources listed twice in the TSV are de-duplicated, not a conflict.

    Robocopy is restartable and skips files that already match on size+timestamp,
    so re-running resumes rather than re-copying. Nothing is ever deleted: no
    /MIR, no /PURGE, no /MOVE.

.EXAMPLE
    # Dry run first - lists what would be copied, transfers nothing.
    .\src\data_deposition\copy_raw_dirs.ps1 -WhatIf

.EXAMPLE
    .\src\data_deposition\copy_raw_dirs.ps1
#>

[CmdletBinding()]
param(
    [string]$Tsv = (Join-Path $PSScriptRoot '..\..\raw_file_locations_normalized.tsv'),
    [string]$Destination = "D:\Data HS\macrophage-pride-submission\raw-metadata-searches",
    [string]$LogDir = ".\robocopy-logs",
    # List only, copy nothing.
    [switch]$WhatIf,
    # Copy even if filenames collide between sources (last writer wins).
    [switch]$Force,
    # QC wash injections - not part of the submission, and the only real source
    # of filename collisions (plain "Wash.raw" appears in most acquisition
    # folders). Pass -IncludeWash to copy them anyway.
    [string[]]$ExcludeFile = @("Wash*.raw"),
    [switch]$IncludeWash
)

$ErrorActionPreference = "Stop"

$rows = Import-Csv -Path $Tsv -Delimiter "`t" |
    Where-Object { $_.raw_file_location -and $_.raw_file_location.Trim() }

if (-not $rows) { throw "No populated raw_file_location rows found in $Tsv" }

if ($IncludeWash) { $ExcludeFile = @() }

# Several experiments can point at the same acquisition folder; copy each once.
$sources = $rows |
    Group-Object { $_.raw_file_location.Trim().TrimEnd('\') } |
    ForEach-Object {
        [pscustomobject]@{
            Path        = $_.Name
            Experiments = ($_.Group.experiment -join ", ")
        }
    }

Write-Host "$($rows.Count) rows -> $($sources.Count) unique source directories" -ForegroundColor Cyan
Write-Host ""

# ---------------------------------------------------------------- pre-flight
$missing   = @()
$nestedAll = @()
$excluded  = 0
$byName    = @{}

foreach ($s in $sources) {
    if (-not (Test-Path -LiteralPath $s.Path)) {
        Write-Warning "MISSING source ($($s.Experiments)): $($s.Path)"
        $missing += $s
        continue
    }

    # .ProviderPath, not .Path: for UNC paths .Path comes back provider-qualified
    # ("Microsoft.PowerShell.Core\FileSystem::\\server\share\..."), which would
    # never match $f.DirectoryName and mark every file as nested.
    $root = (Resolve-Path -LiteralPath $s.Path).ProviderPath.TrimEnd('\')

    foreach ($f in Get-ChildItem -LiteralPath $s.Path -Filter *.raw -Recurse -File -ErrorAction SilentlyContinue) {
        # Mirror robocopy's /XF so the counts and collision check match reality.
        $skip = $false
        foreach ($pat in $ExcludeFile) { if ($f.Name -like $pat) { $skip = $true; break } }
        if ($skip) { $excluded++; continue }

        if ($f.DirectoryName -ne $root) {
            # /LEV:1 would drop these silently, so surface them instead.
            $nestedAll += $f
            continue
        }
        if (-not $byName.ContainsKey($f.Name)) { $byName[$f.Name] = @() }
        $byName[$f.Name] += $root
    }
}

if ($nestedAll) {
    Write-Warning "$($nestedAll.Count) .raw file(s) live in subdirectories and will NOT be copied:"
    $nestedAll | ForEach-Object { Write-Warning "    $($_.FullName)" }
    Write-Host ""
}

$collisions = $byName.GetEnumerator() |
    Where-Object { ($_.Value | Select-Object -Unique).Count -gt 1 } |
    Sort-Object Key

if ($collisions) {
    Write-Host "FILENAME COLLISIONS - the same .raw name comes from more than one source:" -ForegroundColor Red
    foreach ($c in $collisions) {
        Write-Host "  $($c.Key)" -ForegroundColor Red
        $c.Value | Select-Object -Unique | ForEach-Object { Write-Host "      $_" }
    }
    Write-Host ""
    if (-not $Force) {
        throw ("$($collisions.Count) filename collision(s); flattening would overwrite data. " +
               "Resolve them, or re-run with -Force to accept last-writer-wins.")
    }
    Write-Warning "-Force given: proceeding despite collisions."
}

$totalFiles = ($byName.Values | ForEach-Object { $_.Count } | Measure-Object -Sum).Sum
if ($excluded) {
    Write-Host "$excluded .raw file(s) excluded by -ExcludeFile ($($ExcludeFile -join ', '))" -ForegroundColor DarkGray
}
Write-Host "$totalFiles top-level .raw file(s) to copy into $Destination" -ForegroundColor Cyan
Write-Host ""

# --------------------------------------------------------------------- copy
New-Item -ItemType Directory -Force -Path $Destination, $LogDir | Out-Null

$failed = @()

foreach ($s in $sources) {
    if ($missing -contains $s) { continue }

    $leaf = Split-Path $s.Path -Leaf
    $log  = Join-Path $LogDir ("{0}.log" -f $leaf)

    Write-Host "==> $($s.Experiments)" -ForegroundColor Cyan
    Write-Host "    $($s.Path)"

    # *.raw  copy only Thermo raw files
    # /LEV:1 top level only - no recursion, so no subdirectories are created
    #        in the destination and every file lands directly in $Destination
    # /COPY:DAT data+attrs+timestamps           /R:2 /W:5 brief retries on a flaky share
    # /MT:8  multithreaded                      /NP no per-file progress spam
    # /XX    don't report "extra" destination files. Every source copies into the
    #        same flat folder, so each run would otherwise list all 13 other
    #        experiments' files as extras. Reporting only - /XX is the default
    #        behaviour for what gets copied, and nothing is ever deleted.
    # /L     list only (dry run)
    $rcArgs = @($s.Path, $Destination, "*.raw", "/LEV:1", "/COPY:DAT", "/R:2", "/W:5",
                "/MT:8", "/NP", "/XX", "/TEE", "/LOG+:$log")
    if ($ExcludeFile) { $rcArgs += @("/XF") + $ExcludeFile }
    if ($WhatIf)      { $rcArgs += "/L" }

    & robocopy.exe @rcArgs
    $code = $LASTEXITCODE

    # Robocopy: 0-7 are success/informational, >=8 means at least one failure.
    if ($code -ge 8) {
        Write-Warning "robocopy exit $code for $($s.Experiments) - see $log"
        $failed += [pscustomobject]@{ Experiments = $s.Experiments; Reason = "robocopy exit $code"; Source = $s.Path }
    }
}

# ------------------------------------------------------------------ summary
Write-Host ""
if ($missing) {
    Write-Host "Sources not found:" -ForegroundColor Yellow
    $missing | Format-Table Experiments, Path -AutoSize
}
if ($failed) {
    Write-Host "Robocopy failures:" -ForegroundColor Yellow
    $failed | Format-Table -AutoSize
}
if (-not $missing -and -not $failed) {
    Write-Host "All $($sources.Count) source directories copied cleanly." -ForegroundColor Green
}

if (-not $WhatIf) {
    $onDisk = @(Get-ChildItem -LiteralPath $Destination -Filter *.raw -File).Count
    Write-Host "$onDisk .raw file(s) now in $Destination (expected $totalFiles)" -ForegroundColor Cyan
}
