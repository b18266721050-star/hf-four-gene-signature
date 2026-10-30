param(
  [string]$Base = "E:\another data file\GEObulk analysis project"
)
$exe = "D:\R\R-4.6.1\bin\x64\Rscript.exe"
$log = "$Base\_tmp\stage2_trigger.log"
function W($m) {
  "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')  $m" | Out-File -FilePath $log -Append -Encoding utf8
}

$done = "$Base\_tmp\stage2_done.txt"
$lock = "$Base\_tmp\stage2_running.lock"

if (Test-Path $done) { exit 0 }

# a pipeline still running?
if (Test-Path $lock) {
  $lockPid = (Get-Content $lock -TotalCount 1)
  if ($lockPid -and (Get-Process -Id ([int]$lockPid) -ErrorAction SilentlyContinue)) { exit 0 }
  Remove-Item $lock -Force -ErrorAction SilentlyContinue
}

# is one of the stage scripts already running?
$rs = @(Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'" -ErrorAction SilentlyContinue)
foreach ($r in $rs) {
  if ($r.CommandLine -and ($r.CommandLine -like '*41_sc_stage2_merge_norm_pca.R' -or
                            $r.CommandLine -like '*42_sc_stage2_harmony_cluster.R' -or
                            $r.CommandLine -like '*43_sc_stage2_markers.R')) { exit 0 }
}

$os = Get-CimInstance Win32_OperatingSystem
$free = [math]::Round($os.FreePhysicalMemory / 1MB, 2)
$maxOther = 0
foreach ($r in $rs) {
  $g = [math]::Round($r.PrivatePageCount / 1GB, 2)
  if ($g -gt $maxOther) { $maxOther = $g }
}

if (($free -lt 7.5) -or ($maxOther -ge 3)) {
  W "skip: free = $free GB, largest R = $maxOther GB"
  exit 0
}

W "GO: free = $free GB, largest R = $maxOther GB"
$env:SC_SUBSET = ""
$env:SC_RES = "0.6"
$runner = "$Base\scripts\run_pipeline_stage2.ps1"
$p = Start-Process -FilePath "powershell.exe" `
      -ArgumentList "-NoProfile -ExecutionPolicy Bypass -File `"$runner`" -Base `"$Base`"" `
      -WindowStyle Hidden -PassThru
$p.Id | Out-File -FilePath $lock -Encoding utf8
W "pipeline launched, wrapper pid = $($p.Id)"
