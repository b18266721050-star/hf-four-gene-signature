param(
  [int]$MinFreeGB = 11,
  [int]$PollSec = 90,
  [int]$MaxHours = 12,
  [string]$Marker = "E:\another data file\GEObulk analysis project\_tmp\stage2_done.txt"
)
$base = "E:\another data file\GEObulk analysis project"
$exe  = "D:\R\R-4.6.1\bin\x64\Rscript.exe"
$env:SC_SUBSET = ""
$env:SC_RES = "0.6"
$log = "$base\_tmp\stage2_watcher.log"

function W($m) {
  "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')  $m" | Out-File -FilePath $log -Append -Encoding utf8
}

function Get-FreeGB {
  $os = Get-CimInstance Win32_OperatingSystem
  return [math]::Round($os.FreePhysicalMemory / 1MB, 2)
}

function Get-ForeignRCommitGB {
  $me = $PID
  $list = Get-CimInstance Win32_Process -Filter "Name='Rscript.exe' OR Name='Rterm.exe'" -ErrorAction SilentlyContinue
  $max = 0
  foreach ($p in $list) {
    if ($p.ProcessId -eq $me) { continue }
    $g = [math]::Round($p.PrivatePageCount / 1GB, 2)
    if ($g -gt $max) { $max = $g }
  }
  return $max
}

W "watcher v2 started (pid $PID); threshold = $MinFreeGB GB, 3 consecutive polls, poll = ${PollSec}s"
$start = Get-Date
$streak = 0
$started = $false

while (-not $started) {
  $free = Get-FreeGB
  $foreign = Get-ForeignRCommitGB
  if (($free -ge $MinFreeGB) -and ($foreign -lt 3)) {
    $streak = $streak + 1
    W "ok poll ($streak/3): free = $free GB, largest other R = $foreign GB"
  } else {
    $streak = 0
    W "waiting: free = $free GB, largest other R = $foreign GB"
  }
  if ($streak -ge 3) {
    W "conditions met -> launching pipeline"
    $started = $true
    break
  }
  if (((Get-Date) - $start).TotalHours -gt $MaxHours) {
    W "gave up after $MaxHours h"
    "timeout" | Out-File -FilePath $Marker -Encoding utf8
    exit 2
  }
  Start-Sleep -Seconds $PollSec
}

$scripts = @("41_sc_stage2_merge_norm_pca.R",
             "42_sc_stage2_harmony_cluster.R",
             "43_sc_stage2_markers.R")

foreach ($s in $scripts) {
  W "RUN  $s"
  $p = Start-Process -FilePath $exe -ArgumentList "`"$base\scripts\$s`"" -Wait -PassThru -WindowStyle Hidden
  $code = $p.ExitCode
  W "EXIT $s -> $code  (free = $(Get-FreeGB) GB)"
  if ($code -ne 0) {
    W "aborting pipeline"
    "fail:$s" | Out-File -FilePath $Marker -Encoding utf8
    exit 1
  }
}

"ok" | Out-File -FilePath $Marker -Encoding utf8
W "pipeline finished OK"
