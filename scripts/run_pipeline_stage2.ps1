param(
  [string]$Base = "E:\another data file\GEObulk analysis project"
)
$exe = "D:\R\R-4.6.1\bin\x64\Rscript.exe"
$log = "$Base\_tmp\stage2_run.log"
$done = "$Base\_tmp\stage2_done.txt"
function W($m) {
  "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')  $m" | Out-File -FilePath $log -Append -Encoding utf8
}

$env:SC_SUBSET = ""
$env:SC_RES = "0.6"
W "pipeline wrapper start (pid $PID)"

$scripts = @("41_sc_stage2_merge_norm_pca.R",
             "42_sc_stage2_harmony_cluster.R",
             "43_sc_stage2_markers.R")

foreach ($s in $scripts) {
  W "RUN $s"
  $p = Start-Process -FilePath $exe -ArgumentList "`"$Base\scripts\$s`"" -Wait -PassThru -WindowStyle Hidden
  W "EXIT $s -> $($p.ExitCode)"
  if ($p.ExitCode -ne 0) {
    W "ABORT"
    "fail:$s" | Out-File -FilePath $done -Encoding utf8
    exit 1
  }
}
"ok" | Out-File -FilePath $done -Encoding utf8
W "pipeline finished OK"
