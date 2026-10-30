param(
  [string]$Base = "E:\another data file\GEObulk analysis project"
)
$exe = "D:\R\R-4.6.1\bin\x64\Rscript.exe"
$env:SC_SUBSET = ""
$env:SC_RES = "0.6"
$done = "$Base\_tmp\done_43_full.txt"
"running" | Out-File -FilePath $done -Encoding utf8
$p = Start-Process -FilePath $exe -ArgumentList "`"$Base\scripts\43_sc_stage2_markers.R`"" -Wait -PassThru -WindowStyle Hidden
"exit=$($p.ExitCode)" | Out-File -FilePath $done -Encoding utf8
