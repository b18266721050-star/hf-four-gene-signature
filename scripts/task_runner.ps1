param(
  [Parameter(Mandatory=$true)][string]$RScript,
  [string]$Subset = "",
  [string]$Res = "0.6",
  [Parameter(Mandatory=$true)][string]$Marker
)
$ErrorActionPreference = "Continue"
if ($Subset -ne "") { $env:SC_SUBSET = $Subset } else { $env:SC_SUBSET = "" }
$env:SC_RES = $Res
$exe = "D:\R\R-4.6.1\bin\x64\Rscript.exe"
& $exe $RScript
$code = $LASTEXITCODE
"exit=$code" | Out-File -FilePath $Marker -Encoding utf8
