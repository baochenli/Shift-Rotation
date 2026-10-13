$ErrorActionPreference = 'Stop'
$latexExe = 'D:\texlive\2023\texlive\2023\bin\windows\pdflatex.exe'
Push-Location $PSScriptRoot
try {
    for ($pass = 1; $pass -le 2; $pass++) {
        & $latexExe -interaction=nonstopmode -halt-on-error -file-line-error 'ShiftRotation_derivation.tex'
        if ($LASTEXITCODE -ne 0) { throw "pdfLaTeX failed on pass $pass" }
    }
} finally { Pop-Location }
