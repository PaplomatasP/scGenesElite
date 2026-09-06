$ErrorActionPreference = "Stop"

$projectDirectory = Split-Path -Parent $PSScriptRoot
$rscriptCandidates = Get-ChildItem -LiteralPath "C:\Program Files\R" `
    -Filter "Rscript.exe" -Recurse -ErrorAction SilentlyContinue |
    Sort-Object FullName -Descending

if (-not $rscriptCandidates) {
    throw "Rscript.exe was not found under C:\Program Files\R. Install R first."
}

$rscript = $rscriptCandidates[0].FullName
Write-Host "Using $rscript"

$buildToolsCheck = @'
status <- if (
    requireNamespace("pkgbuild", quietly = TRUE) &&
    pkgbuild::has_build_tools(debug = FALSE)
) 0L else 1L
quit(status = status)
'@

& $rscript -e $buildToolsCheck
if ($LASTEXITCODE -ne 0) {
    Write-Host "Rtools was not detected. Installing the official Rtools45 required by R 4.6..."

    $rtoolsDirectoryUrl = "https://cloud.r-project.org/bin/windows/Rtools/rtools45/files/"
    $rtoolsPage = Invoke-WebRequest -Uri $rtoolsDirectoryUrl -UseBasicParsing
    $installerNames = @(
        [regex]::Matches(
            $rtoolsPage.Content,
            'href="(rtools45-[0-9]+-[0-9]+\.exe)"'
        ) | ForEach-Object { $_.Groups[1].Value } | Sort-Object -Unique
    )

    if (-not $installerNames) {
        throw "The current Rtools45 installer was not found at $rtoolsDirectoryUrl"
    }

    $installerName = $installerNames[-1]
    $installerUrl = $rtoolsDirectoryUrl + $installerName
    $temporaryRoot = [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
    $rtoolsInstaller = [IO.Path]::GetFullPath(
        (Join-Path $temporaryRoot "scgenes-$installerName")
    )

    if (-not $rtoolsInstaller.StartsWith(
        $temporaryRoot,
        [StringComparison]::OrdinalIgnoreCase
    )) {
        throw "Refusing to write the Rtools installer outside the temporary directory."
    }

    try {
        Write-Host "Downloading $installerUrl"
        Invoke-WebRequest -Uri $installerUrl -OutFile $rtoolsInstaller -UseBasicParsing

        $signature = Get-AuthenticodeSignature -LiteralPath $rtoolsInstaller
        if ($signature.Status -ne "Valid") {
            throw "The downloaded Rtools installer has an invalid Authenticode signature."
        }

        $process = Start-Process -FilePath $rtoolsInstaller `
            -ArgumentList @("/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART") `
            -WindowStyle Hidden -Wait -PassThru
        if ($process.ExitCode -ne 0) {
            throw "Rtools installation failed with exit code $($process.ExitCode)."
        }
    } finally {
        if (Test-Path -LiteralPath $rtoolsInstaller) {
            Remove-Item -LiteralPath $rtoolsInstaller -Force
        }
    }

    & $rscript -e $buildToolsCheck
    if ($LASTEXITCODE -ne 0) {
        throw "Rtools was installed but R still cannot find the build tools."
    }
}

Push-Location $projectDirectory
try {
    & $rscript "Install_Packages\install_all.R"
    if ($LASTEXITCODE -ne 0) {
        throw "R package installation failed with exit code $LASTEXITCODE."
    }
} finally {
    Pop-Location
}
