$ErrorActionPreference = "Stop"

$Root = "D:\QGenomix_mem1_step1"
$Plink = "$Root\tools\plink2\plink2.exe"
$Raw = "$Root\data\raw\1kg_phase3"
$Out = "$Root\data\interim\qc_by_chr"

New-Item -ItemType Directory -Force -Path $Raw | Out-Null
New-Item -ItemType Directory -Force -Path $Out | Out-Null

if (-not (Test-Path $Plink)) {
    throw "PLINK2 not found at $Plink"
}

$ResourcesUrl = "https://www.cog-genomics.org/plink/2.0/resources"
$ResourcesHtml = (Invoke-WebRequest -Uri $ResourcesUrl -UseBasicParsing).Content

function Get-DownloadUrl {
    param([string]$FileName)

    $pattern = 'href=["'']([^"'']*' + [regex]::Escape($FileName) + '[^"'']*)["'']'
    $match = [regex]::Matches($script:ResourcesHtml, $pattern) | Select-Object -First 1

    if (-not $match) {
        throw "Could not find download link for $FileName"
    }

    $url = [System.Net.WebUtility]::HtmlDecode($match.Groups[1].Value)

    if ($url.StartsWith("//")) {
        $url = "https:$url"
    } elseif ($url.StartsWith("/")) {
        $url = "https://www.cog-genomics.org$url"
    }

    if ($url -match "dropbox.com") {
        if ($url -match "dl=0") {
            $url = $url -replace "dl=0", "dl=1"
        } elseif ($url -notmatch "dl=1") {
            if ($url.Contains("?")) {
                $url = "$url&dl=1"
            } else {
                $url = "$url?dl=1"
            }
        }
    }

    return $url
}

function Download-IfMissing {
    param(
        [string]$FileName,
        [string]$OutFile
    )

    if ((Test-Path $OutFile) -and ((Get-Item $OutFile).Length -gt 1024)) {
        Write-Host "Already downloaded: $OutFile"
        return
    }

    $url = Get-DownloadUrl $FileName
    Write-Host "Downloading $FileName"

    & curl.exe --ssl-no-revoke -L --fail --retry 5 --retry-delay 10 --continue-at - --output $OutFile $url

    if ($LASTEXITCODE -ne 0) {
        throw "Download failed: $FileName"
    }

    if ((Get-Item $OutFile).Length -lt 1024) {
        throw "Downloaded file is suspiciously small: $OutFile"
    }
}

Write-Host "Checking PLINK2..."
& $Plink --version

Download-IfMissing "phase3_corrected.psam" "$Raw\phase3_corrected.psam"

foreach ($Chr in 1..22) {
    $Name = "chr$Chr" + "_phase3"

    $RawPrefix = "$Raw\$Name"
    $QcPrefix = "$Out\chr$Chr" + "_qc"
    $PrunePrefix = "$Out\chr$Chr" + "_prune"
    $FinalPrefix = "$Out\chr$Chr" + "_pruned"

    if ((Test-Path "$FinalPrefix.pgen") -and (Test-Path "$FinalPrefix.pvar.zst") -and (Test-Path "$FinalPrefix.psam")) {
        Write-Host "Already finished chr$Chr. Skipping."
        continue
    }

    Write-Host ""
    Write-Host "===== chr${Chr}: download ====="

    Download-IfMissing "$Name.pgen.zst" "$RawPrefix.pgen.zst"
    Download-IfMissing "$($Name)_noannot.pvar.zst" "$RawPrefix.pvar.zst"

    Copy-Item -Force "$Raw\phase3_corrected.psam" "$RawPrefix.psam"

    Write-Host "===== chr${Chr}: decompress .pgen ====="

    if (-not (Test-Path "$RawPrefix.pgen")) {
        & $Plink --zst-decompress "$RawPrefix.pgen.zst" "$RawPrefix.pgen"
        if ($LASTEXITCODE -ne 0) {
            throw "PLINK decompression failed for chr$Chr"
        }
    }

    Write-Host "===== chr${Chr}: quality control ====="

    & $Plink --pfile $RawPrefix vzs --maf 0.01 --hwe 1e-6 --geno 0.05 --mind 0.1 --max-alleles 2 --snps-only just-acgt --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 200 missing --rm-dup exclude-all --make-pgen vzs --out $QcPrefix
    if ($LASTEXITCODE -ne 0) {
        throw "QC failed for chr$Chr"
    }

    Write-Host "===== chr${Chr}: LD pruning list ====="

    & $Plink --pfile $QcPrefix vzs --indep-pairwise 50 5 0.2 --out $PrunePrefix
    if ($LASTEXITCODE -ne 0) {
        throw "Prune list failed for chr$Chr"
    }

    Write-Host "===== chr${Chr}: final pruned file ====="

    & $Plink --pfile $QcPrefix vzs --extract "$PrunePrefix.prune.in" --make-pgen vzs --out $FinalPrefix
    if ($LASTEXITCODE -ne 0) {
        throw "Final prune extraction failed for chr$Chr"
    }

    Write-Host "===== chr${Chr}: cleanup temporary files ====="

    Remove-Item -Force "$RawPrefix.pgen"

    Remove-Item -Force "$QcPrefix.pgen" -ErrorAction SilentlyContinue
    Remove-Item -Force "$QcPrefix.pvar.zst" -ErrorAction SilentlyContinue
    Remove-Item -Force "$QcPrefix.psam" -ErrorAction SilentlyContinue

    Write-Host "Finished chr$Chr"
}

Write-Host ""
Write-Host "All requested chromosomes finished."
Write-Host "Final files are in: $Out"

