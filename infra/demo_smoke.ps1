param(
    [string]$ApiBase = "http://localhost:18000",
    [string]$WebBase = "http://localhost:13100",
    [string]$ProjectId = "demo-project"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..")).Path
$composeFile = Join-Path $scriptDir "docker-compose.yml"
$outputDir = Join-Path $scriptDir "demo-output"

New-Item -ItemType Directory -Path $outputDir -Force | Out-Null

Write-Host "Starting docker compose stack..."
Push-Location $repoRoot
try {
    docker compose -f $composeFile up -d --build | Out-Host
} finally {
    Pop-Location
}

Write-Host "Waiting for API readiness at $ApiBase/health ..."
$apiReady = $false
for ($i = 0; $i -lt 60; $i++) {
    try {
        $health = Invoke-RestMethod -Method Get -Uri "$ApiBase/health" -TimeoutSec 5
        if ($health.status -eq "ok") {
            $apiReady = $true
            break
        }
    } catch {
        # retry
    }
    Start-Sleep -Seconds 2
}
if (-not $apiReady) {
    throw "API was not ready within 120 seconds."
}

$loginCandidates = @(
    @{ email = "admin@example.com"; password = "changeme" },
    @{ email = "admin@example.com"; password = "password" }
)

$token = ""
foreach ($candidate in $loginCandidates) {
    try {
        $body = $candidate | ConvertTo-Json -Compress
        $login = Invoke-RestMethod -Method Post -Uri "$ApiBase/auth/login" -ContentType "application/json" -Body $body
        if ($login.access_token) {
            $token = $login.access_token
            break
        }
    } catch {
        # try next credential
    }
}

if (-not $token) {
    throw "Unable to login as admin@example.com with known demo passwords."
}

$authHeader = "Authorization: Bearer $token"
$jsonHeaders = @{
    Authorization = "Bearer $token"
    "Content-Type" = "application/json"
}
$authHeaders = @{
    Authorization = "Bearer $token"
}

$vcfPath = Join-Path $scriptDir "example.vcf"
$exprPath = Join-Path $scriptDir "expr.tsv"
$fastqPath = Join-Path $scriptDir "raw_demo.fastq"
$bamPath = Join-Path $scriptDir "raw_demo.bam"
$sampleDataDir = Join-Path $scriptDir "sample_data"
$fastaPath = Join-Path $sampleDataDir "demo_reference.fasta"
$gtfPath = Join-Path $sampleDataDir "demo_annotation.gtf"
$pairR1Path = Join-Path $sampleDataDir "demo_reads_R1.fastq"
$pairR2Path = Join-Path $sampleDataDir "demo_reads_R2.fastq"
$sampleSheetPath = Join-Path $sampleDataDir "demo_samplesheet.tsv"
if (-not (Test-Path $vcfPath)) {
    throw "Missing demo VCF: $vcfPath"
}
if (-not (Test-Path $exprPath)) {
    throw "Missing demo expression table: $exprPath"
}
if (-not (Test-Path $fastqPath)) {
    throw "Missing raw FASTQ demo file: $fastqPath"
}
if (-not (Test-Path $bamPath)) {
    throw "Missing raw BAM demo file: $bamPath"
}
if (-not (Test-Path $fastaPath)) {
    throw "Missing demo FASTA file: $fastaPath"
}
if (-not (Test-Path $gtfPath)) {
    throw "Missing demo GTF file: $gtfPath"
}
if (-not (Test-Path $pairR1Path)) {
    throw "Missing demo FASTQ R1 file: $pairR1Path"
}
if (-not (Test-Path $pairR2Path)) {
    throw "Missing demo FASTQ R2 file: $pairR2Path"
}
if (-not (Test-Path $sampleSheetPath)) {
    throw "Missing demo sample sheet: $sampleSheetPath"
}

Write-Host "Uploading raw FASTQ dataset..."
$uploadFastqRaw = & curl.exe -sS -X POST "$ApiBase/datasets/upload?project_id=$ProjectId" `
    -H $authHeader `
    -F "file=@$fastqPath;type=text/plain"
if ($LASTEXITCODE -ne 0) {
    throw "Raw FASTQ dataset upload command failed."
}
$uploadFastq = $uploadFastqRaw | ConvertFrom-Json
if (-not $uploadFastq.dataset.id) {
    throw "Raw FASTQ dataset upload failed: $uploadFastqRaw"
}

Write-Host "Uploading raw BAM dataset..."
$uploadBamRaw = & curl.exe -sS -X POST "$ApiBase/datasets/upload?project_id=$ProjectId" `
    -H $authHeader `
    -F "file=@$bamPath;type=application/octet-stream"
if ($LASTEXITCODE -ne 0) {
    throw "Raw BAM dataset upload command failed."
}
$uploadBam = $uploadBamRaw | ConvertFrom-Json
if (-not $uploadBam.dataset.id) {
    throw "Raw BAM dataset upload failed: $uploadBamRaw"
}

Write-Host "Uploading bio sample bundle (FASTA/GTF/paired FASTQ/sample sheet)..."
$sampleUploads = @(
    @{ path = $fastaPath; contentType = "text/plain" },
    @{ path = $gtfPath; contentType = "text/plain" },
    @{ path = $pairR1Path; contentType = "text/plain" },
    @{ path = $pairR2Path; contentType = "text/plain" },
    @{ path = $sampleSheetPath; contentType = "text/tab-separated-values" }
)
foreach ($item in $sampleUploads) {
    $uploadRaw = & curl.exe -sS -X POST "$ApiBase/datasets/upload?project_id=$ProjectId" `
        -H $authHeader `
        -F "file=@$($item.path);type=$($item.contentType)"
    if ($LASTEXITCODE -ne 0) {
        throw "Sample dataset upload command failed: $($item.path)"
    }
    $uploadObj = $uploadRaw | ConvertFrom-Json
    if (-not $uploadObj.dataset.id) {
        throw "Sample dataset upload failed: $uploadRaw"
    }
}

Write-Host "Uploading demo VCF..."
$ingestVariantRaw = & curl.exe -sS -X POST "$ApiBase/variants/ingest?project_id=$ProjectId" `
    -H $authHeader `
    -F "file=@$vcfPath;type=text/vcf"
if ($LASTEXITCODE -ne 0) {
    throw "VCF upload command failed."
}
$ingestVariant = $ingestVariantRaw | ConvertFrom-Json
if (-not $ingestVariant.ingested -or [int]$ingestVariant.ingested -lt 1) {
    throw "VCF upload returned no ingested variants: $ingestVariantRaw"
}

Write-Host "Uploading demo expression table..."
$ingestExprRaw = & curl.exe -sS -X POST "$ApiBase/omics/expr/upload?project_id=$ProjectId" `
    -H $authHeader `
    -F "file=@$exprPath;type=text/tab-separated-values"
if ($LASTEXITCODE -ne 0) {
    throw "Expression upload command failed."
}
$ingestExpr = $ingestExprRaw | ConvertFrom-Json
if (-not $ingestExpr.ingested -or [int]$ingestExpr.ingested -lt 1) {
    throw "Expression upload returned no rows: $ingestExprRaw"
}

Write-Host "Running evidence ranking..."
$evidence = Invoke-RestMethod -Method Post `
    -Uri "$ApiBase/runs/evidence?chr=chr1&start=1&end=1000&top_n=5&project_id=$ProjectId" `
    -Headers $authHeaders

Write-Host "Running causal scoring..."
$causal = Invoke-RestMethod -Method Post `
    -Uri "$ApiBase/causal/score?chr=chr1&start=1&end=1000&top_n=5&project_id=$ProjectId" `
    -Headers $authHeaders

Write-Host "Generating research directions..."
$researchBody = @{
    disease = "Alzheimer disease"
    genes = @("GENE1", "GENE2")
    top_n = 2
} | ConvertTo-Json -Compress
$researchDirections = Invoke-RestMethod -Method Post -Uri "$ApiBase/research/directions" -Headers $jsonHeaders -Body $researchBody

$evidenceRunId = $evidence.run.run_id
$reportMdPath = Join-Path $outputDir "report-$evidenceRunId.md"
$reportHtmlPath = Join-Path $outputDir "report-$evidenceRunId.html"

Write-Host "Exporting markdown/html reports..."
Invoke-WebRequest -Uri "$ApiBase/report/$ProjectId/$evidenceRunId" -Headers $authHeaders -OutFile $reportMdPath
Invoke-WebRequest -Uri "$ApiBase/report/$ProjectId/${evidenceRunId}?format=html" -Headers $authHeaders -OutFile $reportHtmlPath

$webStatus = "unreachable"
$webCandidates = @($WebBase)
if ($WebBase -like "*:13000*") {
    $webCandidates += "http://localhost:13100"
} elseif ($WebBase -like "*:13100*") {
    $webCandidates += "http://localhost:13000"
}

foreach ($candidate in $webCandidates) {
    for ($i = 0; $i -lt 20; $i++) {
        try {
            $resp = Invoke-WebRequest -Uri "$candidate/" -Method Get -TimeoutSec 10 -UseBasicParsing
            if ($resp.StatusCode -ge 200 -and $resp.StatusCode -lt 500) {
                $webStatus = "ok ($($resp.StatusCode)) via $candidate"
                break
            }
        } catch {
            # retry next round
        }
        Start-Sleep -Seconds 2
    }
    if ($webStatus -like "ok*") {
        break
    }
}

Write-Host ""
Write-Host "Demo complete."
Write-Host "Project: $ProjectId"
Write-Host "Evidence run: $evidenceRunId"
Write-Host "Causal run: $($causal.run.run_id)"
Write-Host "Top evidence genes: $((($evidence.result.ranked_genes | Select-Object -First 3).gene) -join ', ')"
Write-Host "Research hypotheses: $($researchDirections.hypotheses.Count)"
Write-Host "Report (md): $reportMdPath"
Write-Host "Report (html): $reportHtmlPath"
Write-Host "Web locus status: $webStatus"
