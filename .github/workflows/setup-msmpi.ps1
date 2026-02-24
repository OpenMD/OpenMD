$ErrorActionPreference = "Stop"

$version = "10.1.3"

$rooturl = "https://download.microsoft.com/download/"
$hashurl = @{
 "5.0.0" = "3/7/6/3764A48C-5C4E-4E4D-91DA-68CED9032EDE"
 "6.0.0" = "6/4/A/64A7852A-A8C3-476D-908C-30501F761DF3"
 "7.0.0" = "D/7/B/D7BBA00F-71B7-436B-80BC-4D22F2EE9862"
 "7.1.0" = "E/8/A/E8A080AF-040D-43FF-97B4-065D4F220301"
 "8.0.0" = "B/2/E/B2EB83FE-98C2-4156-834A-E1711E6884FB"
 "8.1.0" = "D/B/B/DBB64BA1-7B51-43DB-8BF1-D1FB45EACF7A"
 "9.0.0" = "2/E/C/2EC96D7F-687B-4613-80F6-E10F670A2D97"
 "9.0.1" = "4/A/6/4A6AAED8-200C-457C-AB86-37505DE4C90D"
"10.0.0" = "A/E/0/AE002626-9D9D-448D-8197-1EA510E297CE"
"10.1.1" = "2/9/e/29efe9b1-16d7-4912-a229-6734b0c4e235"
"10.1.2" = "a/5/2/a5207ca5-1203-491a-8fb8-906fd68ae623"
"10.1.3" = "7/2/7/72731ebb-b63c-4170-ade7-836966263a8f"
}
$baseurl = $rooturl + $hashurl[$version]

$tempdir    = $Env:RUNNER_TEMP
$msmpisdk   = Join-Path $tempdir msmpisdk.msi
$msmpisetup = Join-Path $tempdir msmpisetup.exe

function Download-File($url, $filename) {
  foreach ($i in 1..5) {
    try {
      Write-Host "Downloading ${url}"
      Invoke-WebRequest $url -OutFile $filename
      return
    }
    catch {
      $message = $_
      Write-Warning "${message}"
      Write-Host "Download failed, retrying ..."
      Start-Sleep -Seconds $i
    }
  }
  throw "Failed to download from ${url}"
  return $null
}

Write-Host "Downloading Microsoft MPI $version"
Download-File "$baseurl/msmpisdk.msi"   $msmpisdk
Download-File "$baseurl/msmpisetup.exe" $msmpisetup

Write-Host "Installing Microsoft MPI $version"
Start-Process msiexec.exe -ArgumentList "/quiet /passive /qn /i $msmpisdk" -Wait
Start-Process $msmpisetup -ArgumentList "-unattend" -Wait

if ($Env:GITHUB_ENV) {
  Write-Host 'Adding environment variables to $GITHUB_ENV'
  $envlist = @("MSMPI_BIN", "MSMPI_INC", "MSMPI_LIB32", "MSMPI_LIB64")
  foreach ($name in $envlist) {
    $value = [Environment]::GetEnvironmentVariable($name, "Machine")
    Write-Host "$name=$value"
    Add-Content $Env:GITHUB_ENV "$name=$value"
  }
}

if ($Env:GITHUB_PATH) {
  Write-Host 'Adding $MSMPI_BIN to $GITHUB_PATH'
  $MSMPI_BIN = [Environment]::GetEnvironmentVariable("MSMPI_BIN", "Machine")
  Add-Content $Env:GITHUB_PATH $MSMPI_BIN
}
