# below is the command you should paste into powershell to run this script
# This runs the executable as is. Ensure that make refresh has been run prior to execution 
# powershell -ExecutionPolicy Bypass -File .\run_ensemble.ps1 -NRuns 4

param(
    [int]$NRuns = 4,
    [string]$Exe = ".\VRES.exe",   
    [string]$Xml = ".\config\PhysiCell_settings.xml",
    [string]$OutRoot = ".\output_ensemble"
)

# create an ouput folder for the ensemble runs if they don't exist 
New-Item -ItemType Directory -Force -Path $OutRoot | Out-Null


# loop over each seed run
for ($seed = 1; $seed -le $NRuns; $seed++)
{
    Write-Host "=== Run $seed / $NRuns ==="

    # Update random_seed in XML
    $text = Get-Content $Xml -Raw

    # replace random seed value
    $text = [regex]::Replace(
        $text,
        "<random_seed[^>]*>.*?</random_seed>",
        "<random_seed type=`"int`" units=`"dimensionless`">$seed</random_seed>"
    )

    # write modified XML back to disk
    Set-Content -Path $Xml -Value $text -Encoding UTF8

    # Run simulation
    & $Exe

    # Move output into per run folder
    $runOut = Join-Path $OutRoot ("run_{0:D3}" -f $seed)
    New-Item -ItemType Directory -Force -Path $runOut | Out-Null

    if (Test-Path ".\output")
    {	# moves output folder created by physicell into a per run folder
        Move-Item -Force ".\output" (Join-Path $runOut "output")
    }
    else
    {
        Write-Warning "No output folder found for run $seed"
    }
}

Write-Host "All runs complete. Results in $OutRoot"
