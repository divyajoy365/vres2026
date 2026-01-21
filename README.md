# vres2026

VRES folder
- To run the simulation, place the VRES folder into the sample_projects directory inside PhysiCell.
-  Note: I am using PhysiCell-1.14.2.
The simulation can be run as usual from the Anaconda Prompt:
- Run 'make refresh'
- Then run 'vres'

run_ensemble.ps1 is a PowerShell script to run multitple simulations with different seeds. To run it:
- Open Windows PowerShell (you can do so by typing it into the Start menu)
- Navigate to the main PhysiCell directory
- Ensure 'vres refresh' has already been run in the Anaconda Prompt (as usual)
- In the PowerShell prompt, run the script using 'powershell -ExecutionPolicy Bypass -File .\run_ensemble.ps1 -NRuns 4'
- This will run 4 simulations. You can change the number to the desired number of simulations. 

esemble_output.m is a MATLAB script to visualise outputs. 
- When the PowerShell script is used to run multiple simulations, it creates an output_ensemble folder to store the output of the different simulations.
- You need to point the MATLAB script to this folder: Update line 6 of the script to the correct path on your computer
- Once the PowerShell script has finished running, run the MATLAB script to visualise: total rod cell counts, total cell counts, number of active vs inactive cells
