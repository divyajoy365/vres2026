# vres2026

- To run the simulation, place the VRES folder into the sample_projects directory inside PhysiCell.
-  Note: I am using PhysiCell-1.14.2.
The simulation can be run as usual from the Anaconda Prompt:
- Run 'vres refresh'
- Then run 'vres'

run_ensemble.ps1 is a PowerShell script to run multitple sims with different seeds. 
- Open PowerShell
- Navigate to the main PhysiCell directory
- Ensure vres refresh has already been run in the Anaconda Prompt (as usual)
- In the PowerShell prompt, run the script using 'powershell -ExecutionPolicy Bypass -File .\run_ensemble.ps1 -NRuns 4'
- This will run 4 simulations. You can change the number to changed the number of simulations. 

esemble_output.om is a MATLAB script to visualise outputs. 
- When the PowerShell script is used to run multiple simulations, it creates an output_ensemble folder.
- You need to point the MATLAB script to this folder: Update line 6 of the script to the correct path on your computer
Once the PowerShell script has finished running:
- Run the MATLAB script to visualise: Total rod counts, Total cell counts, Number of active vs inactive cells
