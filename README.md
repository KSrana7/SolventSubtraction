# SolventSubtraction
Tensorial factorisation based approaches for solvent subtraction from spectroscopic data

# Approaches
2 different approaches namely direct and orthogonal is implemented here (for more detail refer to our paper). Appropriate approach can be chosen by selecting corresponding tag for "method" parameter in Solv_extr_*.m file type

# Pre-requisite
The added directory "ToolBoxes" is required to be included in the path.

# Data
The mentioned data in the paper is added in each subdirectory 

# Code execution flow
 First Run Data_pre_process.m for processing the data for artifacts and noise. The processed data is then utilised by the file Solv_extr_*.m. The best of solution is ensembled over N iteratios. If best solution is not converged properly, re-run again.
All other files are supporting fuction to execute these.


# Reference
If you find this repository useful in your research, please consider citing the following papers:
under review (will be updated soon)
