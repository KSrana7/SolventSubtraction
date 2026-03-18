# SolventSubtraction
Tensorial factorisation based approaches for solvent subtraction from spectroscopic data

# Approaches
2 different approaches namely direct and orthogonal is implemented here (for more detail refer to our paper). Appropriate approach can be chosen by selecting corresponding tag for "method" parameter in Solv_extr_*.m file type

# Pre-requisite
The added directory "ToolBoxes" is required to be included in the path.

# Data
The mentioned data in the paper is added in each subdirectory 

# Code execution flow
First Run Data_pre_process.m for processing the data for artifacts and noise. The processed data is then utilised by the file Solv_extr_*.m. The best of solution is ensembled over N iteratios. If best solution is not converged properly, re-run again or tune the paramters.
All other files are supporting fuction to execute these.
For Bitumen data, run solv_extrac_v0.m file only.


# Reference
If you find this repository useful in your research, please consider citing the following papers:
'''
@article{singh2025scalable,
  title={A Scalable and Generalizable Method to Minimize Solvent Interference in Identification of Chemical Reaction Networks from Spectroscopic Data},
  author={Singh, Kuldeep and Srinivasan, Karthik and Sun, Ziting and Liu, Jing and Prasad, Vinay},
  journal={Journal of Chemical Information and Modeling},
  volume={65},
  number={19},
  pages={10068--10092},
  year={2025},
  publisher={ACS Publications}
}
'''
