Project: CPIDM
Desciprtion
- R scripts written for DE analysis
- main.R is the main script for the analysis and reads other scripts used in the analysis
- scripts for 3D PCA is not included in main.R
- scripts for 3D PCA also need manual configuration about the input data to be used
	- `choice_ntop = 0` to use DEGs
	- `choice_ntop = 1` to use ntop
- scripts are generally dependent on each other, so a full run (at least to a point) is required for modular process of parts of main.R

