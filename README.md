# FANGS
Flexible Algorithm for Novel Gene set Simulation (FANGS) is a semi-synthetic simulation software which simulates data based on a user defined gene expression dataset, to test and assist in the development of Gene Set Analysis (GSA) methods.

Below is a description of the files available and there functionality.

caseControl_centered_final.rda
- True expression dataset as a NxM matrix with each row a sample (N) and each column a gene (M). The first column is the class label, where 0 is normal tissue and 1 is a prostate cancer tissue. 
  The column names are the name of each of gene.

bootstrapSampling.R
- main R script to run the simulation software - run this script with the appropriate inputs (see comments in the script)

changeExprs.R
- contains key functions for the simulation software

bootstrap_paramaters.csv
- software parameters to control the reproducibility of the results

c6_KRAS.PROSTATE_UP.V1_DN_genesChange.txt
- comma separated file that contains the gene names of the genes contained in the KRAS pathway

h_TGF_BETA_SIGNALING_genesChange.txt
- comma separated file that contains the gene names of the genes contained in the TGF-Beta pathway

Note: Such txt files as above need to be created when extending this simulation software to target other pathways.


runBootstraping.sh
- example bash shell script to run the R file via command line on a unix based operating system.

startBoostraping.sh
- example bash shell script to run runBootstraping.sh on the brccluster available for BRC affiliates only at NC State University.

exampleResults directory
- example of the result files that are created by the simulation; with the simulationParameters.txt file describing the simulation parameters
