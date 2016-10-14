##### Simulation software based on a real dataset to test Gene Set Analysis (GSA) methods. The simulated data sets are based 
##### on a real dataset (user defined) where specific pathways (user defined) are differentially expressed.
##### Version: 1.0 - September 2016
##### Author: Ravi Mathur (rmathur2@ncsu.edu), Motsinger-Reif Lab, Bioinformatics Research Center, NCSU


#### Parameters files to define the simulation seed and common directory of scripts
fileArgs <- scan("bootstrap_parameters.csv", sep=',', what="")
seed <- as.integer(fileArgs[which(fileArgs == 'seed')+2])
set.seed(seed)
theDir <- fileArgs[which(fileArgs == 'commDir')+2]

#### Import/Source the fuctions used to differentially express the genes contained in the targetted pathway.
source(paste0(theDir, "/changeExprs.R"))



#### Command line arguments 
args <- commandArgs(trailingOnly = TRUE)
whichData <- args[1] # 'prostate' cancer dataset is currently implemented and tested
changeType <- args[2] # 'geneSet' is currently implemented and tested
changeAmount <- as.numeric(args[3]) # tau value
changeProp <- as.numeric(args[4]) # pi value
geneSetType <- args[5] # MSigDB category of the pathway that is being targetted
geneSet <- args[6] # name of the MSigDB pathway that is being targetted
bN <- as.numeric(args[7]) # Number of bootstrap iterations
changeDirection <- args[8] # 'OverExpress' or 'UnderExpress' are currently implemented and tested
dataPath <- args[9] #Full path and filename of expression data r data file
dataObject.name <- args[10] #Object name of expression data as saved in the r data file
resultsDir <- args[11] #Directory where all results should be saved
geneChange.file <- args[12] #Full path and file name of the comma separated text file that specifies that genes in the targetted gene set

print(paste0("Inputted Parameters: seed=", seed, ' theDir=', theDir, " whichData=", whichData, " changeType=", changeType, 
             " changeAmount=", changeAmount, " changeProp=", changeProp, " genesetType=", geneSetType, " geneset=", 
             geneSet, " bN=", bN, ' seed=', seed, ' changeDirection=', changeDirection, ' resultsDir=', resultsDir))



## Load in the Expression Data
## Note: The expression data is expected as a r datafile
## Note: The expression data is expected to be collected by microarray and already be background corrected and normalized
## Note: Each row in the data matrix is expected to be a sample, the first column is the class label (0=control, 1=case), the other columns are the expression of the specific gene as indicated by the column name
## Note: This specific directory structure will have to be defined by the user, thus it will just be easier to have it user 
##       defined
print("Loading in Expression Data")
load(dataPath)
centered.data <- get(dataObject.name)
print("Loaded in Expression Data!!")


## Differentially express the specific genes contained in the pathway
#geneChange.file <- paste0(theDir, 'geneSet_changes/', geneSetType, '_', geneSet, '/', geneSetType, '_', geneSet, "_genesChange.txt")
#resultsDir <- paste0(theDir, 'geneSet_changes/', geneSetType, '_', geneSet, '/withBootstrap/results/bootstrapData/', changeDirection, '/pi_', changeProp, '/seed_', seed, 
#                     '/tau_', changeAmount, '/')
#print(resultsDir)
## Call the function to change the expression of the specified gene set with the appropriate parameters
changed.centered.data <- changeExprs(centered.data, changeType, geneChange.file, geneSetType, geneSet, changeAmount, changeProp, resultsDir, changeDirection)


## Bootstrap Sampling
print("Conducting the Bootstrap Sampling!!")
for (b in 1:bN) {
  
  ## Sample the Dataset with Replacement
  bootstrap.data <- matrix(0, nrow=nrow(changed.centered.data), ncol=ncol(changed.centered.data))
  for (j in 1:nrow(changed.centered.data)) {
    theSample <- sample(1:nrow(changed.centered.data), 1)
    bootstrap.data[j,] <- changed.centered.data[theSample,]
  }
  colnames(bootstrap.data) <- colnames(changed.centered.data)
  bootstrap.path <- resultsDir
  ## Save the bootstrapped data as a R dataset
  save(bootstrap.data, file=paste0(bootstrap.path, "bootstrapData_", changeAmount, "_", changeProp, '_', b, '.rda'))
  
  ## Specific file formatting to be compatable with the GSEA software
  bootstrap.data <- t(bootstrap.data)
  format.bootstrap <- cbind(rownames(bootstrap.data)[2:nrow(bootstrap.data)], rep(NA, nrow(bootstrap.data)-1), bootstrap.data[2:nrow(bootstrap.data),])
  colnames(format.bootstrap) <- c('NAME', 'DESCRIPTION', seq(1,ncol(format.bootstrap)-2))
  gsea.outDir <- paste0(resultsDir, 'withBootstrap/results/bootstrapData/OverExpress/pi_', changeProp, '/diffSeed/tau_', changeAmount)
  write.table(format.bootstrap, file=paste0(bootstrap.path, 'bootstrapData_', changeAmount, '_', changeProp, '_', b, '_gsea.txt'), row.names=F, quote=F, sep='\t')
  
  classFile <- paste0(bootstrap.path, 'bootstrapData_', changeAmount, '_', changeProp, '_', b, '_gsea_labels.cls')
  sink(classFile)
  cat('424 2 1\n')
  if (c(bootstrap.data[1,])[1]==1) {
    cat('# CANC NORM\n')
  } else {
    cat('# NORM CANC\n')
  }
  for (p in 1:length(c(bootstrap.data[1,]))) {
    cat(paste0(c(bootstrap.data[1,])[p], '\t'))
  }
  cat('\n')
  sink()
  
}
## Done!!
print("Finished Bootstrapping the Data!!")

