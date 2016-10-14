##### Contains the function which changes/differentially express select genes within a gene set.
##### Version: 1.0 - Sept. 2016
##### Author: Ravi Mathur (rmathur2@ncsu.edu), Motsinger-Reif Lab, Bioinformatics Research Center, NCSU



library(propagate)

## Read the parameters file for the particular seed
fileArgs <- scan("bootstrap_parameters.csv", sep=',', what="")
seed <- as.integer(fileArgs[which(fileArgs == 'seed')+2])
print(paste0("Seed in changeExprs.R = ", seed))
set.seed(seed)
#set.seed('92187')


#### Function calculates the correlation matrix and saves it for a very large data matrix
#### @param dataMatrix - the data matrix to be utilized
#### @param fileName - the path and name of the file that the correlation will be written to
calcCorrelation <- function(dataMatrix, fileName) {
  dataCorrelation <- bigcor(dataMatrix, fun='cor')
  ffsave(dataCorrelation, file=fileName, safe=T)
}



#### Function to change/differentially express select genes within a gene set
#### @param centered.data - the expression matrix - with each row being a sample (n), the first column being the class label, and the other rows each gene expression value, thus m+1 columns. The column names must be the name of the gene (as it appears in MSigDB) or it will not be changed correctly.
#### @param changeType - "geneSet" is only tested - could generalize this to change specific gene(s)
#### @param geneChange.file - path and name of a comma separated text file where each row is the name of a gene that is to be changed
#### @param geneSetType - the MSigDB category type of the gene set being changed
#### @param geneSet - the name of the MSigDB gene set that is being changed
#### @param changeAmount - the amount of signal or differential expression to add for each gene
#### @param changeProp - the proportion of the given gene set to change
#### @param resultsDir - the directory where the genes that were changed are written to
#### @param changeDirection - "OverExpress" or "UnderExpress" - whether the cases are overexpressed or underexpressed
changeExprs <- function(centered.data, changeType, geneChange.file, geneSetType, geneSet, changeAmount, changeProp, resultsDir, changeDirection) {
  
  #### Change Expression for Select Genes
  #### Note: 
  if (changeType == 'geneSet') {
    print("Changing all Genes in the Given Geneset!!")
    
    ## Load in file  with list of genes to change and process the data
    genes.change <- read.table(geneChange.file, sep=',', header=F)
    genes.change <- c(genes.change)[[1]]
    
    ## Find the index of the genes that are contained in the target gene set
    ## Note: The gene name must be exactly as appears in MSig DB or it will not be changed
    allMatching <- match(genes.change, colnames(centered.data))
    genes.change <- genes.change[which(!is.na(allMatching))]
    
    ## Shuffle the Labels, but still keep the proportion of cases to controls
    ## Note: 1=case, 0=control
    centered.data[,1] <- sample(centered.data[,1])
    
    ## Change a proportion of genes in the geneSet (if changeProp=1, all genes are changed; if changeProp is small enough, a single genes is changed)
    numChange <- ceiling(changeProp*length(genes.change)) #Number of genes to change depending on the proportion entered
    selectGenes <- sample(genes.change, numChange) #select a random set of genes to differentially express
    write.table(selectGenes, file=paste0(resultsDir, geneSetType, '_', geneSet, '_', changeAmount, '_', changeProp, '_', changeDirection, '_', seed, 
                                         '_changedGenes.txt'), col.names = F, row.names = F, quote=F) #write the genes that were changed to a file
    
    ## Actually change the expression values for the genes selected
    for (i in 1:length(selectGenes)) {
      thisMapping <- which(colnames(centered.data) == selectGenes[i]) #Find all occurances of this gene (if multiple exist)
      for (aProbe in thisMapping) {
      	
	## Note: 0=control, 1=cases
        if (changeDirection == 'OverExpress') {
          #Decrease the expression of the controls
          centered.data[which(centered.data[,1]==0),aProbe] <- centered.data[,aProbe][which(centered.data[,1]==0)] - (changeAmount * sd(centered.data[,aProbe]))
          #Increase the expression of the cases
          centered.data[which(centered.data[,1]==1),aProbe] <- centered.data[,aProbe][which(centered.data[,1]==1)] + (changeAmount * sd(centered.data[,aProbe]))
        
	} else if (changeDirection == 'UnderExpress') {
          #Increase the expression of the controls
          centered.data[which(centered.data[,1]==0),aProbe] <- centered.data[,aProbe][which(centered.data[,1]==0)] + (changeAmount * sd(centered.data[,aProbe]))
          #Decrease the expression of the cases
          centered.data[which(centered.data[,1]==1),aProbe] <- centered.data[,aProbe][which(centered.data[,1]==1)] - (changeAmount * sd(centered.data[,aProbe]))
        
	}
      }
    }
    
    #dataPath <- paste0(resultsDir, geneSetType, '_', geneSet, '_', changeAmount, '_', changeProp, '_', changeDirection, '_', seed, '_changedData.rda')
    #save(centered.data, file=dataPath)
    print("Changed and Saved Dataset")
    
    ## Calculate and write the correlation of this specific change
    correlation.dataPath <- paste0(resultsDir, geneSetType, '_', geneSet, '_', changeAmount, '_', changeProp, '_', changeDirection, '_', seed, '_corrMat')
    #calcCorrelation(centered.data[,-1], correlation.dataPath)
    print("Correlation Matrix Calculated and Saved!!")
    
    return(centered.data)
    
  } else {
    
    ## Change specific genes
    print("Code Still To Be Written!!")
  }
  
}
