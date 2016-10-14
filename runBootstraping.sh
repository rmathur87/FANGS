#!/bin/sh

#Rscript bootstrapSampling.R whichData changeType changeAmount changeProp geneSetType geneSet bN changeDirection dataPath dataObject.name resultDir geneChange.file
Rscript bootstrapSampling.R "prostate" "geneSet" $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} 
