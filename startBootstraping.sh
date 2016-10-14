#!/bin/sh

#sbatch -o bootstrap_changeAmount_changeProp_bN_over --mail-type=ALL --mail-user=rmathur2@ncsu.edu runBoostraping.sh changeAmount changeProp bN
theDir="~/networkSimulator/simulateData/"
changeAmount=$1
changeProp=$2
bN=10
#genesetType="c6"
genesetType="h"
#geneset="KRAS.PROSTATE_UP.V1_DN"
geneset="TGF_BETA_SIGNALING"
seed=92187
changeDirection="OverExpress"
dataPath="${theDir}/datasets/prostateCancer/rawData/data/CELfiles/caseControl_centered_final.rda"
dataObjectName="caseControl_centered"
resultDir="./results/"
if [ ! -d "$resultDir" ]; then
    mkdir $resultDir
fi
geneChangeFile="${theDir}/geneSet_changes/${genesetType}_${geneset}/${genesetType}_${geneset}_genesChange.txt"
echo bootstrap_${changeAmount}_${changeProp}_${bN}_${changeDirection}_${geneset}_${seed}

sbatch -o bootstrap_${changeAmount}_${changeProp}_${bN}_${changeDirection}_${geneset}_${seed} -p short --mail-type=ALL --mail-user=rmathur2@ncsu.edu runBootstraping.sh $changeAmount $changeProp $genesetType $geneset $bN $changeDirection $dataPath $dataObjectName $resultDir $geneChangeFile
#./runBootstraping.sh $changeAmount $changeProp $genesetType $geneset $bN $changeDirection
