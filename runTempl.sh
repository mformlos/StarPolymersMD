#!/bin/bash

#$ -N star_polymers
#$ -P p70679
#$ -pe mpich 16
#$ -V 
#$ -M maud.formanek@univie.ac.at
#$ -m beas
#$ -v PATH
#$ -v LD_LIBRARY_PATH

jobInd=

pwd
echo $TMP
echo ${SGE_O_WORKDIR}
echo $NSLOTS

cd runIn

python runExecuter.py $jobInd

echo finished
