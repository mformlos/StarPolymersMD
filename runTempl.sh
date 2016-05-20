#!/bin/bash

#$ -N star_polymers
#$ -P p70679
#$ -pe mpich 16
#$ -V 
#$ -l h_rt=00:10:00
#$ -l s_rt=00:08:00
#$ -M maud.formanek@univie.ac.at
#$ -m beas
#$ -v PATH
#$ -v LD_LIBRARY_PATH

jobInd=

pwd
echo $TMP
echo ${SGE_O_WORKDIR}
echo $NSLOTS

exit_handler() {
    date
    echo "exit caught by shell script" 
    kill -SIGUSR1 "${child}" 2>/dev/null
    echo "child killed"
    exit 255
}

trap exit_handler SIGUSR1

python runExecuter.py $jobInd &

child=$!
echo "child = "${child}""
wait "${child}"
echo finished
