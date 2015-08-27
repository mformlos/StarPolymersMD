#!/bin/bash

#$ -N star_polymers
#$ -P p70679
#$ -pe mpich 16
#$ -V 
#$ -l h_rt=00:05:00
#$ -l s_rt=00:01:00
#$ -M maud.formanek@univie.ac.at
#$ -m beas
#$ -v PATH
#$ -v LD_LIBRARY_PATH

jobInd=

pwd
echo $TMP
echo ${SGE_O_WORKDIR}
echo $NSLOTS

function exit_handler() {
    date 
    echo "exit caught by shell script"
    kill -12 $child 2>/dev/null
    echo "child killed"
    exit 255
}

trap exit_handler SIGUSR2 12

cd runIn

python runExecuter.py $jobInd &

child = $!
wait $child
wait $child
echo finished
