#!/usr/bin/env python

__author__ = 'adminuser'
import subprocess
import params
import sys

runTemplScript= "runTempl.sh"
currRunScript= "currRun.sh"
queue= "all.q"

currJobInd= 0
for jobParametersPart in params.jobParametersParts:
    currRunScript = "currRun"+"%i"%currJobInd+".sh"
    print "Submitting run" \
          " for the following jobParametersPart:"
    with open(runTemplScript, "r") as inF:
        lines = inF.readlines()
    with open(currRunScript, "w") as outF:
        lineInd = 0
        for currLine in lines:
            if lineInd== 13:
                outF.write("jobInd="+  str(currJobInd)+"\n")
            else:
                outF.write(currLine)
            lineInd += 1

    currExec= "qsub.py -q "+queue+" "+currRunScript
    print (currExec)
    #subprocess.Popen(currExec, shell= True).wait()
    subprocess.Popen("python runExecuter.py " + str(currJobInd), shell=True).wait()
    currJobInd+=1

sys.exit(0)
