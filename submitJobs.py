__author__ = 'adminuser'
import subprocess
import glob
import sys
import os

runTemplScript= "scriptTempl.sh"
currRunScript= "currRun.sh"
queue= "all.q"

for runscript in glob.glob("currRun*.sh"):
    os.remove(runscript)

currJobInd= 0
for param in sorted(glob.glob("params*.dat")): 
    currJobInd = param.split("params")
    currJobInd = currJobInd[1].split(".dat")
    currJobInd = currJobInd[0]
    currRunScript = "currRun"+currJobInd+".sh"
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
    subprocess.Popen(currExec, shell= True).wait()
    #subprocess.Popen("python Executer.py " + str(currJobInd), shell=True).wait()

sys.exit(0)

