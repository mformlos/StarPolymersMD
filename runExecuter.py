#!/usr/bin/env python

import subprocess
import time
import params
import os
import shutil
import sys

__author__ = 'maud'

"""Das Skript hat eine Serie von Simulationen z.B. am Cluster auszufuehren:
Es startet soviele Jobs, wie cores verfuegbar stehen. Danach prueft es alle paar Sekunden, ob Cores frei sind
und startet neue jobs, wenn notwendig."""

#By jobInd the script knows which job bundle it needs to execute
jobInd= int(sys.argv[1])
print "Doing job with jobInd= "+ str(jobInd)
jobParamList= params.jobParametersParts[jobInd]

runningList= []
procList= []
outputFileList=[]
jobsTodoList= []

nTask = 0
jobcount = 0
for currJobParam in jobParamList:
    runningList.append(False)
    procList.append(None)
    outputFileList.append(None)
    #jobsTodoList.append(params.nJobPerTask)
    jobsTodoList.append(1)
    nTask= nTask+1
nFreeCores= params.coresPerRun

while True:
    noRun= True
    for ind in range(0,nTask):
        if runningList[ind]:
            noRun= False
            if procList[ind].poll()!= None:
                print "job "+str(ind)+" finished"
                nFreeCores+=1
                runningList[ind]= False

           
    noNewWork= True

    #By using myRange as the range to iterate through the indices, we make sure that we
    #rather do jobs for tasks, where there are still a lot of jobs to do:
    sp= sorted(zip(range(0,nTask),jobsTodoList), key=lambda sp: sp[1], reverse= True)
    myRange= zip(*sp)[0]
    for ind in myRange:
        if not runningList[ind] and jobsTodoList[ind]>0 and nFreeCores>=1:
            noNewWork= False

            currExec = "./star-polymers/build/src/star-polymers %s %s %s %s %s %s %s %s %s %s %s %s %s %s" %(jobParamList[ind].TypeA, jobParamList[ind].TypeB, jobParamList[ind].Arms, jobParamList[ind].Lambda, jobParamList[ind].Temperature, jobParamList[ind].Lx, jobParamList[ind].Ly, jobParamList[ind].Lz, jobParamList[ind].step_size, jobParamList[ind].step_warm, jobParamList[ind].step_total, jobParamList[ind].step_output, jobParamList[ind].MPC, jobParamList[ind].Shear)

            print (currExec)
            procList[ind]= subprocess.Popen(currExec, stdout =open("output.txt","w"), shell= True)
	    jobcount+= 1
            runningList[ind]= True
            jobsTodoList[ind]-=1
            nFreeCores-=1
    if(noRun and noNewWork):
        print "breaking"
        break
    #The following command makes it possible to read the stdout during the execution of the script
    sys.stdout.flush()
    time.sleep(2)

