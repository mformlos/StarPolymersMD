

import subprocess
import time
import params
import os
import shutil
import sys
import signal


def signal_handler(signum, frame): 
    print("exiting with sigint")
    for i in range(len(jobParamList)):
       print("killing child %i"%i)
       os.killpg(procList[i].pid, signal.SIGINT)
       #time.sleep(5)
    time.sleep(30)
    os.system("python submit.py")
    #subprocess.Popen("python submit.py", stdout=open("subout.txt","w"), shell=True, preexec_fn=os.setsid)
    time.sleep(60)
    sys.exit(0)

signal.signal(signal.SIGUSR1, signal_handler) 


__author__ = 'maud'

"""Das Skript hat eine Serie von Simulationen z.B. am Cluster auszufuehren:
Es startet soviele Jobs, wie cores verfuegbar stehen. Danach prueft es alle paar Sekunden, ob Cores frei sind
und startet neue jobs, wenn notwendig."""

#By jobInd the script knows which job bundle it needs to execute
jobInd= int(sys.argv[1])
print "Doing job with jobInd= "+ str(jobInd)
jobParamList= params.jobParametersParts[jobInd]
print "number of jobs: %i" %len(jobParamList)
runningList= []
procList= []
outputFileList=[]
jobsTodoList= []

print os.getpid()

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
            if (type(jobParamList[ind]) == params.ParamSet): 
                currExec = "./star-polymers/build/src/star-polymers %s %s %s %s %s %s %s %s %s %s %s %s %s %s" %(jobParamList[ind].TypeA, jobParamList[ind].TypeB, jobParamList[ind].Arms, jobParamList[ind].Lambda, jobParamList[ind].Temperature, jobParamList[ind].Lx, jobParamList[ind].Ly, jobParamList[ind].Lz, jobParamList[ind].step_size, jobParamList[ind].step_warm, jobParamList[ind].step_total, jobParamList[ind].step_output, jobParamList[ind].MPC, jobParamList[ind].Shear)
            elif (type(jobParamList[ind]) == params.ParamSetContinue): 
                currExec = "./star-polymers/build/src/star-polymers %s %s %s %s %s %s %s %s %s"%(jobParamList[ind].File, jobParamList[ind].step_size, jobParamList[ind].step_total, jobParamList[ind].step_output, jobParamList[ind].Lx, jobParamList[ind].Ly, jobParamList[ind].Lz, jobParamList[ind].MPC, jobParamList[ind].Shear)
            print (currExec)
            print ind
            procList[ind]= subprocess.Popen(currExec, shell=True, preexec_fn=os.setsid)
            #procList[ind]= subprocess.Popen(currExec, shell=True, stdout =open("output_%i"%ind+".txt","w"), preexec_fn=os.setsid)
            #procList[ind]=subprocess.Popen(currExec, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	    #out[ind], err[ind] = procList[ind].communicate()
	    jobcount+= 1
            runningList[ind]= True
            jobsTodoList[ind]-=1
            nFreeCores-=1
    sys.stdout.flush()
    for p in procList: 
        p.wait()      
    if(noRun and noNewWork):
        print "breaking"
        break
    #The following command makes it possible to read the stdout during the execution of the script
    time.sleep(2)

