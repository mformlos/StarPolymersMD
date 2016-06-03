

import subprocess
import time
import os
import shutil
import sys
import signal


def signal_handler(signum, frame): 
    print("exiting with signal %i" %signum)
    for i in range(len(procList)):
       print("killing child %i"%i)
       if procList[i].poll()==None: 
           os.killpg(procList[i].pid, signal.SIGTERM)
           #time.sleep(5)
    time.sleep(30)
    #os.system("python submit.py")
    #subprocess.Popen("python submit.py", stdout=open("subout.txt","w"), shell=True, preexec_fn=os.setsid)
    #time.sleep(60)
    sys.exit(0)

signal.signal(signal.SIGUSR1, signal_handler) 


__author__ = 'maud'

"""Das Skript hat eine Serie von Simulationen z.B. am Cluster auszufuehren:
Es startet soviele Jobs, wie cores verfuegbar stehen. Danach prueft es alle paar Sekunden, ob Cores frei sind
und startet neue jobs, wenn notwendig."""

#By jobInd the script knows which job bundle it needs to execute
jobInd= int(sys.argv[1])
paramfilename = "params"+str(jobInd)+".dat"
paramfile = open(paramfilename, "r")
print "Doing job with jobInd= "+ str(jobInd)
runningList= []
procList= []
outputFileList=[]
jobsTodoList= []

print os.getpid()

nTask = 0
jobcount = 0

for ind, line in enumerate(paramfile):
    print line
    procList.append(subprocess.Popen(line, shell=True, preexec_fn=os.setsid))
    runningList.append(True)
    nTask += 1 
 

while True: 
    noRun = True
    for ind in range(nTask): 
        if runningList[ind]:
            noRun = False 
            if procList[ind].poll() != None: 
                print "job" + str(ind) + " finished"
                runningList[ind] = False
    if noRun: 
        print "breaking"
        break 
    time.sleep(50)

