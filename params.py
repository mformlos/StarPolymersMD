# -*- coding: utf-8 -*-
__author__ = 'maud'
import numpy as np
import glob
import os
import sys
from collections import namedtuple
#from src.jobParamLib import JobParam

#if we want to change the g File after a few jobs it's not sufficient to change the name of the gFile here, as this would only change the g File used for
#completely newly started simulations. The name of the gFile is part of the jobParams.in File in the jobFolder. Therefore we need to rename gL6.npz to
#g.npz if we want to start using it in the middle of the run.
def gPath(iTask):
    gDir= "gDir/kU2.5,5.0NR128NCosAlpha116NCosAlpha216NPhi16"
    if iTask==0:
        return "%s/g.npz"%(gDir)
    if iTask==1:
        return "%s/gL4.npz"%(gDir)
    if iTask==2:
        return "%s/gL6.npz"%(gDir)
    assert False,"iTask=%d is not valid"%iTask

gPaths=[]
for iTask in xrange(3):
    gPaths.append(gPath(iTask))

ParamSetMixed=namedtuple("ParamSet","Type, Arms, Lambda, Temperature, Lx, Ly, Lz, step_size, step_warm, step_total, step_output, MPC")
ParamSet=namedtuple("ParamSet","TypeA, TypeB, Arms, Lambda, Temperature, Lx, Ly, Lz, step_size, step_warm, step_total, step_output, MPC, Shear")
Type = namedtuple("Type", ["TypeA", "TypeB"])
MPC = namedtuple("MPC", ["Status", "Shear"])
ParamSetContinue = namedtuple("ParamSetContinue", "File, step_size, step_total, step_output, Lx, Ly, Lz, MPC, Shear")

#paramSets=[ParamSetMixed([Type(3,3), Type(5,5)],[3],[1.1],1.0, 50, 50, 50, 0.01, 1E3, 1E4, 1E3,MPC=[MPC("No", [0.0]), MPC("MPC", [0.0,0.5])])]
paramSets=[ParamSetMixed([Type(3,3), Type(5,5)],[3],[1.0, 1.1],1.0, 50, 50, 50, 0.01, 1E3, 1E7, 1E3,MPC=[MPC("No", [0.0])])]

nJobPerTask=1
coresPerRun=4
tasksPerRun=coresPerRun*nJobPerTask


pythonCmd="python"
stateOutName= "state.out"
stateInName= "state.in"
jobParamInName="jobParam.in"
jobParamOutName="jobParam.out"
simsDirPath="sims"
simDirName= "sim"
statesOutName="states.out"
jobNumFile= "jobNum.txt"

resultsDirName="results"



jobParametersTotal=[]
jobParametersParts=[]
currJobParametersPart=[]


for paramSet in paramSets:
    Temperature = paramSet.Temperature
    Lx = paramSet.Lx
    Ly = paramSet.Ly
    Lz = paramSet.Lz
    step_size = paramSet.step_size
    step_warm = paramSet.step_warm
    step_total = paramSet.step_total
    step_output = paramSet.step_output

    for Types in paramSet.Type: 
        for Arms in paramSet.Arms:
            for Lambda in paramSet.Lambda:
                for Hydrodynamic in paramSet.MPC: 
	            for Shear in Hydrodynamic.Shear: 
		        TypeA = Types.TypeA
		        TypeB = Types.TypeB
                        jobParam= ParamSet(TypeA, TypeB, Arms, Lambda, Temperature, Lx, Ly, Lz, step_size, step_warm, step_total, step_output, Hydrodynamic.Status, Shear)
                        mpc_string = ""
                        if (jobParam.MPC == "MPC"): 
                            mpc_string = "MPCON_Shear%.2f.pdb" %jobParam.Shear
                        else: 
                            mpc_string = "MPCOFF.pdb"
                        file_name = "./results/end_config_Star_A%i_B%i_Arms%i_Lx%i_Ly%i_Lz%i_Lambda%.2f_T%.2f_run*_" %(jobParam.TypeA, jobParam.TypeB, jobParam.Arms, jobParam.Lx, jobParam.Ly, jobParam.Lz, jobParam.Lambda, jobParam.Temperature)
                        file_name += mpc_string 
			files = glob.glob(file_name)
                        files.sort()
         	
                        if (len(files) >= 1):
                            continueJobParam = ParamSetContinue(files[0], jobParam.step_size, jobParam.step_total, jobParam.step_output, jobParam.Lx, jobParam.Ly, jobParam.Lz, jobParam.MPC, jobParam.Shear)
                            run = files[0].split('run')
                            run = run[1].split('_')
                            run = run[0]
                            if (float(jobParam.step_total) > float(run)): 
                                jobParametersTotal.append(continueJobParam)
                                currJobParametersPart.append(continueJobParam)
                        else: 
                            jobParametersTotal.append(jobParam)
                            currJobParametersPart.append(jobParam)
	 	 
                        if len(currJobParametersPart)==tasksPerRun:
                            jobParametersParts.append(currJobParametersPart)
                            currJobParametersPart=[]
print(len(jobParametersTotal))
if (not len(currJobParametersPart)== 0):
    jobParametersParts.append(currJobParametersPart)




