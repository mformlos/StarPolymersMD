# -*- coding: utf-8 -*-
__author__ = 'maud'
import numpy as np
import glob
import os
import sys
import shutil
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

#ParamSetMixed=namedtuple("ParamSet","Type, Arms, Lambda, Temperature, Lx, Ly, Lz, step_size, step_warm, step_total, step_output, MPC, pdb, pdb_out, fluid, fluid_out")
ParamSetMixed=namedtuple("ParamSet","Type, Arms, Mass, Lambda, Temperature, Size, step_size, step_warm, step_total, step_output, MPC, pdb, pdb_out, fluid, fluid_out")

ParamSet=namedtuple("ParamSet","TypeA, TypeB, Arms, Mass, Lambda, Temperature, Lx, Ly, Lz, step_size, step_warm, step_total, step_output, MPC, Shear, pdb, pdb_out, fluid, fluid_out")
Type = namedtuple("Type", ["TypeA", "TypeB"])
MPC = namedtuple("MPC", ["Status", "Shear"])
Size = namedtuple("Size", ["Lx", "Ly", "Lz"]) 
ParamSetContinue = namedtuple("ParamSetContinue", "File, step_size, step_total, step_output, Lx, Ly, Lz, MPC, Shear, pdb, pdb_out, fluid, fluid_out")

#paramSets=[ParamSetMixed([Type(3,3)],[3],10,[1.1],0.5, 50, 50, 50, 0.01, 1E3, 1E4, 1E3, [MPC("MPC", [0.3])], "pdb", 1E6, "No",1E6)]
#paramSets=[ParamSetMixed([Type(21,9)],[6,9,15],10,[0.5,0.8,1.05,1.15],0.5,[Size(120, 120, 120)], 0.001, 1E6, 5E8, 1E4,[MPC("No", [0.0])],"pdb", 1E5, "No",1E6)] 
paramSets=[ParamSetMixed([Type(21,9)],[6,9,15],10,[0.5, 0.8, 1.05, 1.15],0.5, [Size(50, 50, 50)], 0.001, 1E5, 5E8, 1E4, [MPC("MPC", [0.0])],"pdb", 1E5, "No",1E6)] 
#paramSets=[ParamSetMixed([Type(21,9)],[6],10,[0.5, 0.8, 1.05, 1.15],0.5, [Size(50, 50, 50)], 0.001, 1E5, 5E8, 1E3, [MPC("MPC", [0.0001, 0.0002, 0.0005])],"pdb", 1E5, "No",1E6)]
#paramSets=[ParamSetMixed([Type(21,9)],[6],10,[0.5, 0.8, 1.05, 1.15],0.5, [Size(80, 50, 50)], 0.001, 1E5, 5E8, 1E3, [MPC("MPC", [0.001, 0.002, 0.005])],"pdb", 1E5, "No",1E6)]

#paramSets=[ParamSetMixed([Type(21,9)],[15],10,[1.05],0.5, [Size(60,60,60), Size(70,70,70),Size(80,80,80),Size(90,90,90),Size(100,100,100)], 0.001, 1E5, 5E8, 1E3, [MPC("MPC", [0.005])],"pdb", 1E5, "No",1E6)] 


#paramSets=[ParamSetMixed([Type(2,6)],[3],10,[1.1],1.0, 50, 50, 60, 0.01, 1E3, 1E6, 1E3,[MPC("No", [0.0])],"pdb", 1E3, "No",1E3)] 

coresPerJob=16
coresPerRun=16
tasksPerRun=16



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
    #Lx = paramSet.Lx
    #Ly = paramSet.Ly
    #Lz = paramSet.Lz
    step_size = paramSet.step_size
    step_warm = paramSet.step_warm
    step_total = paramSet.step_total
    step_output = paramSet.step_output
    pdb = paramSet.pdb
    pdb_out = paramSet.pdb_out
    fluid = paramSet.fluid
    fluid_out = paramSet.fluid_out
    Mass = paramSet.Mass
    for Types in paramSet.Type: 
        for Arms in paramSet.Arms:
            for Lambda in paramSet.Lambda:
                for Hydrodynamic in paramSet.MPC: 
	            for Shear in Hydrodynamic.Shear: 
                        for Sizes in paramSet.Size:
                            Lx = Sizes.Lx
                            Ly = Sizes.Ly
                            Lz = Sizes.Lz  
		            TypeA = Types.TypeA
		            TypeB = Types.TypeB
                            jobParam= ParamSet(TypeA, TypeB, Arms, Mass, Lambda, Temperature, Lx, Ly, Lz, step_size, step_warm, step_total, step_output, Hydrodynamic.Status, Shear, pdb, pdb_out, fluid, fluid_out)
                            mpc_string = ""
                            if (jobParam.MPC == "MPC"): 
                                mpc_string = "MPCON_Shear%.4f.pdb" %jobParam.Shear
                            else: 
                                mpc_string = "MPCOFF.pdb"
                            file_name = "./results/end_config_Star_A%i_B%i_Arms%i_Mass%i_Lx%i_Ly%i_Lz%i_Lambda%.2f_T%.2f_run*_" %(jobParam.TypeA, jobParam.TypeB, jobParam.Arms, jobParam.Mass, jobParam.Lx, jobParam.Ly, jobParam.Lz, jobParam.Lambda, jobParam.Temperature)
                            file_name += mpc_string 
			    files = glob.glob(file_name)
                            overwrite = True
                            if len(files) < 1 : 
                                overwrite = False
                                files.sort()
                        #file_name = file_name.replace("Lx%i_Ly%i_Lz%i" %(jobParam.Lx, jobParam.Ly, jobParam.Lz), "Lx*_Ly*_Lz*")
                        #files += glob.glob(file_name)
                            if len(files) < 1 and jobParam.MPC == "MPC": 
                                file_to_copy = "./results/end_config_Star_A%i_B%i_Arms%i_Mass%i_Lx*_Ly*_Lz*_Lambda0.50_T%.2f_run*_MPCOFF.pdb" %(jobParam.TypeA, jobParam.TypeB, jobParam.Arms, jobParam.Mass, jobParam.Temperature) #file_name.replace(mpc_string, "MPCOFF*.pdb")
                                files_to_copy = glob.glob(file_to_copy)
                                file_to_generate = file_name.replace("*", "0") 
                                if len(files_to_copy) >= 1: 
                                    shutil.copyfile(files_to_copy[0], file_to_generate)
                                    print "generating start config file ", file_to_generate    
                                files += glob.glob(file_name) 
                            
         	
                            if (len(files) >= 1):
                                continueJobParam = ParamSetContinue(files[0], jobParam.step_size, jobParam.step_total, jobParam.step_output, jobParam.Lx, jobParam.Ly, jobParam.Lz, jobParam.MPC, jobParam.Shear, jobParam.pdb, jobParam.pdb_out, jobParam.fluid, jobParam.fluid_out)
                                run = files[0].split('run')
                                run = run[1].split('_')
                                run = run[0]
                            
                                if (jobParam.MPC == "MPC" and files[0].find("MPCON")==-1):
                                    overwrite = False
                                if (float(jobParam.step_total) > float(run) or overwrite == False): 
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




