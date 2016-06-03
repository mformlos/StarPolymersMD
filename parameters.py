# -*- coding: utf-8 -*-
__author__ = 'maud'
import numpy as np
import glob
import os
import sys
import shutil
from collections import namedtuple
#from src.jobParamLib import JobParam

for filename in glob.glob("params*.dat"):
    os.remove(filename) 

ParamSetMixed=namedtuple("ParamSet","Type, Arms, Mass, Lambda, Temperature, Size, step_size, step_warm, step_total, step_output, MPC, pdb, pdb_out, fluid, fluid_out")

ParamSet=namedtuple("ParamSet","TypeA, TypeB, Arms, Mass, Lambda, Temperature, Lx, Ly, Lz, step_size, step_warm, step_total, step_output, MPC, Shear, pdb, pdb_out, fluid, fluid_out")
Type = namedtuple("Type", ["TypeA", "TypeB"])
MPC = namedtuple("MPC", ["Status", "Shear"])
Size = namedtuple("Size", ["Lx", "Ly", "Lz"]) 
ParamSetContinue = namedtuple("ParamSetContinue", "File, step_size, step_total, step_output, Lx, Ly, Lz, MPC, Shear, pdb, pdb_out, fluid, fluid_out")

paramSets = []

#paramSets.append([ParamSetMixed([Type(70,30)],[6,9,15],5,[0.8,0.9,1.0,1.1,1.15],0.5, [Size(400, 400, 400)], 0.001, 1E7, 5E8, 1E4,[MPC("No", [0.0])],"No", 1E5, "No",1E6)]) 
paramSets.append([ParamSetMixed([Type(28,12)],[6,9,15],5,[0.5, 0.8, 1.05, 1.15],0.5, [Size(50, 50, 50)], 0.001, 1E7, 5E8, 1E4, [MPC("MPC", [0.0])],"pdb", 1E6, "No",1E6)]) 
paramSets.append([ParamSetMixed([Type(28,12)],[6,9,15],5,[0.5,0.8, 1.05, 1.15],0.5, [Size(50, 50, 50)], 0.001, 1E5, 5E8, 1E4, [MPC("MPC", [0.00001,0.00002, 0.0001, 0.0002])],"pdb", 1E6, "No",1E6)])
paramSets.append([ParamSetMixed([Type(28,12)],[6,9,15],5,[0.5, 0.8, 1.05, 1.15],0.5, [Size(80, 50, 50)], 0.001, 1E5, 5E8, 1E4, [MPC("MPC", [0.0005, 0.001, 0.002, 0.005])],"pdb", 1E6, "No",1E6)])
paramSets.append([ParamSetMixed([Type(28,12)],[6,9,15],5,[0.5, 0.8, 1.05, 1.15],0.5, [Size(100, 50, 50)], 0.001, 1E5, 5E8, 1E4, [MPC("MPC", [0.01,0.02])],"pdb", 1E6, "No",1E6)]) 

#coresPerJob=16
#coresPerRun=16
tasksPerRun=16 #160



pythonCmd="python"

resultsDirName="results"



jobParametersTotal=[]
jobParametersParts=[]
currJobParametersPart=[]


for paramSettings in paramSets:
    for paramSet in paramSettings:
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
                                    mpc_string = "MPCON_Shear%.5f.pdb" %jobParam.Shear
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

submitInd = 0
currExec = " "
for jobParamList in jobParametersParts: 
    submitfilename = "params" + str(submitInd)+".dat" 
    submitfile = open(submitfilename, "w") 
    for ind in range(len(jobParamList)): 
        if (type(jobParamList[ind]) == ParamSet): 
            currExec = "./star-polymers/build/src/star-polymers %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s profile" %(jobParamList[ind].TypeA, jobParamList[ind].TypeB, jobParamList[ind].Arms, jobParamList[ind].Mass, jobParamList[ind].Lambda, jobParamList[ind].Temperature, jobParamList[ind].Lx, jobParamList[ind].Ly, jobParamList[ind].Lz, jobParamList[ind].step_size, jobParamList[ind].step_warm, jobParamList[ind].step_total, jobParamList[ind].step_output, jobParamList[ind].MPC, jobParamList[ind].Shear, jobParamList[ind].pdb, jobParamList[ind].pdb_out, jobParamList[ind].fluid, jobParamList[ind].fluid_out)
        elif (type(jobParamList[ind]) == ParamSetContinue): 
            currExec = "./star-polymers/build/src/star-polymers %s %s %s %s %s %s %s %s %s %s %s %s %s profile"%(jobParamList[ind].File, jobParamList[ind].step_size, jobParamList[ind].step_total, jobParamList[ind].step_output, jobParamList[ind].Lx, jobParamList[ind].Ly, jobParamList[ind].Lz, jobParamList[ind].MPC, jobParamList[ind].Shear, jobParamList[ind].pdb, jobParamList[ind].pdb_out, jobParamList[ind].fluid, jobParamList[ind].fluid_out)
        submitfile.write(currExec + "\n") 
    submitfile.close()
    submitInd += 1
    
