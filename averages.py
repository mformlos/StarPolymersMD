#computes the averages of all currently available files in the results/ directory. If an input number is given, it will only start calculating averages after this timestep: 

import glob
import numpy as np
import sys 

f_output = open("./results/averages.dat", "w") 
filenames = "./results/statistics*"
equil = 0
if len(sys.argv) > 1: 
    equil = float(sys.argv[1])

for file in sorted(glob.glob(filenames)): 
    averages = np.zeros(16)
    gyration = np.zeros((3,3))
    N = 0
    N_patch = 0
    S = 0
    stat = open(file, "r")
    print file
    for i, line in enumerate(stat): 
        data = line.split(" ")
        if float(data[0])>equil:   
            for j in range(15): 
                if j == 4 and float(data[j+1]) > 0.0: 
                    N_patch += 1   
                averages[j] += float(data[j+1])
            N += 1
            gyration[0,0] = float(data[7])
            gyration[0,1] = float(data[8])  
            gyration[0,2] = float(data[9])
            gyration[1,0] = gyration[0,1]
            gyration[1,1] = float(data[11])
            gyration[1,2] = float(data[12])
            gyration[2,0] = gyration[0,2]
            gyration[2,1] = gyration[1,2]
            gyration[2,2] = float(data[15])
            eigenvalues, eigenvectors = np.linalg.eig(gyration)
            eigenvalues = np.sort(eigenvalues) 
            ev = eigenvalues[::-1]
            I = ev[0]+ev[1]+ev[2]
            averages[15] += (3.*ev[0]-I)*(3.*ev[1]-I)*(3.*ev[2]-I)/pow(I,3)

    stat.close()

    for j in range(16):
        if j == 5: 
            averages[j] /= float(N_patch) 
        else:  
            averages[j] /= float(N)

    gyration[0,0] = averages[6]
    gyration[0,1] = averages[7]  
    gyration[0,2] = averages[8]
    gyration[1,0] = gyration[0,1]
    gyration[1,1] = averages[10]
    gyration[1,2] = averages[11]
    gyration[2,0] = gyration[0,2]
    gyration[2,1] = gyration[1,2]
    gyration[2,2] = averages[14]

    eigenvalues, eigenvectors = np.linalg.eig(gyration)
    A = file.split('_A')
    A = A[1].split('_')
    A = float(A[0])
    B = file.split('_B')
    B = B[1].split('_')
    B = float(B[0])
    alpha = str(B/(A+B))
    Arms = file.split('Arms')
    Arms = Arms[1].split('_')
    Arms = Arms[0]
    Lambda = file.split('_Lambda')
    Lambda = Lambda[1].split('_')
    Lambda = Lambda[0]

    f_output.write("%s %s %s" %(Arms, alpha, Lambda)) 
    for j in range(16):
        f_output.write(" %f" %averages[j])
    for j in range(3): 
        f_output.write(" %f" %eigenvalues[j])
    f_output.write("\n")

    
    
    

