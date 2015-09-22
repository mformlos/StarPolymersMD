import glob
import numpy as np
import sys 

f_output = open("./results/averages.dat", "w") 
filenames = "./results/statistics*"
equil = 0
if len(sys.argv) > 1: 
    equil = float(sys.argv[1])

for file in glob.glob(filenames): 
    averages = np.zeros(15)
    gyration = np.zeros((3,3))
    N = 0
    stat = open(file, "r")
    print file
    for i, line in enumerate(stat): 
        data = line.split(" ")
        if float(data[0])>equil:   
            for j in range(15): 
                averages[j] += float(data[j+1])

            N += 1
    stat.close()

    for j in range(15): 
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

    f_output.write("%s" %file) 
    for j in range(15):
        f_output.write(" %f" %averages[j])
    for j in range(3): 
        f_output.write(" %f" %eigenvalues[j])
    f_output.write("\n")

    
    
    

