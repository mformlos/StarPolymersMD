#computes the averages of all currently available files in the directory given as first input parameter. If an input number is given as second parameter, it will only start calculating averages after this timestep: 

import glob
import numpy as np
import sys 

def sort_files(file_name):
    f = file_name.split('Arms')
    f = f[1].split('_')
    f = f[0]
    f = float(f)
    lam = file_name.split('Lambda')
    lam = lam[1].split('_')
    lam = lam[0]
    lam = float(lam)
    s = file_name.split('Shear')
    s = s[1].split('.dat')
    s = s[0]
    s = float(s)
    print s

    return f+lam+s



directory = str(sys.argv[1])
directory = directory.replace('/', '')
f_output = open("./averages_" + directory+".dat", "w") 
f_output.write("#1Arms  2Lambda  3Shear  4Epot  6Ekin  8Temp  10N_patch  12S_patch 14R_gyr  16G_xx  18G_xy  20G_xz  22G_yx  24G_yy  26G_yz  28G_zx  30Gzy  32Gzz  34wx  36wy  38wz 40S 42delta 44c  46m_g 48ev1 50ev2 52ev3 \n")
filenames = "./"+str(sys.argv[1])+"statistics*"
equil = 100000000

print filenames

if len(sys.argv) > 2: 
    equil = float(sys.argv[2])

for file in sorted(glob.glob(filenames), key=sort_files):
    averages = np.zeros(22)
    stdev = np.zeros(22)
    gyration = np.zeros((3,3))
    N = 0
    N_patch = 0
    S = 0
    stat = open(file, "r")
    print file
    for i, line in enumerate(stat): 
        if i == 0:
            continue
        data = line.split(" ")
        if float(data[0])>equil:   
            for j in range(18):
                if j == 3 and float(data[j+1]) > 0.0: 
                    N_patch += 1   
                if float(data[j+1])>1000: 
                    print i, line 
                averages[j] += float(data[j+1])

                   
                #stdev[j] += np.power(float(data[j+1]),2)
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
            I_a = ev[0]+ev[1]+ev[2]
            I_b = ev[0]*ev[1]+ev[1]*ev[2]+ev[0]*ev[2]
            S = (3.*ev[0]-I_a)*(3.*ev[1]-I_a)*(3.*ev[2]-I_a)/pow(I_a,3)
            averages[18] += S  #prolateness S
            averages[19] += I_b/np.power(I_a,2) #anisotropy delta
            averages[20] += (ev[1]-ev[2])/I_a #acylindricity c
            averages[21] += 2.*gyration[0,1]/(gyration[1,1]-gyration[0,0]) #orientational resistance
            #stdev[18] += np.power(temp, 2)

    for j in range(22):
        if j == 4: 
            #stdev[j] = np.sqrt((1/(N_patch-1))*(stdev[j] - (1/N_patch)*np.power(averages[j],2)))
            averages[j] /= float(N_patch)


        else:
            #print stdev[j], np.power(averages[j], 2)/N 
            #stdev[j] = np.sqrt((1/(N-1))*(stdev[j] - (1/N)*np.power(averages[j],2)))
            averages[j] /= float(N)
            #print stdev[j]


    stat.seek(0)
    for i, line in enumerate(stat): 
        if i == 0:
            continue
        data = line.split(" ")
        if float(data[0])>equil:   
            for j in range(18):
                if (j < 4 and j > 5) or float(data[j+1]) > 0.0:  
                    stdev[j] += np.power(float(data[j+1]) - averages[j],2)
  
    for j in range(19):
        if j == 4:
            stdev[j] = np.sqrt(stdev[j]/((N_patch-1)))

        else:
            #print stdev[j]
            stdev[j] = np.sqrt(stdev[j]/((N-1)))
            #print stdev[j]               
    
    averages[19] = 1 - 3*averages[19]
    stdev[19] *= 3
    #averages[21] *= 27856
    stat.close()

   

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
    Arms = Arms.replace('/statistics','')
    Lambda = file.split('_Lambda')
    Lambda = Lambda[1].split('_')
    Lambda = Lambda[0]
    Shear = file.split('_Shear')
    if len(Shear) > 1:
        Shear = Shear[1].split('.dat')
        Shear = Shear[0]
    else:
        Shear = "MPCOFF"


    f_output.write("%s %s %s" %(Arms, Lambda, Shear)) 
    for j in range(22):
        f_output.write(" %f %f" %(averages[j], stdev[j]))
    for j in range(3): 
        f_output.write(" %f" %eigenvalues[j])
    f_output.write("\n")

    
    
    

