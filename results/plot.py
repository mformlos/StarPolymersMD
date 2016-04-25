import glob
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size':24})
fig = plt.figure(figsize=(10,8), dpi =80)
plt.subplot(1,1,1)

filename = str(sys.argv[1])
Arms = filename.split("Arms")[1]
Arms = Arms.split(".dat")[0]

Lambda, Shear, E_pot, N_patch, S_patch, D_patch, R_gyr = np.loadtxt(filename, unpack=True, usecols=(2, 3, 4, 10, 12, 14, 18))

Lamda_datapoints = np.zeros(3)
Lambdas = np.zeros(3)
current_lambda = Lambda[0]
Lambdas[0] = current_lambda
i = 0
for value in Lambda:
    if value == current_lambda:
        Lamda_datapoints[i] += 1
    else:
        i += 1
        Lamda_datapoints[i] += 1
        current_lambda == value
        Lambdas[i] = current_lambda

counter = 0
for i in range(3):
    Shear_temp = np.zeros(Lamda_datapoints[i])
    E_pot_temp = np.zeros(Lamda_datapoints[i])
    for j in range(Lamda_datapoints[i]):
        Shear_temp[j] = Shear[counter]
        E_pot_temp[j] = E_pot[counter]
        counter += 1
    labellam = str(Lambdas[i])
    plt.plot(Shear_temp, E_pot_temp, label = labellam)






