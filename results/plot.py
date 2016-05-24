import glob
import sys
import numpy as np
import matplotlib.pyplot as plt


###COLOURS###
# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

#setting up the figure
plt.rcParams.update({'font.size':20})
plt.rcParams.update({'text.latex.preamble': ['\usepackage{amssymb}', '\usepackage{amsmath}']})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


filename = str(sys.argv[1])
Arms = filename.split("Arms")[1]
Arms = Arms.split(".dat")[0]

Lambda, Shear, N_patch, S_patch, R_gyr, S, delta, c = np.loadtxt(filename, unpack=True, usecols=(1, 2, 9, 11, 13, 39, 41, 43))

characteristics = {'N_patch':[N_patch, r"$N_p$"], 'S_patch':[S_patch, r"$S_p$"], 'R_gyr':[R_gyr, r"R_{gyr}"], 'S':[S, r"$S$"], 'delta':[delta, r"$\delta$"], 'c':[c, r"$c$"]}

Lamda_datapoints = np.zeros(3)
Lambdas = np.zeros(3)

Lambda_color = [12, 16, 18]
Lambda_marker = ['s', 'D', 'o']

current_lambda = Lambda[0]
Lambdas[0] = current_lambda
i = 0

for value in Lambda:
    if value == current_lambda:
        Lamda_datapoints[i] += 1
    else:
        i += 1
        Lamda_datapoints[i] += 1
        current_lambda = value
        Lambdas[i] = current_lambda


figno = 0
for key, value in characteristics.iteritems():
    fig = plt.figure(figno, figsize(15,10), dpi=80)
    plt.subplot(1,1,1)
    counter = 0
    for i in range(3):
        Shear_temp = np.zeros(Lamda_datapoints[i])
        data = np.zeros(Lamda_datapoints[i])
        for j in range(int(Lamda_datapoints[i])):
            Shear_temp[j] = Shear[counter]
            data[j] = value[0][counter]
            counter += 1
        labellam = str(Lambdas[i])
        plt.plot(Shear_temp, data, marker=Lambda_marker[i], mew=1.5, markersize=10, lw=1.5, color=tableau20[Arm_color[i]], mec=tableau20[Arm_color[i]], label =r"$\lambda = "+labellam+"$")
    plt.xlabel(r"$W_i$")
    plt.ylabel(value[1])
    plt.legend(loc='upper right', numpoints=1, prop={'size':20})
    fig.savefig(key+'_f'+Arms+'.pdf')
    figno += 1




