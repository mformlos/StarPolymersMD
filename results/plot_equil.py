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
plt.rcParams.update({'font.size':36})
plt.rcParams.update({'text.latex.preamble': ['\usepackage{amssymb}', '\usepackage{amsmath}']})
plt.rcParams.update({'xtick.labelsize':24})
plt.rcParams.update({'ytick.labelsize':24})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')



filename = str(sys.argv[1])


Arm, Lambda, N_patch, N_patch_dev, S_patch, S_patch_dev, R_gyr, R_gyr_dev, S, S_dev, delta, delta_dev, c, c_dev = np.loadtxt(filename, unpack=True, usecols=(0,1, 9, 10, 11, 12, 13, 14, 39, 40, 41, 42, 43, 44))
print Arm, Lambda

characteristics = {'N_patch':[N_patch,N_patch_dev, r"$N_p$"],'S_patch':[S_patch, S_patch_dev, r"$S_p$"], 'R_gyr':[R_gyr, R_gyr_dev, r"$R_{gyr}$"], 'S':[S, S_dev, r"$S$"], 'delta':[delta, delta_dev, r"$\delta$"], 'c':[c, c_dev, r"$c$"]}

Arms_datapoints = np.zeros(3)
Arms = np.zeros(3)
Arm_color = [12, 16, 18]
Arm_marker = ['s', 'D', 'o']
current_arm = Arm[0]
Arms[0] = current_arm
i = 0
for value in Arm:
    if value == current_arm:
        Arms_datapoints[i] += 1
    else:
        i += 1
        Arms_datapoints[i] += 1
        current_arm = value
        Arms[i] = current_arm

print Arms, Arms_datapoints

figno = 0
for key, value in characteristics.iteritems():
    fig = plt.figure(figno, figsize=(15,10), dpi =80)
    plt.subplot(1,1,1)
    counter = 0
    for i in range(3):
        Lambda_temp = np.zeros(Arms_datapoints[i])
        data = np.zeros(Arms_datapoints[i])
        data_dev = np.zeros(Arms_datapoints[i])
        for j in range(int(Arms_datapoints[i])):
            Lambda_temp[j] = Lambda[counter]
            data[j] = value[0][counter]
            data_dev[j] = value[1][counter]
            counter += 1
        labelarm = str(Arms[i])
        plt.errorbar(Lambda_temp, data, yerr=data_dev, marker=Arm_marker[i], mew=1.5, markersize=10, lw=1.5, color=tableau20[Arm_color[i]], mec=tableau20[Arm_color[i]], label =r"$f= "+labelarm+"$")
    plt.xlabel(r"$\lambda$")
    plt.ylabel(value[2])
    plt.legend(loc='upper right', numpoints = 1, prop={'size':24}, frameon=False)
    fig.savefig(key+'_equil.pdf')
    plt.show()
    figno += 1






