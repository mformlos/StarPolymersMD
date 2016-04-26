import sys
import glob
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
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(0, figsize=(10,8), dpi =80)
plt.subplot(1,1,1)

#getting parameters 
profile_filename = str(sys.argv[1]) 
N = profile_filename.split("N")[1]
N = N.split("_")[0]
shear = profile_filename.split("Shear")[1]
shear = shear.split("_")[0]
checkfile = open(profile_filename, "r")
linenumber = 0
for i, line in enumerate(checkfile):
    if line in ['\n', '\r\n']:
        linenumber = i
checkfile.close()  
param_string = profile_filename.split("fluid_")[1]
param_string = param_string.split(".dat")[0]




###SHEAR VELOCITY###
#extract data from file
x, vxx, vxy = np.genfromtxt(profile_filename, usecols=(0,1,2), skip_footer=linenumber, unpack=True)
y, vyx, vyy = np.genfromtxt(profile_filename, usecols=(0,1,2), skip_header=linenumber, unpack=True)

#setting plot limits and labels
plt.xlim(np.round(x[0])-0.5, np.round(x[len(x)-1])+0.5)
plt.xlabel(r"$y$ coordinate")
plt.ylabel(r"average velocity $v_x$ in $x$-direction")

#plot
plt.plot(y, vyx, marker="o", color=tableau20[16], linestyle="None")
plt.plot(y, y*float(shear), color=tableau20[16])

#save and show
fig.savefig("velocity_profile_"+param_string+".pdf")
plt.show()


###MAXWELL DISTRIBUTION###

fig = plt.figure(1, figsize=(10,8), dpi =80)

maxwell_filename = profile_filename.replace("fluid", "maxwell") 

v_abs, prob = np.genfromtxt(maxwell_filename, unpack=True)

plt.xlabel(r"$\vert \mathbf{v} \vert$")
plt.ylabel(r"$P(\vert \mathbf{v} \vert)$")

plt.plot(v_abs, prob, marker="o", markersize=3, color=tableau20[16], linestyle="None")
plt.plot(v_abs, ((1./(2.*np.pi))**(3./2.))*4.*np.pi*(v_abs**2)*np.exp(-v_abs**2/2), color="black")

fig.savefig("maxwell_distribution_"+param_string+".pdf")
plt.show()


###Transverse velocity autocorrelation###

fig = plt.figure(2, figsize=(10,8), dpi=80)
ax = fig.add_subplot(1,1,1)
#ax.set_yscale('log')

transverse_filename = profile_filename.replace("fluid", "fourier")

time, ck1, ck2, ck3 = np.genfromtxt(transverse_filename, unpack=True)
ck1 /= ck1[0]
ck2 /= ck2[0]
ck3 /= ck3[0]

plt.xlabel(r"$t/\sqrt{ma^2/k_BT}$")
plt.ylabel(r"$C_v^T (k,t)$")

plt.semilogy(time, ck1, marker='o', color=tableau20[12], linestyle="None", label=r"$\mathbf{k} = (1,0,0)$" )
plt.semilogy(time, ck2, marker='o', color=tableau20[16], linestyle="None", label=r"$\mathbf{k} = (0,1,0)$")
plt.semilogy(time, ck3, marker='o', color=tableau20[18], linestyle="None", label=r"$\mathbf{k} = (0,0,1)$")
plt.semilogy(time, np.exp(-(2.*np.pi/10.)**2*0.8705*time), color="black")
fig.savefig("transverse_autocorr_"+param_string+".pdf")
plt.show()


###Velocity Autocorrelation in real space###

fig = plt.figure(3, figsize=(10,8), dpi=80)
ax = fig.add_subplot(1,1,1)

autocorr_filename = profile_filename.replace("fluid", "velocity")

time, cv = np.genfromtxt(autocorr_filename, unpack=True)

plt.xlabel(r"$t/\sqrt{ma^2/k_BT}$")
plt.ylabel(r"$C_v(t)$")

plt.loglog(time, cv, marker='o', color=tableau20[16], linestyle="None")

fig.savefig("velocity_autocorr_"+param_string+".pdf")
plt.show()

