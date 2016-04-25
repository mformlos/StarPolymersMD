import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size':24})
fig = plt.figure(figsize=(10,8), dpi =80)
plt.subplot(1,1,1)

N = str(sys.argv[1])
filename = "fluid_N"+N+".dat"


y, vx = np.loadtxt(filename, unpack=True)
plt.xlabel(r"$y$")
plt.ylabel(r"$v_x$")

plt.plot(y, vx)

fig.savefig("velocity_profile_N"+N+".pdf")
plt.show()



