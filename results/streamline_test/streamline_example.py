import numpy as np
from matplotlib import pyplot as plt

Lx = 10
Ly = 10

x = np.linspace(0,10,10)
y = np.linspace(0,10,10)

xm,ym = np.meshgrid(x,y)

u = 0.1*xm
v = -0.1*ym

plt.quiver(xm,ym,u,v)

plt.show()

f_output = open("example_field.dat", "w")

for i in range(len(x)):
    for j in range(len(y)):
        f_output.write("%f %f " %(u[i][j], v[i][j]))
    f_output.write("\n")
f_output.close()