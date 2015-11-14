import sys
import numpy as np 
import matplotlib.pyplot as plt
import linecache as lc

print np.version.version

L=50
filename = sys.argv[1]
d = 3


line = lc.getline(filename, d)
line = line.split(" ")
line.pop()
x = np.array(line, dtype=float)

line = lc.getline(filename, d+1)
line = line.split(" ")
line.pop()
y = np.array(line, dtype=float)

print x
print y

vx = np.zeros((L,L)) 
vy = np.zeros((L,L)) 

for i in range(0,L): 
    line = lc.getline(filename, d+2+i) 
    line = line.split(" ")
    line.pop()
    for j in range(0, L): 
       vx[i][j] = float(line[j*2])
       vy[i][j] = float(line[j*2+1]) 

print vx
print "\n"
print vy


plt.streamplot(x, y, vx, vy)

plt.show()
