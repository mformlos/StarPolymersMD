#calculates the average after each av_step steps; input parameters: file name, av_step

import Gnuplot
import sys 
import numpy as np

if len(sys.argv) < 2 : 
    print "Please provide a filename and at which steps you want to calculate the averages"
    exit()

filename = str(sys.argv[1])
av_step = float(sys.argv[2])
n_output = str.replace(filename, "statistics", "cumulative_average")
f_output = open(n_output, "w")

# e_pot, r_gyr, number of patches, patch size
average = np.zeros(4)
N = 0
f_input = open(filename, "r")
for i, line in enumerate(f_input): 
    data = line.split(" ") 
    average[0] += float(data[1])
    average[1] += float(data[4]) 
    average[2] += float(data[5])
    average[3] += float(data[6])
    N += 1
    if (data[0] >= av_step): 
        f_output.write("%i %f %f %f %f \n" %(int(data[0]), average[0]/N, average[1]/N, average[2]/N, average[3]/N)) 
        av_step += av_step

f_output.close()
g = Gnuplot.Gnuplot(persist = 1)
g('set title "' + n_output + '"') 
g('set ylabel "cumulative average"') 
g('set xlabel "step"') 

title = ["e_pot", "r_gyr", "n_patch", "s_patch"]

for i in range(4): 
    plot_command = 'plot "' + n_output + '" u 1:' + str(i+2) + ' title "' + title[i] + '" '
    g('set term x11 ' + str(i))
    g(plot_command)
 
    n_plot = str.replace(n_output, ".dat", ".png")
    n_plot = n_plot.replace("cumulative_average", "cumulative_average_" + title[i])  
    g('set terminal png size 1024, 768') 
    g('set output "' + n_plot + '"')
    g('repl')
    print "Das Bild '" + n_plot + "' wurde erstellt.\n"




       

