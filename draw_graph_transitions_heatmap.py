"""
Usage: draw_graph_vary_L.py <dir_name> <state>

Arguments
    dir_name    : name of the directory from which to read in data files
    state       : state we want to draw a graph of, e.g. (0, 0)
Options
    -h          : displays this help file
"""
from __future__ import division
import os
from csv import writer
import yaml
import docopt
import matplotlib.pyplot as plt
import numpy as np

arguments = docopt.docopt(__doc__)
dirname = arguments['<dir_name>']
state = arguments['<state>']

r11s = [i/40.0 for i in range(41)]

mean_times_to_deadlock = []
for r11 in r11s:
    r12s = [i/40.0 if (i/40.0)+r11<=1.0 else 'None' for i in range(41)]
    root = os.getcwd()
    directory = os.path.join(root, dirname)
    a_list = []
    for r12 in r12s:
        if r12 != 'None':
            parameter_file_name = directory + str(r11) + '/theoretical_results_%s.yml' % str(r12)
            parameter_file = open(parameter_file_name, 'r')
            parameters = yaml.load(parameter_file)
            parameter_file.close()
            a_list.append(parameters[state])
        else:
            a_list.append(np.nan)
    mean_times_to_deadlock.append(a_list)


x = r11s
y = r11s
X,Y = np.meshgrid(x, y)

data = np.array(mean_times_to_deadlock)
masked_data = np.ma.masked_where(np.isnan(data),data)

plt.pcolor(X, Y, masked_data, cmap='PuBuGn_r')
plt.colorbar()
plt.xlabel(r'$r_{12}$')
plt.ylabel(r'$r_{11}$')
plt.title('Mean Times to Deadlock from (0, 0) - Varying + ' r'$r_{11}$' + ' and ' + r'$r_{12}$')
plt.show()