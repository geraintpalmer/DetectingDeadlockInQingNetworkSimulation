"""
Usage: draw_graph_vary_L.py <dir_name> <state> <ran>

Arguments
    dir_name    : name of the directory from which to read in data files
    state       : state we want to draw a graph of, e.g. (0, 0)
    ran         : range of values
Options
    -h          : displays this help file
"""
from __future__ import division
import os
from csv import reader
import yaml
import docopt
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set(style="whitegrid")


arguments = docopt.docopt(__doc__)
dirname = arguments['<dir_name>']
state = arguments['<state>']
ran = int(arguments['<ran>'])

n1s = [i for i in range(ran)]

mean_times_to_deadlock = []
for n1 in n1s:
    n2s = [i for i in range(ran)]
    root = os.getcwd()
    directory = os.path.join(root, dirname)
    a_list = []
    for n2 in n2s:
        parameter_file_name = directory + 'n1_' + str(n1) + '/theoretical_results_%s.yml' % str(n2)
        parameter_file = open(parameter_file_name, 'r')
        parameters = yaml.load(parameter_file)
        parameter_file.close()
        a_list.append(parameters[state])
    mean_times_to_deadlock.append(a_list)


x = n2s
y = n1s
X,Y = np.meshgrid(x, y)


data = np.array(mean_times_to_deadlock)
masked_data = np.ma.masked_where(np.isnan(data),data)

# plt.pcolor(X, Y, masked_data, cmap='YlOrBr', linewidth=1, linestyle='solid', edgecolor='black')
plt.pcolor(X, Y, masked_data, cmap='hot_r', linewidth=1, linestyle='solid'gecolor='black')
# plt.pcolor(X, Y, masked_data, cmap='BrBG')
# plt.pcolor(X, Y, masked_data, linewidth=1, linestyle='solid', edgecolor='black')
plt.colorbar()
plt.xlabel(r'$n_2$')
plt.ylabel(r'$n_1$')
plt.title('Mean Times to Deadlock from (0, 0) - Varying ' + r'$n_1$' + ' and ' + r'$n_2$', fontsize=18)
plt.show()