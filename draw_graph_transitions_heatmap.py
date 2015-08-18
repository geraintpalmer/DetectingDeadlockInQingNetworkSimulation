"""
Usage: draw_graph_vary_L.py <dir_name> <state> <statesim> <anorsim>

Arguments
    dir_name    : name of the directory from which to read in data files
    state       : state we want to draw a graph of, e.g. (0, 0)
    statesim    : state we want to draw a graph of, e.g. ((0, 0), (0, 0))
    anorsim     : an if analytical, sim if simulation, err if error
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

arguments = docopt.docopt(__doc__)
dirname = arguments['<dir_name>']
state = arguments['<state>']
statesim = arguments['<statesim>']
anorsim = arguments['<anorsim>']

r11s = [i/40.0 for i in range(41)]

if anorsim == 'an' or anorsim == 'err':
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

if anorsim == 'sim' or anorsim == 'err':
    mean_times_to_deadlock_sim = []
    for r11 in r11s:
        r12s = [i/40.0 if (i/40.0)+r11<=1.0 else 'None' for i in range(41)]
        root = os.getcwd()
        directory = os.path.join(root, dirname)
        a_list = []
        for r12 in r12s:
            if r12 != 'None':
                parameter_file_name = directory + str(r11) + '/deadlocking_times_%s.csv' % str(r12)
                parameter_file = open(parameter_file_name, 'r')
                rdr = reader(parameter_file)
                for row in rdr:
                    if row[0] == statesim:
                        a_list.append(sum([float(obs) for obs in row[1:]])/len([float(obs) for obs in row[1:]]))
                        break
                parameter_file.close()
            else:
                a_list.append(np.nan)
        mean_times_to_deadlock_sim.append(a_list)

errors = []
if anorsim == 'err':
    for i in range(len(mean_times_to_deadlock)):
        errors.append([])
        for j in range(len(mean_times_to_deadlock[i])):
            errors[i].append((mean_times_to_deadlock_sim[i][j] - mean_times_to_deadlock[i][j])/mean_times_to_deadlock[i][j])

x = r11s
y = r11s
X,Y = np.meshgrid(x, y)

if anorsim == 'an':
    data = np.array(mean_times_to_deadlock)
if anorsim == 'sim':
    data = np.array(mean_times_to_deadlock_sim)
if anorsim == 'err':
    data = np.array(errors)
masked_data = np.ma.masked_where(np.isnan(data),data)

# plt.pcolor(X, Y, masked_data, cmap='YlOrBr')
plt.pcolor(X, Y, masked_data, cmap='hot_r')
# plt.pcolor(X, Y, masked_data, cmap='BrBG')
# plt.pcolor(X, Y, masked_data)
plt.colorbar()
plt.xlabel(r'$r_{12}$')
plt.ylabel(r'$r_{11}$')
if anorsim == 'err':
    plt.title('Error in Times to Deadlock')
else:
    plt.title('Mean Times to Deadlock from (0, 0) - Varying ' + r'$r_{11}$' + ' and ' + r'$r_{12}$')
plt.show()