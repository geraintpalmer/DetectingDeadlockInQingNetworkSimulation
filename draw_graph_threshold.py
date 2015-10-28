"""
Usage: draw_graph_threshold.py <dir_name> <var> <legend_loc> <line_colour>

Arguments
    dir_name    : name of the directory from which to read in data files
    legend_loc  : location of the legend

Options
    -h          : displays this help file
"""
from __future__ import division
import os
from csv import writer
import yaml
import shutil
import docopt
import matplotlib.pyplot as plt
from csv import reader
import numpy as np
import pylab
import seaborn as sns
sns.set(style="whitegrid")


arguments = docopt.docopt(__doc__)
dirname = arguments['<dir_name>']
var = str(arguments['<var>'])
leg_loc = str(arguments['<legend_loc>'])
line_colour = str(arguments['<line_colour>'])

root = os.getcwd()
directory = os.path.join(root, dirname)
threshold_results = []

parameter_file_name = directory + 'thresholds_%s.csv' % var
parameter_file = open(parameter_file_name, 'r')
rdr = reader(parameter_file)
for row in rdr:
    threshold_results.append([float(obs) for obs in row])
parameter_file.close()



fig, ax = plt.subplots()

plt.plot(threshold_results[1], threshold_results[0], line_colour, linewidth=3, label='Threshold')

if var[0] == 'L':
    ax.set_xlabel(r'$\Lambda_{'+var[1:]+'}$')
if var[0] == 'm':
    ax.set_xlabel(r'$\mu_{'+var[2:]+'}$')
if var[0] == 'n':
    ax.set_xlabel(r'$n_{'+var[1:]+'}$')
if var[0] == 'c':
    ax.set_xlabel(r'$c_{'+var[1:]+'}$')
if var[0] == 'r':
    ax.set_xlabel(r'$r_{'+var[1:]+'}$')

ax.set_ylabel('Threshold')
ax.set_title('Effect of varying parameters on the ' + r'$\mu$' + ' Threshold.', fontsize=20)

if leg_loc == 'l':
    plt.legend(loc=2, prop={'size':16})
if leg_loc == 'r':
    plt.legend(loc=1, prop={'size':16})
plt.show()