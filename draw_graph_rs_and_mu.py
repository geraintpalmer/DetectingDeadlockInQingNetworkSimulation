"""
Usage: draw_graph_threshold.py <dir_name> <state> <legend_loc>

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
state = arguments['<state>']
leg_loc = str(arguments['<legend_loc>'])

root = os.getcwd()
directory = os.path.join(root, dirname)
threshold_results = []

r11s = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
mus = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]

data = []
for r in r11s:
    data.append([])
    for mu in mus:
        parameter_file_name = directory + '%s/theoretical_results_%s.yml' % (str(r), str(mu))
        parameter_file = open(parameter_file_name, 'r')
        parameters = yaml.load(parameter_file)
        parameter_file.close()
        data[-1].append(parameters[state])


fig, ax = plt.subplots()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.78, box.height])

ax.set_xlabel(r'$\mu$')
ax.set_ylabel('Mean Time to Deadlock')
plt.title('Effect of ' + r'$r_{11}$' + ' and ' r'$\mu$' + ' on Time to Deadlock', fontsize=18)


colormap = plt.cm.gray
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(data))])



for i in range(len(data)):
    plt.plot(mus, data[i], linewidth=1, label=r'$r_{11}$ = %s' % str(r11s[i]))

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':10})


plt.show()
