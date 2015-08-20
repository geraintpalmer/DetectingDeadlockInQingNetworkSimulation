"""
Usage: dual_axis_ns.py <dir_name> <state> <ran>

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

arguments = docopt.docopt(__doc__)
directory = arguments['<dir_name>']
state = arguments['<state>']
ran = int(arguments['<ran>'])

node = directory[-2]

mean_times_to_deadlock = []
for i in range(ran):
    parameter_file_name = directory + 'theoretical_results_%s.yml' % str(int(i))
    parameter_file = open(parameter_file_name, 'r')
    parameters = yaml.load(parameter_file)
    parameter_file.close()
    mean_times_to_deadlock.append(parameters[state])

absorbtion_probs = []
for i in range(ran):
    parameter_file_name = directory + 'absorbtion_probabilities_%s.yml' % str(int(i))
    parameter_file = open(parameter_file_name, 'r')
    parameters = yaml.load(parameter_file)
    parameter_file.close()
    absorbtion_probs.append(parameters[state])





fig, ax1 = plt.subplots()

ax1.set_xlabel(r'$n_{'+node+'}$')
ax1.set_ylabel('Mean Time to Deadlock', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')
ax1.plot(range(ran), mean_times_to_deadlock, lw=1, color='b')

if node == '1':
	indx = 1
if node == '2':
	indx = 0

ax2 = ax1.twinx()
# ax2.bar(range(ran), [absorbtion_probs[i][indx] for i in range(len(absorbtion_probs))], width=0.5, color='g', alpha=0.6)
ax2.plot(range(ran), [absorbtion_probs[i][indx] for i in range(len(absorbtion_probs))], lw=1, color='g')
ax2.set_ylabel('Probability to Deadlock State (-' + str(indx+1) + ')', color='g')
for tl in ax2.get_yticklabels():
    tl.set_color('g')

ax1.set_title('Comparison of Time to Deadlock and Absorption Probabiliy - Varying ' + r'$n_{'+node+'}$')

plt.show()