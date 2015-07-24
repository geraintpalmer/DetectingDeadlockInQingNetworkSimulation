"""
Usage: draw_graph_vary_L.py <dir_name> <state> <begin> <end> <step> <var>

Arguments
    dir_name    : name of the directory from which to read in data files
    state       : state we want to draw a graph of, e.g. (0, 0)
    begin       : number that the varying parameter starts on
    end         : number that the varying parameter ends on
    step        : the step size of the varying parameter
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


arguments = docopt.docopt(__doc__)
dirname = arguments['<dir_name>']
state = arguments['<state>']
begin = float(arguments['<begin>'])
end = float(arguments['<end>'])
step = float(arguments['<step>'])
var = str(arguments['<var>'])
root = os.getcwd()
directory = os.path.join(root, dirname)

params = [i*step + begin for i in range(int((end-begin)/step)+1)]

absorbtion_probs = []
for i in params:
    if var[0] == 'n':
        parameter_file_name = directory + 'absorbtion_probabilities_%s.yml' % str(int(i))
    else:
        parameter_file_name = directory + 'absorbtion_probabilities_%s.yml' % str(i)
    parameter_file = open(parameter_file_name, 'r')
    parameters = yaml.load(parameter_file)
    parameter_file.close()
    absorbtion_probs.append(parameters[state])


todeadlock1 = [ap[0] for ap in absorbtion_probs]
todeadlock2 = [ap[1] for ap in absorbtion_probs]
todeadlock3 = [ap[2] for ap in absorbtion_probs]

fig, ax = plt.subplots()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.78, box.height])

p1 = plt.bar(params, todeadlock1, width=step, color='red', label='Deadlock (-1)')
p2 = plt.bar(params, todeadlock2, width=step, color='cyan', bottom=todeadlock1, label='Deadlock (-2)')
p3 = plt.bar(params, todeadlock3, width=step, color='yellow', bottom=[todeadlock2[i]+todeadlock1[i] for i in range(len(todeadlock3))], label='Deadlock (-3)')

plt.title("Absorbtion Probabilities of Deadlocking from State (0, 0)")
ax.set_ylabel('Probability')
if var[0] == 'L':
    ax.set_xlabel(r'$\Lambda_{'+var[1:]+'}$')
if var[0] == 'm':
    ax.set_xlabel(r'$\mu_{'+var[2:]+'}$')
if var[0] == 'n':
    ax.set_xlabel(r'$n_{'+var[1:]+'}$')
if var[0] == 'r':
    ax.set_xlabel(r'$r_{'+var[1:]+'}$')
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

pylab.ylim([0,1])
pylab.xlim([begin,end+step])

plt.show()
