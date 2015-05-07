"""
Usage: draw_graph_vary_L.py <dir_name> <state> <begin> <end> <step> <var>

Arguments
    dir_name_th : name of the directory from which to read in data files (theoretical results)
    dir_name_s  : name of the directory from which to read in data files (simulation results)
    state       : state we want to draw a graph of, e.g. (0, 0)
    begin       : number that the varying parameter starts on
    end         : number that the varying parameter ends on
    step        : the step size of the varying parameter

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

arguments = docopt.docopt(__doc__)
dirname = arguments['<dir_name>']
state = arguments['<state>']
begin = float(arguments['<begin>'])
end = float(arguments['<end>'])
step = float(arguments['<step>'])
var = str(arguments['<var>'])


Ls = [i*step + begin for i in range(int((end-begin)/step)+1)]

times_to_deadlock = []
simulation_results = []

root = os.getcwd()
directoryth = os.path.join(root, dirname)
directorys = os.path.join(root, dirname)

for i in Ls:
    parameter_file_name = directoryth + 'theoretical_results_%s.yml' % str(i)
    parameter_file = open(parameter_file_name, 'r')
    parameters = yaml.load(parameter_file)
    parameter_file.close()
    times_to_deadlock.append(parameters[state])

for i in Ls:
    parameter_file_name = directorys + 'deadlocking_times_%s.csv' % str(i)
    parameter_file = open(parameter_file_name, 'r')
    rdr = reader(parameter_file)
    for row in rdr:
        if row[0] == "((0, 0), (0, 0))":
            simulation_results.append(row[1:])
            break
    parameter_file.close()

i = 0
simulation_results_to_plot = []
for row in simulation_results:
    simulation_results_to_plot.append([])
    for obs in row:
        simulation_results_to_plot[i].append(float(obs))
    i += 1



fig, ax = plt.subplots()
plt.plot(Ls, times_to_deadlock)
plt.boxplot(simulation_results_to_plot, positions=Ls, widths=step/2, showmeans=True, sym='')
ax.set_xlabel(var)
ax.set_ylabel('Time to Deadlock from (0, 0)')
ax.set_title('Expected Time to Deadlock From State (0, 0)')
fig.tight_layout()
plt.show()
