"""
Usage: draw_graph_vary_L.py <dir_name> <stateth> <state> <begin> <end> <step> <var> <legend_loc>

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
stateth = arguments['<stateth>']
state = arguments['<state>']
begin = float(arguments['<begin>'])
end = float(arguments['<end>'])
step = float(arguments['<step>'])
var = str(arguments['<var>'])
leg_loc = str(arguments['<legend_loc>'])

Ls = [i*step + begin for i in range(int((end-begin)/step)+1)]

mean_times_to_deadlock = []
median_times_to_deadlock = []
simulation_results = []


root = os.getcwd()
directoryth = os.path.join(root, dirname)
directorys = os.path.join(root, dirname)

for i in Ls:
    if var[0] == 'n':
        parameter_file_name = directoryth + 'theoretical_results_%s.yml' % str(int(i))
    else:
        parameter_file_name = directoryth + 'theoretical_results_%s.yml' % str(i)
    parameter_file = open(parameter_file_name, 'r')
    parameters = yaml.load(parameter_file)
    parameter_file.close()
    mean_times_to_deadlock.append(parameters[stateth])

# for i in Ls:
#     parameter_file_name = directoryth + 'theoretical_results_median_%s.yml' % str(i)
#     parameter_file = open(parameter_file_name, 'r')
#     parameters = yaml.load(parameter_file)
#     parameter_file.close()
#     median_times_to_deadlock.append(parameters[state])

for i in Ls:
    if var[0] == 'n':
        parameter_file_name = directorys + 'deadlocking_times_%s.csv' % str(int(i))
    else:
        parameter_file_name = directorys + 'deadlocking_times_%s.csv' % str(i)
    parameter_file = open(parameter_file_name, 'r')
    rdr = reader(parameter_file)
    for row in rdr:
        if row[0] == state:
            simulation_results.append([float(obs) for obs in row[1:]])
            break
    parameter_file.close()


simulation_results_no_outliers = []
for data_series in simulation_results:
    q75, q25 = np.percentile(data_series, [75 ,25])
    iqr = q75 - q25
    upper = q75 + (1.5*iqr)
    lower = q25 - (1.5*iqr)
    new_data = [obs for obs in data_series if obs > lower and obs < upper]
    simulation_results_no_outliers.append(new_data)



fig, ax = plt.subplots()
plt.plot(Ls, mean_times_to_deadlock, linewidth=2, label='Analytical Mean')
plt.plot([], [], 'r', linewidth=2, label='Simulation Results')

bp = plt.boxplot(simulation_results, positions=Ls, widths=step/2, showmeans=True, sym='')
pylab.setp(bp['boxes'], color='red')
pylab.setp(bp['whiskers'], color='red')

vp = plt.violinplot(simulation_results_no_outliers, widths=step/1.5, positions=Ls, showmeans=False, showmedians=False, showextrema=False)
pylab.setp(vp['bodies'], color='red')

if var[0] == 'L':
    ax.set_xlabel(r'$\Lambda_{'+var[1:]+'}$')
if var[0] == 'm':
    ax.set_xlabel(r'$\mu_{'+var[2:]+'}$')
if var[0] == 'n':
    ax.set_xlabel(r'$n_{'+var[1:]+'}$')
if var[0] == 'r':
    ax.set_xlabel(r'$r_{'+var[1:]+'}$')

ax.set_ylabel('Time to Deadlock from (0, 0)')
ax.set_title('Expected Time to Deadlock From State (0, 0)')

if leg_loc == 'l':
    plt.legend(loc=2)
if leg_loc == 'r':
    plt.legend(loc=1)
plt.show()