"""
Usage: draw_graph_vary_L.py <dir_name> <state> <begin> <end> <step> <var>

Arguments
    dir_name    : name of the directory from which to read in data files
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
import numpy as np

arguments = docopt.docopt(__doc__)
dirname = arguments['<dir_name>']
state = arguments['<state>']
begin = float(arguments['<begin>'])
end = float(arguments['<end>'])
step = float(arguments['<step>'])
var = str(arguments['<var>'])


Ls = [i*step + begin for i in range(int((end-begin)/step)+1)]

mean_times_to_deadlock = []
median_times_to_deadlock = []
simulation_results = []

root = os.getcwd()
directoryth = os.path.join(root, dirname)
directorys = os.path.join(root, dirname)

for i in Ls:
    parameter_file_name = directoryth + 'theoretical_results_mean_%s.yml' % str(i)
    parameter_file = open(parameter_file_name, 'r')
    parameters = yaml.load(parameter_file)
    parameter_file.close()
    mean_times_to_deadlock.append(parameters[state])

for i in Ls:
    parameter_file_name = directoryth + 'theoretical_results_median_%s.yml' % str(i)
    parameter_file = open(parameter_file_name, 'r')
    parameters = yaml.load(parameter_file)
    parameter_file.close()
    median_times_to_deadlock.append(parameters[state])

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



simulation_results_to_plot_no_outliers = []
for data_series in simulation_results_to_plot:
    q75, q25 = np.percentile(data_series, [75 ,25])
    iqr = q75 - q25
    upper = q75 + (1.5*iqr)
    lower = q25 - (1.5*iqr)
    new_data = [obs for obs in data_series if obs > lower and obs < upper]
    simulation_results_to_plot_no_outliers.append(new_data)



fig, ax = plt.subplots()
plt.plot(Ls, mean_times_to_deadlock, 'g')
plt.plot(Ls, median_times_to_deadlock, 'g')
plt.boxplot(simulation_results_to_plot, positions=Ls, widths=step/2, showmeans=True, sym='')
plt.violinplot(simulation_results_to_plot_no_outliers, widths=step/1.5, positions=Ls, showmeans=False, showmedians=False, showextrema=False)
ax.set_xlabel(var)
ax.set_ylabel('Time to Deadlock from (0, 0)')
ax.set_title('Expected Time to Deadlock From State (0, 0)')
fig.tight_layout()
plt.show()
