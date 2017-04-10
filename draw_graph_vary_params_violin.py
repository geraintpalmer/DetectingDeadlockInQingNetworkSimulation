"""
Usage: draw_graph_vary_L.py <dir_name> <stateth> <state> <begin> <end> <step> <var> <legend_loc> <line_colour> <plot_colour>

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
import seaborn as sns
sns.set(style="whitegrid")


arguments = docopt.docopt(__doc__)
dirname = arguments['<dir_name>']
stateth = arguments['<stateth>']
state = arguments['<state>']
begin = float(arguments['<begin>'])
end = float(arguments['<end>'])
step = float(arguments['<step>'])
var = str(arguments['<var>'])
leg_loc = str(arguments['<legend_loc>'])
line_colour = str(arguments['<line_colour>'])
plot_colour = str(arguments['<plot_colour>'])

Ls = [round(i*step + begin, 2) for i in range(int((end-begin)/step)+1)]

mean_times_to_deadlock = []
median_times_to_deadlock = []
simulation_results = []


root = os.getcwd()
directoryth = os.path.join(root, dirname)
directorys = os.path.join(root, dirname)

for i in Ls:
    if var[0] == 'n' or var[0] == 'c':
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
    if var[0] == 'n' or var[0] == 'c':
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

if var[0] == 'c' or var[0] == 'n':
    Ls_show = [int(obs) for obs in Ls]
else:
    Ls_show = Ls


fig, ax = plt.subplots()
plt.plot(Ls_show, mean_times_to_deadlock, linewidth=2, label='Analytical Mean', color=line_colour)
plt.plot([], [], 'r', linewidth=2, label='Simulation Results', color=plot_colour)
meanpointprops = dict(marker='D', markeredgecolor='black',
                      markerfacecolor=plot_colour)
bp = ax.boxplot(simulation_results, positions=Ls_show, widths=step/2, meanprops=meanpointprops, meanline=False, showmeans=True, sym='')
for median in bp['medians']:
    median.set(color=plot_colour, linewidth=2)
pylab.setp(bp['boxes'], color=plot_colour)
pylab.setp(bp['whiskers'], color=plot_colour)


vp = plt.violinplot(simulation_results_no_outliers, widths=step/1.5, positions=Ls, showmeans=False, showmedians=False, showextrema=False)
pylab.setp(vp['bodies'], color=plot_colour)

if var[0] == 'L':
    ax.set_xlabel(r'$\Lambda_{'+var[1:]+'}$', fontsize=18)
if var[0] == 'm':
    ax.set_xlabel(r'$\mu_{'+var[2:]+'}$', fontsize=18)
if var[0] == 'n':
    ax.set_xlabel(r'$n_{'+var[1:]+'}$', fontsize=18)
if var[0] == 'c':
    ax.set_xlabel(r'$c_{'+var[1:]+'}$', fontsize=18)
if var[0] == 'r':
    ax.set_xlabel(r'$r_{'+var[1:]+'}$', fontsize=18)

ax.set_ylabel('Time to Deadlock from ' + stateth, fontsize=18)
# ax.set_title('Expected Time to Deadlock From State (0, 0)', fontsize=20)

ax.set_xlim(begin - step/2, end + step/2)
ax.set_ylim(bottom=0)

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

if leg_loc == 'l':
    loc=2
if leg_loc == 'r':
    loc=1

legend = plt.legend(loc=loc, prop={'size':20}, frameon=1)
frame = legend.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('grey')


plt.savefig('newplots_paper/2Nmsfb_vary' + var + '.pdf')