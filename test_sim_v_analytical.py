"""
Usage: test_sim_v_analytical.py <dir_name> <state> <begin> <end> <step>

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
from scipy import stats



class Statistical_Tests:
	"""
	A class to perform statistical tests on data
	"""
	def __init__(self, dirname, state, begin, end, step):
		"""
		Initialises the class

			>>> ST = Statistical_Tests('data_for_graphs/2NodeSimple/run_10000_itrs/vary_L1/', '(0, 0)', 2.5, 14.0, 0.5)
			>>> ST.dirname
			'data_for_graphs/2NodeSimple/run_10000_itrs/vary_L1/'
			>>> ST.state
			'(0, 0)'
			>>> ST.begin
			2.5
			>>> ST.end
			14.0
			>>> ST.step
			0.5
			>>> ST.var
			'L1'
			>>> ST.Ls
			[2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0]
			>>> ST.alpha
			0.05
		"""
		self.dirname = dirname
		self.state = state
		self.begin = begin
		self.end = end
		self.step = step
		self.var = self.dirname.split('_')[-1][:-1]
		self.Ls = [i*step + begin for i in range(int((end-begin)/step)+1)]
		self.alpha = 0.05
		self.load_data()

	def load_data(self):
		"""
		Loads the data into the code

			>>> ST = Statistical_Tests('data_for_graphs/2NodeSimple/run_10000_itrs/vary_L1/', '(0, 0)', 2.5, 14.0, 0.5)
			>>> ST.times_to_deadlock[0]
			796.0852043581275
			>>> ST.simulation_results_floats[0][0]
			2170.229871430381
		"""
		self.times_to_deadlock = []
		self.simulation_results = []
		root = os.getcwd()
		directory = os.path.join(root, self.dirname)

		for i in self.Ls:
			parameter_file_name = directory + 'theoretical_results_%s.yml' % str(i)
			parameter_file = open(parameter_file_name, 'r')
			parameters = yaml.load(parameter_file)
			parameter_file.close()
			self.times_to_deadlock.append(parameters[self.state])

		for i in self.Ls:
			parameter_file_name = directory + 'deadlocking_times_%s.csv' % str(i)
			parameter_file = open(parameter_file_name, 'r')
			rdr = reader(parameter_file)
			for row in rdr:
				if row[0] == "((0, 0), (0, 0))":
					self.simulation_results.append(row[1:])
					break
			parameter_file.close()

		i = 0
		self.simulation_results_floats = []
		for row in self.simulation_results:
			self.simulation_results_floats.append([])
			for obs in row:
				self.simulation_results_floats[i].append(float(obs))
			i += 1

	def compare_means(self):
		"""
		Compares the means of the simulation results and the means derived analytically
		"""
		self.stats_tests_results = {}
		for i in range(len(self.Ls)):
			my_key = var + str(self.Ls[i])
			tstat, pval = stats.ttest_1samp(self.simulation_results_floats[i], self.times_to_deadlock[i])
			ttest = pval > self.alpha			
			self.stats_tests_results[my_key] = [ttest, pval]


if __name__ == '__main__':
	arguments = docopt.docopt(__doc__)
	dirname = arguments['<dir_name>']
	state = arguments['<state>']
	begin = float(arguments['<begin>'])
	end = float(arguments['<end>'])
	step = float(arguments['<step>'])
	var = dirname.split('_')[-1][:-1]
	ST = Statistical_Tests(dirname, state, begin, end, step)
	ST.compare_means()
	print ST.stats_tests_results

