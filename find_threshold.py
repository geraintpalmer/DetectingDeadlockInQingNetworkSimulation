import ExpectedStepsToAbsorbtion_1node
from csv import writer

def loop_through(startmu, step, n1, r11, L1):
	deadlock_list = ['Inf']
	prevmu = 0.0
	i = startmu
	go = True
	while go:
		mu1 = i
		Q = ExpectedStepsToAbsorbtion_1node.Network(n1, mu1, r11, L1)
		Q.find_mean_time_to_absorbtion()
		deadlock_list.append(Q.mean_time_to_absorbtion['0'])
		if deadlock_list[-1] > deadlock_list[-2]:
			go = False
			actualmu = prevmu
		else:
			prevmu = i
			i += step
	return actualmu



def find_threshold_to_5dc(start, n1, r11, L1):
	s1 = loop_through(start, 0.1, n1, r11, L1)
	s2 = loop_through(s1, 0.01, n1, r11, L1)
	s3 = loop_through(s2, 0.001, n1, r11, L1)
	s4 = loop_through(s3, 0.0001, n1, r11, L1)
	s5 = loop_through(s4, 0.00001, n1, r11, L1)
	s6 = loop_through(s5, 0.000001, n1, r11, L1)
	return s5

def write_to_file(threshold_list, var_list, directory, sffx):
	data_file = open('%sthresholds_%s.csv' % (directory, sffx), 'w')
	csv_wrtr = writer(data_file)
	csv_wrtr.writerow(threshold_list)
	csv_wrtr.writerow(var_list)
	data_file.close()


if __name__ == '__main__':
	directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/1Node/thresholds/'
	L1 = 4.0
	n1 = 3
	r11 = 0.25
	start = 0.0001
	threshold_list = []
	var_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	for n1 in var_list:
		threshold_list.append(find_threshold_to_5dc(start, n1, r11, L1))
	write_to_file(threshold_list, var_list, directory, 'n1')