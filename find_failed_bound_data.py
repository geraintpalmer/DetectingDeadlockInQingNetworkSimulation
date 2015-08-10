from __future__ import division
import numpy as np
from csv import reader, writer

directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/bound_data/'

huge_list = []
data_file = open('%sbound_data.csv' % directory, 'r')
rdr = reader(data_file)
for row in rdr:
	huge_list.append(row)
data_file.close()


failed_list = []
for row in huge_list:
	omega = float(row[-1])
	minomega = min(float(row[-2]), max(float(row[-3]), float(row[-4])), max(float(row[-5]), float(row[-6])))
	if omega > minomega:
		failed_list.append(row)


data_file = open('%sfailed_bound_data.csv' % directory, 'w')
csv_wrtr = writer(data_file)
for row in failed_list:
    csv_wrtr.writerow(row)
data_file.close()

