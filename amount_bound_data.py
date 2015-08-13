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

print len(huge_list)