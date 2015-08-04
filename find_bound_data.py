import ExpectedStepsToAbsorbtion_1node
import ExpectedStepsToAbsorbtion_2NodeSimple
import ExpectedStepsToAbsorbtion_2NodeFeedback

directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/bound_data/'

def write_to_file(big_list, directory):
    """
    Writes the records for times to deadlock to a csv file
    """
    data_file = open('%sbound_data.csv' % directory, 'w')
    csv_wrtr = writer(data_file)
    for row in big_list:
        csv_wrtr.writerow(row)
    data_file.close()


LandMus = [float(i)/2.0 for i in range(31)][1:]
ns = range(11)
rs = [i/20.0 for i in range(41)][1:]

a_list = []
for L1 in LandMus:
	print "Beginning on L1 = " + str(L1)
	for L2 in LandMus:
		for mu1 in LandMus:
			print "Beginnning on mu1 = " + str(mu1) + " (L1 = " + str(L1) + ")"
			for mu2 in LandMus:
				for n1 in ns:
					for n2 in ns:
						for r11 in rs:
							for r12 in rs:
								for r21 in rs:
									for r22 in rs:
										if r11+r12 <= 1.0:
											if r21+r22 <= 1.0:
												Q11 = ExpectedStepsToAbsorbtion_1node.Network(n1, mu1, r11, L1)
												Q11.find_mean_time_to_absorbtion()
												Q12 = ExpectedStepsToAbsorbtion_1node.Network(n2, mu2, r22, L2)
												Q12.find_mean_time_to_absorbtion()
												Q2 = ExpectedStepsToAbsorbtion_2NodeSimple.Network(n1, n2, mu1, mu2, r12, r21, L1, L2)
												Q2.find_mean_time_to_absorbtion()
												Q3 = ExpectedStepsToAbsorbtion_2NodeFeedback.Network(n1, n2, mu1, mu2, r11, r12, r21, r22, L1, L2)
												Q3.find_mean_time_to_absorbtion()
												# order of row
												# [L1, L2, mu1, mu2, n1, n2, r11, r12, r21, r22, Q11, Q12, Q2, Q3]
												the_row = [L1, L2, mu1, mu2, n1, n2, r11, r12, r21, r22, Q11.mean_time_to_absorbtion['0'], Q12.mean_time_to_absorbtion['0'], Q2.mean_time_to_absorbtion['(0, 0)'], Q3.mean_time_to_absorbtion['(0, 0)']]
												a_list.append(the_row)

write_to_file(directory, a_list)