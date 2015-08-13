"""
Usage: experiment_parallelprocessing.py <LandMs> <nis> <ris>

Arguments
    LandMs 	: list of lambdas and mus
    nis		: list of ns
    ris 	: list of rs

Options
    -h          : displays this help file
"""
import ExpectedStepsToAbsorbtion_1node
import ExpectedStepsToAbsorbtion_2NodeSimple
import ExpectedStepsToAbsorbtion_2NodeFeedback
import docopt
from csv import writer, reader
from math import exp

directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/bound_data/'


arguments = docopt.docopt(__doc__)
LandMs = arguments['<LandMs>']
nis = arguments['<nis>']
ris = arguments['<ris>']


def write_to_file(big_list, directory):
    """
    Writes the records for times to deadlock to a csv file
    """
    huge_list = []
    data_file = open('%sbound_data.csv' % directory, 'r')
    rdr = reader(data_file)
    for row in rdr:
    	huge_list.append(row)
    data_file.close()

    for row in big_list:
    	huge_list.append(row)

    data_file = open('%sbound_data.csv' % directory, 'w')
    csv_wrtr = writer(data_file)
    for row in huge_list:
        csv_wrtr.writerow(row)
    data_file.close()

def find_p_full(L, mu, n):
 	"""
 	Finds the probability of a ode beign full
 	"""
 	n += 1
 	if L == mu:
 		return 1.0/(n+1)
 	else:
 		rho = L/mu
 		return (((1-rho)*(rho**n))/(1-(rho**(n+1))))

def find_p_empty(L, mu, n):
 	"""
 	Finds the probability of a ode beign full
 	"""
 	n += 1
 	if L == mu:
 		return 1.0/(n+1)
 	else:
 		rho = L/mu
 		return ((1-rho)/(1-(rho**(n+1))))


LandMus = eval(LandMs)
ns = eval(nis)
rs = eval(ris)

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
												mu1_dash = (mu2/(mu1+mu2))/((2.0/mu1)+(1.0/mu2))
												mu2_dash = (mu1/(mu2+mu1))/((1.0/mu1)+(2.0/mu2))
												Q11_low = ExpectedStepsToAbsorbtion_1node.Network(n1, mu1_dash, r11, L1)
												Q11_low.find_mean_time_to_absorbtion()
												Q12_low = ExpectedStepsToAbsorbtion_1node.Network(n2, mu2_dash, r22, L2)
												Q12_low.find_mean_time_to_absorbtion()
												Q11_high = ExpectedStepsToAbsorbtion_1node.Network(n1, mu1, r11, L1)
												Q11_high.find_mean_time_to_absorbtion()
												Q12_high = ExpectedStepsToAbsorbtion_1node.Network(n2, mu2, r22, L2)
												Q12_high.find_mean_time_to_absorbtion()
												Q2 = ExpectedStepsToAbsorbtion_2NodeSimple.Network(n1, n2, mu1, mu2, r12, r21, L1, L2)
												Q2.find_mean_time_to_absorbtion()
												Q3 = ExpectedStepsToAbsorbtion_2NodeFeedback.Network(n1, n2, mu1, mu2, r11, r12, r21, r22, L1, L2)
												Q3.find_mean_time_to_absorbtion()
												# order of row
												# [L1, L2, mu1, mu2, n1, n2, r11, r12, r21, r22, Q11_low, Q11_high Q12_low, Q12_high, Q2, Q3]
												the_row = [L1, L2, mu1, mu2, n1, n2, r11, r12, r21, r22, Q11_low.mean_time_to_absorbtion['0'], Q11_high.mean_time_to_absorbtion['0'], Q12_low.mean_time_to_absorbtion['0'], Q12_high.mean_time_to_absorbtion['0'], Q2.mean_time_to_absorbtion['(0, 0)'], Q3.mean_time_to_absorbtion['(0, 0)']]
												a_list.append(the_row)

write_to_file(a_list, directory)