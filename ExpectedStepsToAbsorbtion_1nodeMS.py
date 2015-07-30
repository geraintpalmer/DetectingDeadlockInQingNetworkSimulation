from __future__ import division
import numpy as np
import yaml


# Input parameters here
n1 = 3
mu1 = 2.0
r11 = 0.5
L1 = 6.0
c1 = 1
directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/1Node_MultiServer/vary_c1/'

class Network:
    """
    A class to hold the queueing network object
    """

    def __init__(self, n1, mu1, r11, L1, c1):
        """
        Initialises the Network object

            >>> Q = Network(3, 10.0, 0.25, 7.0, 2)
            >>> Q.n1
            3
            >>> Q.mu1
            10.0
            >>> Q.r11
            0.25
            >>> Q.L1
            7.0
            >>> Q.c1
            2
            >>> Q.State_Space
            [0, 1, 2, 3, 4, 5, -1]
        """
        self.n1 = n1
        self.mu1 = mu1
        self.r11 = r11
        self.L1 = L1
        self.c1 = c1
        self.State_Space = [(i) for i in range(self.n1+(2*self.c1)+1)]
        self.write_transition_matrix()
        self.discretise_transition_matrix()

    def find_transition_rates(self, state1, state2):
        """
        Finds the transition rates for given state transition

            Testing all possible types of transitions
            >>> Q = Network(3, 10.0, 0.25, 7.0, 1)
            >>> Q.find_transition_rates(-1, -1)
            0
            >>> Q.find_transition_rates(-1, 3)
            0
            >>> Q.find_transition_rates(1, -1)
            0
            >>> Q.find_transition_rates(4, -1)
            2.5
            >>> Q.find_transition_rates(1, 2)
            7.0
            >>> Q.find_transition_rates(3, 2)
            7.5
        """
        delta = (state2-state1)
        if state1 < self.c1+self.n1:
            if delta==1:
                return self.L1
            if delta==-1:
                return (1-self.r11)*self.mu1*min(state1, self.c1)
        for j in range(self.c1+1):
            if state1==self.n1+self.c1+j:
                if delta==1:
                    return (self.c1-j)*self.r11*self.mu1
                if delta==-j-1:
                    return (1-self.r11)*(self.c1-j)*self.mu1
        return 0



    def write_transition_matrix(self):
        """
        Writes the transition matrix for the markov chain

            >>> Q = Network(3, 10.0, 0.25, 7.0, 2)
            >>> Q.write_transition_matrix()
            >>> Q.transition_matrix
            array([[ -7. ,   7. ,   0. ,   0. ,   0. ,   0. ,   0. ],
                   [  7.5, -14.5,   7. ,   0. ,   0. ,   0. ,   0. ],
                   [  0. ,  15. , -22. ,   7. ,   0. ,   0. ,   0. ],
                   [  0. ,   0. ,  15. , -22. ,   7. ,   0. ,   0. ],
                   [  0. ,   0. ,   0. ,  15. , -22. ,   7. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,  15. , -20. ,   5. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,  -0. ]])
        """
        self.transition_matrix = [[self.find_transition_rates(s1, s2) for s2 in self.State_Space] for s1 in self.State_Space]
        for i in range(len(self.transition_matrix)):
            a = sum(self.transition_matrix[i])
            self.transition_matrix[i][i] = -a
            self.transition_matrix = np.array(self.transition_matrix)

    def discretise_transition_matrix(self):
        """
        Disctetises the transition matrix

            >>> Q = Network(3, 10.0, 0.25, 7.0, 1)
            >>> Q.write_transition_matrix()
            >>> Q.discretise_transition_matrix()
            >>> Q.discrete_transition_matrix
            array([[ 0.51724138,  0.48275862,  0.        ,  0.        ,  0.        ,
                     0.        ],
                   [ 0.51724138,  0.        ,  0.48275862,  0.        ,  0.        ,
                     0.        ],
                   [ 0.        ,  0.51724138,  0.        ,  0.48275862,  0.        ,
                     0.        ],
                   [ 0.        ,  0.        ,  0.51724138,  0.        ,  0.48275862,
                     0.        ],
                   [ 0.        ,  0.        ,  0.        ,  0.51724138,  0.31034483,
                     0.17241379],
                   [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                     1.        ]])
        """
        self.time_step = 1 / max([abs(self.transition_matrix[i][i]) for i in range(len(self.transition_matrix))])
        self.discrete_transition_matrix = self.transition_matrix*self.time_step + np.identity(len(self.transition_matrix))

    def find_mean_time_to_absorbtion(self):
        """
        Finds the mean time to absorbtion

            >>> Q = Network(3, 10.0, 0.25, 7.0, 1)
            >>> Q.find_mean_time_to_absorbtion()
            >>> Q.mean_steps_to_absorbtion
            {'1': 53.625541441066225, '0': 55.696970012494788, '3': 42.665993336109956, '2': 49.334725114535608, '4': 33.449495002082472}
            >>> Q.mean_time_to_absorbtion
            {'1': 3.6983132028321535, '0': 3.841170345689296, '3': 2.942482299042066, '2': 3.4023948354852145, '4': 2.30686172428155}
        """
        T = self.discrete_transition_matrix[:-1, :-1]
        S = np.linalg.inv(np.identity(len(T)) - T)
        steps2absorb = [sum([S[i,j] for j in range(len(S))]) for i in range(len(S))]
        time2absorb = [s*self.time_step for s in steps2absorb]
        self.mean_steps_to_absorbtion = {str(self.State_Space[i]): steps2absorb[i] for i in range(len(steps2absorb))}
        self.mean_time_to_absorbtion = {str(self.State_Space[i]): float(time2absorb[i]) for i in range(len(time2absorb))}

    def write_results_to_file(self, param):
		"""
		Takes the summary statistics and writes them into a .yml file
		"""
		results_file = open('%stheoretical_results_%s.yml' % (directory, str(param)), 'w')
		results_file.write(yaml.dump(self.mean_time_to_absorbtion, default_flow_style=False))
		results_file.close()

if __name__ == '__main__':
    c1s = [1, 2, 3, 4, 5]
    for c1 in c1s:
        print "Starting c1 = " + str(c1)
        Q = Network(n1, mu1, r11, L1, c1)
        Q.find_mean_time_to_absorbtion()
        Q.write_results_to_file(c1)
