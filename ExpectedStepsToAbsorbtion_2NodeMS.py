from __future__ import division
import numpy as np
import yaml
import random
import csv


# # Input parameters here
n1, n2 = 2, 1
c1, c2 = 2, 2
mu1, mu2 = 5.5, 6.5
r12, r21 = 0.7, 0.6
L1, L2 = 9.0, 7.5
directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/2NodeSimple_MultiServer/run_10000_itrs/vary_c2/'


class Network:
    """
    A class to hold the queueing network object
    """

    def __init__(self, n1, n2, mu1, mu2, r12, r21, L1, L2, c1, c2):
        """
        Initialises the Network object

            >>> Q = Network(3, 2, 10.0, 8.0, 0.4, 0.3, 6.0, 5.0, 2, 2)
            >>> Q.n1
            3
            >>> Q.n2
            2
            >>> Q.mu1
            10.0
            >>> Q.mu2
            8.0
            >>> Q.r12
            0.4
            >>> Q.r21
            0.3
            >>> Q.L1
            6.0
            >>> Q.L2
            5.0
            >>> Q.c1
            2
            >>> Q.c2
            2
            >>> Q.State_Space
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (4, 0), (4, 1), (4, 2), (4, 3), (5, 0), (5, 1), (5, 2), -1]
        """
        self.n1 = n1
        self.n2 = n2
        self.mu1 = mu1
        self.mu2 = mu2
        self.r12 = r12
        self.r21 = r21
        self.L1 = L1
        self.L2 = L2
        self.c1 = c1
        self.c2 = c2
        self.State_Space = [(i, j) for i in range(self.n1+self.c1+self.c2+1) for j in range(self.n2+self.c2+self.c1+1) if i<=self.n1+self.c1+j if j<=self.n2+c2+i]
        self.write_transition_matrix()
        self.discretise_transition_matrix()

    def find_transition_rates(self, state1, state2):
        """
        Finds the transition rates for given state transition

            Testing all possible types of transitions
            >>> Q = Network(3, 2, 10.0, 8.0, 0.4, 0.3, 6.0, 5.0)
            >>> Q.find_transition_rates((1, 1), -1)
            0
            >>> Q.find_transition_rates(-1, (4, 2))
            0
            >>> Q.find_transition_rates((3, 4), -1)
            2.4
            >>> Q.find_transition_rates((5, 2), -1)
            4.0
            >>> Q.find_transition_rates((2, 2), (3, 2))
            6.0
            >>> Q.find_transition_rates((2, 1), (2, 2))
            5.0
            >>> Q.find_transition_rates((4, 1), (5, 1))
            0
            >>> Q.find_transition_rates((1, 3), (1, 4))
            0
            >>> Q.find_transition_rates((1, 2), (2, 2))
            6.0
            >>> Q.find_transition_rates((3, 2), (3, 1))
            5.6
            >>> Q.find_transition_rates((2, 2), (1, 3))
            4.0
            >>> Q.find_transition_rates((2, 2), (3, 1))
            2.4
            >>> Q.find_transition_rates((2, 4), (1, 5))
            0
            >>> Q.find_transition_rates((5, 2), (6, 1))
            0
        """
        delta = (state2[0] - state1[0], state2[1] - state1[1])
        k1 = max(0, state1[0]-(self.n1+self.c1))
        k2 = max(0, state1[1]-(self.n2+self.c2))
        if state1[0] < self.n1 + self.c1:
            if state1[1] < self.n2 + self.c2:
                if delta == (1, 0):
                    return self.L1
                if delta == (0, 1):
                    return self.L2
                if delta == (-1, 1):
                    return self.r12*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (1, -1):
                    return self.r21*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (-1, 0):
                    return (1-self.r12)*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (0, -1):
                    return (1-self.r21)*self.mu2*min(state1[1]-k1, self.c2-k1)
            if state1[1] == self.n2 + self.c2:
                if delta == (1, 0):
                    return self.L1
                if delta == (1, -1):
                    return self.r21*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (0, 1):
                    return self.r12*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (-1, 0):
                    return (1-self.r12)*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (0, -1):
                    return (1-self.r21)*self.mu2*min(state1[1]-k1, self.c2-k1)
            if state1[1] > self.n2 + self.c2:
                if delta == (1, 0):
                    return self.L1
                if delta == (0, -1):
                    return self.r21*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (0, 1):
                    return self.r12*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (-1, -1):
                    return (1-self.r21)*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (-1, 0):
                    return (1-self.r12)*self.mu1*min(state1[0]-k2, self.c1-k2)
        if state1[0] == self.n1 + self.c1:
            if state1[1] < self.n2 + self.c2:
                if delta == (0, 1):
                    return self.L2
                if delta == (-1, 1):
                    return self.r12*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (1, 0):
                    return self.r21*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (-1, 0):
                    return (1-self.r12)*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (0, -1):
                    return (1-self.r21)*self.mu2*min(state1[1]-k1, self.c2-k1)
            if state1[1] == self.n2 + self.c2:
                if delta == (0, 1):
                    return self.r12*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (1, 0):
                    return self.r21*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (-1, 0):
                    return (1-self.r12)*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (0, -1):
                    return (1-self.r21)*self.mu2*min(state1[1]-k1, self.c2-k1)
            if state1[1] > self.n2 + self.c2:
                if delta == (0, 1):
                    return self.r12*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (1, 0):
                    return self.r21*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (-1, 0):
                    return (1-self.r12)*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (-1, -1):
                    return (1-self.r21)*self.mu2*min(state1[1]-k1, self.c2-k1)
        if state1[0] > self.n1 + self.c1:
            if state1[1] < self.n2 + self.c2:
                if delta == (0, 1):
                    return self.L2
                if delta == (-1, 0):
                    return self.r12*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (1, 0):
                    return self.r21*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (-1, -1):
                    return (1-self.r12)*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (0, -1):
                    return (1-self.r21)*self.mu2*min(state1[1]-k1, self.c2-k1)
            if state1[1] == self.n2 + self.c2:
                if delta == (0, 1):
                    return self.r12*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (1, 0):
                    return self.r21*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (0, -1):
                    return (1-self.r21)*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (-1, -1):
                    return (1-self.r12)*self.mu1*min(state1[0]-k2, self.c1-k2)
            if state1[1] > self.n2 + self.c2:
                if delta == (1, 0):
                    return self.r21*self.mu2*min(state1[1]-k1, self.c2-k1)
                if delta == (0, 1):
                    return self.r12*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (-min(k1+1,k2+1), -min(k1,k2+1)):
                    return (1-self.r12)*self.mu1*min(state1[0]-k2, self.c1-k2)
                if delta == (-min(k1+1,k2), -min(k1+1,k2+1)):
                    return (1-self.r21)*self.mu2*min(state1[1]-k1, self.c2-k1)
        return 0

    def write_transition_matrix(self):
        """
        Writes the transition matrix for the markov chain

            >>> Q = Network(1, 1, 6.0, 7.0, 0.1, 0.1, 3.0, 4.0)
            >>> Q.write_transition_matrix()
            >>> Q.transition_matrix
            array([[ -7. ,   4. ,   0. ,   0. ,   3. ,   0. ,   0. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ],
                   [  6.3, -14. ,   4. ,   0. ,   0.7,   3. ,   0. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ],
                   [  0. ,   6.3, -10. ,   0. ,   0. ,   0.7,   3. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ],
                   [  0. ,   0. ,   6.3, -10. ,   0. ,   0. ,   0.7,   3. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ],
                   [  5.4,   0.6,   0. ,   0. , -13. ,   4. ,   0. ,   0. ,   3. ,
                      0. ,   0. ,   0. ,   0. ,   0. ],
                   [  0. ,   5.4,   0.6,   0. ,   6.3, -20. ,   4. ,   0. ,   0.7,
                      3. ,   0. ,   0. ,   0. ,   0. ],
                   [  0. ,   0. ,   5.4,   0.6,   0. ,   6.3, -16. ,   0. ,   0. ,
                      0.7,   3. ,   0. ,   0. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   6.3,  -7.7,   0. ,
                      0. ,   0.7,   0. ,   0. ,   0.7],
                   [  0. ,   0. ,   0. ,   0. ,   5.4,   0.6,   0. ,   0. , -10. ,
                      4. ,   0. ,   0. ,   0. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   5.4,   0.6,   0. ,   6.3,
                    -17. ,   4. ,   0.7,   0. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   5.4,   0.6,   0. ,
                      6.3, -13. ,   0. ,   0.7,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   5.4,
                      0.6,   0. , -10. ,   4. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,
                      5.4,   0.6,   0. ,  -6.6,   0.6],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,  -0. ]])
        """
        self.transition_matrix = [[self.find_transition_rates(s1, s2) for s2 in self.State_Space] for s1 in self.State_Space]
        for i in range(len(self.transition_matrix)):
            a = sum(self.transition_matrix[i])
            self.transition_matrix[i][i] = -a
            self.transition_matrix = np.array(self.transition_matrix)

    def discretise_transition_matrix(self):
        """
        Disctetises the transition matrix

            >>> Q = Network(1, 1, 6.0, 7.0, 0.1, 0.1, 3.0, 4.0)
            >>> Q.write_transition_matrix()
            >>> Q.discretise_transition_matrix()
            >>> Q.discrete_transition_matrix
            array([[ 0.65 ,  0.2  ,  0.   ,  0.   ,  0.15 ,  0.   ,  0.   ,  0.   ,
                     0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
                   [ 0.315,  0.3  ,  0.2  ,  0.   ,  0.035,  0.15 ,  0.   ,  0.   ,
                     0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
                   [ 0.   ,  0.315,  0.5  ,  0.   ,  0.   ,  0.035,  0.15 ,  0.   ,
                     0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
                   [ 0.   ,  0.   ,  0.315,  0.5  ,  0.   ,  0.   ,  0.035,  0.15 ,
                     0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
                   [ 0.27 ,  0.03 ,  0.   ,  0.   ,  0.35 ,  0.2  ,  0.   ,  0.   ,
                     0.15 ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
                   [ 0.   ,  0.27 ,  0.03 ,  0.   ,  0.315,  0.   ,  0.2  ,  0.   ,
                     0.035,  0.15 ,  0.   ,  0.   ,  0.   ,  0.   ],
                   [ 0.   ,  0.   ,  0.27 ,  0.03 ,  0.   ,  0.315,  0.2  ,  0.   ,
                     0.   ,  0.035,  0.15 ,  0.   ,  0.   ,  0.   ],
                   [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.315,  0.615,
                     0.   ,  0.   ,  0.035,  0.   ,  0.   ,  0.035],
                   [ 0.   ,  0.   ,  0.   ,  0.   ,  0.27 ,  0.03 ,  0.   ,  0.   ,
                     0.5  ,  0.2  ,  0.   ,  0.   ,  0.   ,  0.   ],
                   [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.27 ,  0.03 ,  0.   ,
                     0.315,  0.15 ,  0.2  ,  0.035,  0.   ,  0.   ],
                   [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.27 ,  0.03 ,
                     0.   ,  0.315,  0.35 ,  0.   ,  0.035,  0.   ],
                   [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                     0.27 ,  0.03 ,  0.   ,  0.5  ,  0.2  ,  0.   ],
                   [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                     0.   ,  0.27 ,  0.03 ,  0.   ,  0.67 ,  0.03 ],
                   [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                     0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  1.   ]])
        """
        self.time_step = 1 / max([abs(self.transition_matrix[i][i]) for i in range(len(self.transition_matrix))])
        self.discrete_transition_matrix = self.transition_matrix*self.time_step + np.identity(len(self.transition_matrix))

    def find_mean_time_to_absorbtion(self):
        """
        Finds the mean time to absorbtion

            >>> Q = Network(1, 1, 6.0, 7.0, 0.1, 0.1, 3.0, 4.0)
            >>> Q.find_mean_time_to_absorbtion()
            >>> Q.mean_steps_to_absorbtion
            {'(0, 1)': 3574.6221167527738, '(3, 0)': 3433.1460274057613, '(0, 0)': 3576.9015539371817, '(1, 0)': 3573.2741368497277, '(3, 1)': 3230.0029866619047, '(1, 1)': 3568.1378355791135, '(2, 1)': 3552.5015527639657, '(1, 2)': 3557.8317616860791, '(2, 0)': 3566.6569251391866, '(0, 3)': 3471.0401728417651, '(1, 3)': 3233.9311662561877, '(0, 2)': 3571.1311105506093, '(2, 2)': 3524.1855450719122}
            >>> Q.mean_time_to_absorbtion
            {'(0, 1)': 178.7311058376387, '(3, 0)': 171.65730137028808, '(0, 0)': 178.8450776968591, '(1, 0)': 178.6637068424864, '(3, 1)': 161.50014933309524, '(1, 1)': 178.4068917789557, '(2, 1)': 177.6250776381983, '(1, 2)': 177.89158808430398, '(2, 0)': 178.33284625695933, '(0, 3)': 173.55200864208825, '(1, 3)': 161.6965583128094, '(0, 2)': 178.55655552753046, '(2, 2)': 176.20927725359562}
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

    def write_matrix_to_file(self):
        """
        Takes the summary statistics and writes them into a .yml file
        """
        results_file = open('%stransmat.csv' % directory, 'w')
        np.savetxt(results_file, Q.transition_matrix, delimiter=",")


if __name__ == '__main__':
    # mu2s = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0]
    c2s = [1, 2, 3, 4, 5]
    # L2s = [4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0]
    # r12s = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    for c2 in c2s:
        Q = Network(n1, n2, mu1, mu2, r12, r21, L1, L2, c1, c2)
        Q.find_mean_time_to_absorbtion()
        Q.write_results_to_file(c2)
        print "Now starting c2 = " + str(c2)
    # Q = Network(n1, n2, mu1, mu2, r12, r21, L1, L2, c1, c2)
    # for state in Q.State_Space:
    #     print state
    # # for row in Q.transition_matrix:
    # #     print row
    # # print Q.transition_matrix[-1:]
    # # Q.write_matrix_to_file()