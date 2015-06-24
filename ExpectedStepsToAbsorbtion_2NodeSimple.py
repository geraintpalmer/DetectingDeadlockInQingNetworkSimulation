"""
This file is a simulation
"""
from __future__ import division
import numpy as np
import yaml
import random


# Input parameters here
n1, n2 = 3, 2
mu1, mu2 = 10.0, 8.0
r12, r21 = 0.25, 0.15
L1, L2 = 4.0, 5.0
directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/2NodeSimple/run_10000_itrs/vary_r12/'

class Network:
    """
    A class to hold the queueing network object
    """

    def __init__(self, n1, n2, mu1, mu2, r12, r21, L1, L2):
        """
        Initialises the Network object

            >>> Q = Network(3, 2, 10.0, 8.0, 0.4, 0.3, 6.0, 5.0)
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
        self.State_Space = [(i, j) for i in range(self.n1+3) for j in range(self.n2+3) if i+j<=self.n1+self.n2+2] + [-1]
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
        if state1 == -1:
            return 0
        if state2 == -1:
            if state1[0] == self.n1 and state1[1] == self.n2 + 2:
                return self.r21*self.mu2
            if state1[0] == self.n1 + 2 and state1[1] == self.n2:
                return self.r12*self.mu1
            else:
                return 0
        else:
            delta = (state2[0]-state1[0], state2[1]-state1[1])
            if delta == (1, 0):
                if state1[0]<=self.n1:
                    return self.L1
                return 0
            if delta == (0, 1):
                if state1[1]<=self.n2:
                    return self.L2
                return 0
            if delta == (-1, 0):
                if state1[1] <= self.n2+1:
                    return (1-self.r12)*self.mu1
                return 0
            if delta == (0, -1):
                if state1[0] <= self.n1+1:
                    return (1-self.r21)*self.mu2
                return 0
            if delta == (-1, 1):
                if state1[1] <= self.n2+1:
                    return self.r12*self.mu1
                return 0
            if delta == (1, -1):
                if state1[0] <= self.n1 + 1:
                    return self.r21*self.mu2
                return 0
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


    def find_median_time_to_absorption(self):
        """
        Finds the median time to absorption

            >>> Q = Network(1, 1, 6.0, 7.0, 0.1, 0.1, 3.0, 4.0)
            >>> Q.find_median_time_to_absorption()
        """
        points = [False for state in self.State_Space[:-1]]
        step = 1
        p = np.matrix(self.discrete_transition_matrix)

        probs = [0 for i in range(len(p)-1)]
        while not all(points):
            new_probs = [p[i,-1] for i in range(len(p)-1)]
            for j in range(len(probs)):
                if points[j] == False and new_probs[j] >= 0.5:
                    points[j] = [step-1, probs[j], step, new_probs[j]]
            probs = new_probs
            p = p*np.matrix(self.discrete_transition_matrix)
            step += 1

        steps2absorb = [self.interpolate(pnt) for pnt in points]
        time2absorb = [s*self.time_step for s in steps2absorb]
        self.median_steps_to_absorbtion = {str(self.State_Space[i]): steps2absorb[i] for i in range(len(steps2absorb))}
        self.median_time_to_absorbtion = {str(self.State_Space[i]): float(time2absorb[i]) for i in range(len(time2absorb))}

    def interpolate(self, pnt):
        """
        Takes in list of form [y1, x1, y2, x2] and interpolates to find y value that corresponds to x=0.5
        """
        y1 = pnt[0]
        y2 = pnt[2]
        x1 = pnt[1]
        x2 = pnt[3]
        m = (x2-x1)/(y2-y1)
        c = y1 - (x1*m)
        return (0.5*m) + c

    def write_results_to_file(self):
        """
        Takes the summary statistics and writes them into a .yml file
        """
        results_file = open('%stheoretical_results_mean_%s.yml' % (directory, str(r12)), 'w')
        results_file.write(yaml.dump(self.mean_time_to_absorbtion, default_flow_style=False))
        results_file.close()

        results_file = open('%stheoretical_results_median_%s.yml' % (directory, str(r12)), 'w')
        results_file.write(yaml.dump(self.median_time_to_absorbtion, default_flow_style=False))
        results_file.close()

if __name__ == '__main__':
    r12s = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    for r12 in r12s:
        Q = Network(n1, n2, mu1, mu2, r12, r21, L1, L2)
        Q.find_mean_time_to_absorbtion()
        Q.find_median_time_to_absorption()
        Q.write_results_to_file()
