from __future__ import division
import numpy as np
import yaml
import random

class Network:
    """
    A class to hold the queueing network object
    """

    def __init__(self, n1, n2, mu1, mu2, r11, r12, r21, r22, L1, L2):
        """
        Initialises the Network object

            >>> Q = Network(3, 2, 10.0, 8.0, 0.1, 0.4, 0.3, 0.1, 6.0, 5.0)
            >>> Q.n1
            3
            >>> Q.n2
            2
            >>> Q.mu1
            10.0
            >>> Q.mu2
            8.0
            >>> Q.r11
            0.1
            >>> Q.r12
            0.4
            >>> Q.r21
            0.3
            >>> Q.r22
            0.1
            >>> Q.L1
            6.0
            >>> Q.L2
            5.0
            >>> Q.State_Space
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (4, 0), (4, 1), (4, 2), (4, 3), (5, 0), (5, 1), (5, 2), -1, -2, -3]
        """
        self.n1 = n1
        self.n2 = n2
        self.mu1 = mu1
        self.mu2 = mu2
        self.r11 = r11
        self.r12 = r12
        self.r21 = r21
        self.r22 = r22
        self.L1 = L1
        self.L2 = L2
        self.State_Space = [(i, j) for i in range(self.n1+3) for j in range(self.n2+3) if i+j<=self.n1+self.n2+2] + [-1, -2, -3]
        self.write_transition_matrix()
        self.discretise_transition_matrix()

    def find_transition_rates(self, state1, state2):
        """
        Finds the transition rates for given state transition

            Testing all possible types of transitions
            >>> Q = Network(3, 2, 10.0, 8.0, 0.1, 0.4, 0.3, 0.1, 6.0, 5.0)
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
            4.8
            >>> Q.find_transition_rates((2, 2), (1, 3))
            4.0
            >>> Q.find_transition_rates((2, 2), (3, 1))
            2.4
            >>> Q.find_transition_rates((2, 4), (1, 5))
            0
            >>> Q.find_transition_rates((5, 2), (6, 1))
            0
        """
        if state1 in [-1, -2, -3]:
            return 0
        if state2 == -3:
            if state1[0] == self.n1 and state1[1] == self.n2 + 2:
                return self.r21 * self.mu2
            if state1[0] == self.n1 + 2 and state1[1] == self.n2:
                return self.r12 * self.mu1
            else:
                return 0
        elif state2 == -1:
            if state1[0] >= self.n1+1 and state1[1] < self.n2+2:
                return self.r11*self.mu1
            else:
                return 0
        elif state2 == -2:
            if state1[1] >= self.n2+1 and state1[0] < self.n1+2:
                return self.r22*self.mu2
            else:
                return 0
        else:
            delta = (state2[0] - state1[0], state2[1] - state1[1])
            if delta == (1, 0):
                if state1[0] < self.n1 + 1:
                    return self.L1
                return 0
            if delta == (0, 1):
                if state1[1] < self.n2 + 1:
                    return self.L2
                return 0
            if delta == (-1, 0):
                if state1[1] < self.n2 + 2:
                    return (1 - self.r12 - self.r11) * self.mu1
                return 0
            if delta == (0, -1):
                if state1[0] < self.n1 + 2:
                    return (1 - self.r21 - self.r22) * self.mu2
                return 0
            if delta == (-1, 1):
                if state1[1] < self.n2 + 2 and (state1[0], state1[1]) != (self.n1+2, self.n2):
                # if state1[1] < self.n2 + 2:
                    return self.r12 * self.mu1
                return 0
            if delta == (1, -1):
                if state1[0] < self.n1 + 2 and (state1[0], state1[1]) != (self.n1, self.n2+2):
                # if state1[0] < self.n1 + 2:
                    return self.r21 * self.mu2
                return 0
            return 0

    def write_transition_matrix(self):
        """
        Writes the transition matrix for the markov chain

            >>> Q = Network(1, 1, 6.0, 7.0, 0.1, 0.1, 0.1, 0.1, 3.0, 4.0)
            >>> Q.write_transition_matrix()
            >>> Q.transition_matrix
            array([[ -7. ,   4. ,   0. ,   0. ,   3. ,   0. ,   0. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ],
                   [  5.6, -13.3,   4. ,   0. ,   0.7,   3. ,   0. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ],
                   [  0. ,   5.6, -10. ,   0. ,   0. ,   0.7,   3. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0.7],
                   [  0. ,   0. ,   5.6, -10. ,   0. ,   0. ,   0.7,   3. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0.7],
                   [  4.8,   0.6,   0. ,   0. , -12.4,   4. ,   0. ,   0. ,   3. ,
                      0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ],
                   [  0. ,   4.8,   0.6,   0. ,   5.6, -18.7,   4. ,   0. ,   0.7,
                      3. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ],
                   [  0. ,   0. ,   4.8,   0.6,   0. ,   5.6, -15.4,   0. ,   0. ,
                      0.7,   3. ,   0. ,   0. ,   0. ,   0. ,   0.7],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   5.6,  -7. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0.7,   0. ,   0.7],
                   [  0. ,   0. ,   0. ,   0. ,   4.8,   0.6,   0. ,   0. , -10. ,
                      4. ,   0. ,   0. ,   0. ,   0. ,   0.6,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   4.8,   0.6,   0. ,   5.6,
                    -16.3,   4. ,   0.7,   0. ,   0. ,   0.6,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   4.8,   0.6,   0. ,
                      5.6, -13. ,   0. ,   0.7,   0. ,   0.6,   0.7],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   4.8,
                      0.6,   0. , -10. ,   4. ,   0. ,   0.6,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,
                      4.8,   0. ,   0. ,  -6. ,   0.6,   0.6,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,  -0. ,   0. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ,  -0. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,
                      0. ,   0. ,   0. ,   0. ,   0. ,   0. ,  -0. ]])
        """
        self.transition_matrix = [[self.find_transition_rates(s1, s2) for s2 in self.State_Space] for s1 in self.State_Space]
        for i in range(len(self.transition_matrix)):
            a = sum(self.transition_matrix[i])
            self.transition_matrix[i][i] = -a
            self.transition_matrix = np.array(self.transition_matrix)

    def discretise_transition_matrix(self):
        """
        Disctetises the transition matrix

            >>> Q = Network(0, 0, 6.0, 7.0, 0.1, 0.1, 0.1, 0.1, 3.0, 4.0)
            >>> Q.write_transition_matrix()
            >>> Q.discretise_transition_matrix()
            >>> Q.discrete_transition_matrix
            array([[  4.61538462e-01,   3.07692308e-01,   0.00000000e+00,
                      2.30769231e-01,   0.00000000e+00,   0.00000000e+00,
                      0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
                   [  4.30769231e-01,   2.30769231e-01,   0.00000000e+00,
                      5.38461538e-02,   2.30769231e-01,   0.00000000e+00,
                      0.00000000e+00,   0.00000000e+00,   5.38461538e-02],
                   [  0.00000000e+00,   4.30769231e-01,   4.61538462e-01,
                      0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      5.38461538e-02,   0.00000000e+00,   5.38461538e-02],
                   [  3.69230769e-01,   4.61538462e-02,   0.00000000e+00,
                      2.30769231e-01,   3.07692308e-01,   0.00000000e+00,
                      0.00000000e+00,   4.61538462e-02,   0.00000000e+00],
                   [  0.00000000e+00,   3.69230769e-01,   4.61538462e-02,
                      4.30769231e-01,   1.11022302e-16,   5.38461538e-02,
                      0.00000000e+00,   4.61538462e-02,   5.38461538e-02],
                   [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      3.69230769e-01,   0.00000000e+00,   5.38461538e-01,
                      4.61538462e-02,   4.61538462e-02,   0.00000000e+00],
                   [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      1.00000000e+00,   0.00000000e+00,   0.00000000e+00],
                   [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      0.00000000e+00,   1.00000000e+00,   0.00000000e+00],
                   [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                      0.00000000e+00,   0.00000000e+00,   1.00000000e+00]])
        """
        self.time_step = 1 / max([abs(self.transition_matrix[i][i]) for i in range(len(self.transition_matrix))])
        self.discrete_transition_matrix = self.transition_matrix*self.time_step + np.identity(len(self.transition_matrix))

    def find_mean_time_to_absorbtion(self):
        """
        Finds the mean time to absorbtion

            >>> Q = Network(1, 1, 6.0, 7.0, 0.1, 0.1, 0.1, 0.1, 3.0, 4.0)
            >>> Q.find_mean_time_to_absorbtion()
            >>> Q.mean_steps_to_absorbtion
            {'(0, 1)': 74.486012826212615, '(3, 0)': 60.484428199084014, '(0, 0)': 77.176019347127678, '(1, 0)': 74.529361375014403, '(3, 1)': 55.195821120657477, '(1, 1)': 71.597265173978755, '(2, 1)': 65.098943067488506, '(1, 2)': 65.366675650696777, '(2, 0)': 67.979506597441031, '(0, 3)': 61.129325949581869, '(1, 3)': 54.964769091985993, '(0, 2)': 68.203978440066606, '(2, 2)': 59.125389272379401}
            >>> Q.mean_time_to_absorbtion
            {'(0, 1)': 3.983209242043455, '(3, 0)': 3.2344614010205355, '(0, 0)': 4.12705985813517, '(1, 0)': 3.9855273462574545, '(3, 1)': 2.951648188270453, '(1, 1)': 3.8287307579667784, '(2, 1)': 3.4812269020047326, '(1, 2)': 3.495544152443678, '(2, 0)': 3.6352677324834777, '(0, 3)': 3.268947911742346, '(1, 3)': 2.9392924648120853, '(0, 2)': 3.6472715743351127, '(2, 2)': 3.1617855225871336}
        """
        T = self.discrete_transition_matrix[:-3, :-3]
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

    def find_absorpion_probabilities(self):
        """
        Finds the absorbtion probabilities of the queueing network
        """
        T = self.discrete_transition_matrix[:-3, :-3]
        S = np.linalg.inv(np.identity(len(T)) - T)
        B = self.discrete_transition_matrix[:-3,-3:]
        A = np.matrix(S)*np.matrix(B)
        self.absorbtion_probabilities = {str(self.State_Space[i]): [A[0,j] for j in range(3)] for i in range(len(A))}

    def write_absorb_results_to_file(self, param):
        """
        Takes the summary statistics and writes them into a .yml file
        """
        results_file = open('%sabsorbtion_probabilities_%s.yml' % (directory, str(param)), 'w')
        results_file.write(yaml.dump(self.absorbtion_probabilities, default_flow_style=False))
        results_file.close()


if __name__ == '__main__':
    n1, n2 = 3, 2
    mu1, mu2 = 10.0, 8.0
    r11, r12, r21, r22 = 0.1, 0.25, 0.15, 0.1
    L1, L2 = 4.0, 5.0
    directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/2NodeFeedback/run_10000_itrs/vary_L2/'

    L2s = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0]
    # n1s = [0, 1, 2, 3, 4, 5, 6]
    # r22s = [0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0]
    for L2 in L2s:
        Q = Network(n1, n2, mu1, mu2, r11, r12, r21, r22, L1, L2)
        # Q.find_mean_time_to_absorbtion()
        Q.find_absorpion_probabilities()
        # Q.write_results_to_file(L1)
        Q.write_absorb_results_to_file(L2)
        print "Now starting L2 = " + str(L2)