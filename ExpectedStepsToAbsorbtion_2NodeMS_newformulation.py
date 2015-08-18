from __future__ import division
import numpy as np
import yaml
import random


# Input parameters here
# n1, n2 = 3, 2
# c1, c2 = 2, 2
# mu1, mu2 = 8.0, 6.0
# r12, r21 = 0.5, 0.5
# L1, L2 = 5.0, 6.0
# directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/2NodeSimple_MultiServer/run_10000_itrs/vary_c1/'
# n1, n2 = 3, 2
# c1, c2 = 1, 1
# mu1, mu2 = 10.0, 8.0
# r12, r21 = 0.25, 0.15
# L1, L2 = 4.0, 5.0
# directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/test/vary_mu1/'
# n1, n2 = 2, 3
# c1, c2 = 4, 2
# mu1, mu2 = 6.0, 7.0
# r12, r21 = 0.7, 0.6
# L1, L2 = 7.5, 6.0
# directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/2NodeSimple_MultiServer/run_10000_itrs/vary_L1/'
n1, n2 = 0, 0
c1, c2 = 4, 2
mu1, mu2 = 6.0, 7.0
r12, r21 = 0.7, 0.6
L1, L2 = 7.5, 7.0
directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/2NodeSimple_MultiServer/Results_test/'



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
            [((0, 0), (0, 0)), ((0, 0), (1, 0)), ((0, 0), (2, 0)), ((0, 0), (3, 0)), ((0, 0), (4, 0)), ((1, 0), (0, 0)), ((1, 0), (1, 0)), ((1, 0), (2, 0)), ((1, 0), (3, 0)), ((1, 0), (4, 0)), ((1, 0), (4, 1)), ((2, 0), (0, 0)), ((2, 0), (1, 0)), ((2, 0), (2, 0)), ((2, 0), (3, 0)), ((2, 0), (4, 0)), ((2, 0), (4, 1)), ((2, 0), (4, 2)), ((3, 0), (0, 0)), ((3, 0), (1, 0)), ((3, 0), (2, 0)), ((3, 0), (3, 0)), ((3, 0), (4, 0)), ((3, 0), (4, 1)), ((3, 0), (4, 2)), ((4, 0), (0, 0)), ((4, 0), (1, 0)), ((4, 0), (2, 0)), ((4, 0), (3, 0)), ((4, 0), (4, 0)), ((4, 0), (4, 1)), ((4, 0), (4, 2)), ((5, 0), (0, 0)), ((5, 0), (1, 0)), ((5, 0), (2, 0)), ((5, 0), (3, 0)), ((5, 0), (4, 0)), ((5, 0), (4, 1)), ((5, 0), (4, 2)), ((5, 1), (1, 0)), ((5, 1), (2, 0)), ((5, 1), (3, 0)), ((5, 1), (4, 0)), ((5, 1), (4, 1)), ((5, 1), (4, 2)), ((5, 2), (2, 0)), ((5, 2), (3, 0)), ((5, 2), (4, 0)), ((5, 2), (4, 1)), ((5, 2), (4, 2))]
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
        self.State_Space = [((k1, b1), (k2, b2)) for k1 in range(self.n1+self.c1+1) for k2 in range(self.n2+self.c2+1) for b1 in range(self.c2+1) for b2 in range(self.c1+1) if b1<=k2 if b2<=k1]
        self.prune_state_space()
        self.write_transition_matrix()
        self.discretise_transition_matrix()

    def prune_state_space(self):
        """
        Prunes the state space so that only reachable states are left
        """
        infeasible_states1 = [state for state in self.State_Space if state[0][1] > 0 and state[0][0] < self.n1 + self.c1]
        infeasible_states2 = [state for state in self.State_Space if state[1][1] > 0 and state[1][0] < self.n2 + self.c2]
        feasible_states = [state for state in self.State_Space if state not in infeasible_states1 if state not in infeasible_states2]
        self.State_Space = sorted(feasible_states, key=lambda state: (state[0][0]+state[0][1],state[1][0]+state[1][1]))


    def find_transition_rates(self, state1, state2):
        """
        Finds the transition rates for given state transition

            >>> Q = Network(0, 0, 6.0, 5.0, 0.4, 0.6, )
        """
        k11, b11, k21, b21 = state1[0][0], state1[0][1], state1[1][0], state1[1][1]
        k12, b12, k22, b22 = state2[0][0], state2[0][1], state2[1][0], state2[1][1]
        delta = ((k12-k11, b12-b11), (k22-k21, b22-b21))
        S1 = min(k11, self.c1) - b21
        S2 = min(k21, self.c2) - b11

        if b11 == 0 and k11 < self.n1 + self.c1:
            if b21 == 0 and k21 < self.n2 + self.c2:
                if delta == ((1, 0), (0, 0)):
                    return self.L1
                if delta == ((0, 0), (1, 0)):
                    return self.L2
                if delta == ((-1, 0), (0, 0)):
                    return (1-self.r12)*S1*self.mu1
                if delta == ((0, 0), (-1, 0)):
                    return (1-self.r21)*S2*self.mu2
                if delta == ((-1, 0), (1, 0)):
                    return self.r12*S1*self.mu1
                if delta == ((1, 0), (-1, 0)):
                    return self.r21*S2*self.mu2
            if b21 == 0 and k21 == self.n2 + self.c2:
                if delta == ((1, 0), (0, 0)):
                    return self.L1
                if delta == ((-1, 0), (0, 0)):
                    return (1-self.r12)*S1*self.mu1
                if delta == ((0, 0), (-1, 0)):
                    return (1-self.r21)*S2*self.mu2
                if delta == ((0, 0), (0, 1)):
                    return self.r12*S1*self.mu1
                if delta == ((1, 0), (-1, 0)):
                    return self.r21*S2*self.mu2
            if b21 > 0:
                if delta == ((1, 0), (0, 0)):
                    return self.L1
                if delta == ((-1, 0), (0, 0)):
                    return (1-self.r12)*S1*self.mu1
                if delta == ((-1, 0), (0, -1)):
                    return (1-self.r21)*S2*self.mu2
                if delta == ((0, 0), (0, 1)):
                    return self.r12*S1*self.mu1
                if delta == ((0, 0), (0, -1)):
                    return self.r21*S2*self.mu2
        if b11 == 0 and k11 == self.n1 + self.c1:
            if b21 == 0 and k21 < self.n2 + self.c2:
                if delta == ((0, 0), (1, 0)):
                    return self.L2
                if delta == ((-1, 0), (0, 0)):
                    return (1-self.r12)*S1*self.mu1
                if delta == ((0, 0), (-1, 0)):
                    return (1-self.r21)*S2*self.mu2
                if delta == ((-1, 0), (1, 0)):
                    return self.r12*S1*self.mu1
                if delta == ((0, 1), (0, 0)):
                    return self.r21*S2*self.mu2
            if b21 == 0 and k21 == self.n2 + self.c2:
                if delta == ((-1, 0), (0, 0)):
                    return (1-self.r12)*S1*self.mu1
                if delta == ((0, 0), (-1, 0)):
                    return (1-self.r21)*S2*self.mu2
                if delta == ((0, 0), (0, 1)):
                    return self.r12*S1*self.mu1
                if delta == ((0, 1), (0, 0)):
                    return self.r21*S2*self.mu2
            if b21 > 0:
                if delta == ((-1, 0), (0, 0)):
                    return (1-self.r12)*S1*self.mu1
                if delta == ((-1, 0), (0, -1)):
                    return (1-self.r21)*S2*self.mu2
                if delta == ((0, 0), (0, 1)):
                    return self.r12*S1*self.mu1
                if delta == ((0, 1), (0, 0)):
                    return self.r21*S2*self.mu2
        if b11 > 0:
            if b21 == 0 and k21 < self.n2 + self.c2:
                if delta == ((0, 0), (1, 0)):
                    return self.L2
                if delta == ((0, -1), (-1, 0)):
                    return (1-self.r12)*S1*self.mu1
                if delta == ((0, 0), (-1, 0)):
                    return (1-self.r21)*S2*self.mu2
                if delta == ((0, -1), (0, 0)):
                    return self.r12*S1*self.mu1
                if delta == ((0, 1), (0, 0)):
                    return self.r21*S2*self.mu2
            if b21 == 0 and k21 == self.n2 + self.c2:
                if delta == ((0, -1), (-1, 0)):
                    return (1-self.r12)*S1*self.mu1
                if delta == ((0, 0), (-1, 0)):
                    return (1-self.r21)*S2*self.mu2
                if delta == ((0, 0), (0, 1)):
                    return self.r12*S1*self.mu1
                if delta == ((0, 1), (0, 0)):
                    return self.r21*S2*self.mu2
            if b21 > 0:
                w1 = -1 if b11 <= b21 else 0
                w2 = -1 if b21 > b11 else 0
                y1 = -1 if b11 > b21 else 0
                y2 = -1 if b21 <= b11 else 0
                x1 = -min(b11, b21+1)
                x2 = -min(b11, b21)
                z1 = -min(b11, b21)
                z2 = -min(b11+1, b21)
                if delta == ((0, 0), (0, 1)):
                    return self.r12*S1*self.mu1
                if delta == ((0, 1), (0, 0)):
                    return self.r21*S2*self.mu2
                if delta == ((w1, x1), (y1, z1)):
                    return (1-self.r12)*S1*self.mu1
                if delta == ((w2, x2), (y2, z2)):
                    return (1-self.r21)*S2*self.mu2
        return 0

    def write_transition_matrix(self):
        """
        Writes the transition matrix for the markov chain
        """
        self.transition_matrix = [[self.find_transition_rates(s1, s2) for s2 in self.State_Space] for s1 in self.State_Space]
        for i in range(len(self.transition_matrix)):
            a = sum(self.transition_matrix[i])
            self.transition_matrix[i][i] = -a
            self.transition_matrix = np.array(self.transition_matrix)

    def discretise_transition_matrix(self):
        """
        Disctetises the transition matrix
        """
        self.time_step = 1 / max([abs(self.transition_matrix[i][i]) for i in range(len(self.transition_matrix))])
        self.discrete_transition_matrix = self.transition_matrix*self.time_step + np.identity(len(self.transition_matrix))

    def find_mean_time_to_absorbtion(self):
        """
        Finds the mean time to absorbtion
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
        results_file = open('%stransmat_nf4.csv' % directory, 'w')
        np.savetxt(results_file, Q.transition_matrix, delimiter=",")


if __name__ == '__main__':
    # mu1s = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0]
    # c1s = [1, 2, 3, 4, 5]
    L1s = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
    for L1 in L1s:
        Q = Network(n1, n2, mu1, mu2, r12, r21, L1, L2, c1, c2)
        Q.find_mean_time_to_absorbtion()
        Q.write_results_to_file(L1)
        print "Now starting L1 = " + str(L1)
    # Q = Network(n1, n2, mu1, mu2, r12, r21, L1, L2, c1, c2)
    # print Q.State_Space
    # Q.write_matrix_to_file()
    # # # print Q.State_Space
    # # for state in Q.State_Space:
    # #     print (state[0][0]+state[0][1], state[1][0]+state[1][1])