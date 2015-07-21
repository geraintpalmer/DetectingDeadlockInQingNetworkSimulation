from __future__ import division
import numpy as np
import yaml
import random


# Input parameters here
n1, n2 = 3, 2
mu1, mu2 = 10.0, 8.0
r11, r12, r21, r22 = 0.1, 0.25, 0.15, 0.2
L1, L2 = 4.0, 5.0
directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/2NodeSimple/run_10000_itrs/vary_n2/'
class Network:
    """
    A class to hold the queueing network object
    """

    def __init__(self, n1, n2, mu1, mu2, r11, r12, r21, r22, L1, L2):
        """
        Initialises the Network object

            >>> Q = Network(0, 0, 10.0, 8.0, 0.1, 0.4, 0.1, 0.3, 6.0, 5.0)
            >>> Q.n1
            0
            >>> Q.n2
            0
            >>> Q.mu1
            10.0
            >>> Q.mu2
            8.0
            >>> Q.r11
            0.1
            >>> Q.r12
            0.4
            >>> Q.r21
            0.1
            >>> Q.r22
            0.3
            >>> Q.L1
            6.0
            >>> Q.L2
            5.0
            >>> Q.State_Space
            [((0, 0, 0), (0, 0, 0)), ((0, 0, 0), (1, 0, 0)), ((0, 0, 1), (1, 0, 0)), ((1, 0, 0), (0, 0, 0)), ((1, 0, 0), (0, 1, 0)), ((1, 0, 0), (1, 0, 0)), -1, -2, -3]
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
        self.State_Space = self.find_state_space()
        self.write_transition_matrix()
        self.discretise_transition_matrix()

    def find_state_space(self):
        """
        Returns the state space of the Queueing network.
        """
        all_states =  [((k1, b11, b12), (k2, b21, b22)) for k1 in range(self.n1+2) for k2 in range(self.n2+2) for b11 in range(2) for b12 in range(2) for b21 in range(2) for b22 in range(2)]
        the_states = []
        for state in all_states:
            if self.select_states(state):
                the_states.append(state)
        the_states.append(-1)
        the_states.append(-2)
        the_states.append(-3)
        return the_states

    def select_states(self, state):
        """
        Returns True is state in state space
        """
        if state[0][0]+state[0][1]+state[0][2] > self.n1+1:
            return False
        if state[1][0]+state[1][1]+state[1][2] > self.n2+1:
            return False
        if state[0][1]+state[0][2] > 1:
            return False
        if state[1][1]+state[1][2] > 1:
            return False
        if state[0][1]==1 and state[0][0]!=self.n1:
            return False
        if state[1][2]==1 and state[1][0]!=self.n2:
            return False
        if state[0][2]==1 and state[1][0]+state[1][1]+state[1][2]!=self.n2+1:
            return False
        if state[1][1]==1 and state[0][0]+state[0][1]+state[0][2]!=self.n1+1:
            return False
        if state[0] == (self.n1, 1, 0):
            return False
        if state[1] == (self.n2, 0, 1):
            return False
        if state == ((self.n1, 0, 1), (self.n2, 1, 0)):
            return False
        return True


    def find_transition_rates(self, state1, state2):
        """
        Finds the transition rates for given state transition

            Testing all possible types of transitions
            >>> Q = Network(3, 2, 10.0, 8.0, 0.1, 0.4, 0.3, 0.2, 6.0, 5.0)
            >>> Q.find_transition_rates(((0, 0, 0), (0, 0, 0)), -1)
            0
            >>> Q.find_transition_rates(((0, 0, 0), (0, 0, 0)), -2)
            0
            >>> Q.find_transition_rates(((0, 0, 0), (0, 0, 0)), -3)
            0
            >>> Q.find_transition_rates(((3, 1, 0), (1, 0, 0)), -1)
            1.0
            >>> Q.find_transition_rates(((1, 0, 0), (2, 0, 1)), -2)
            1.6
            >>> Q.find_transition_rates(((3, 0, 1), (3, 0, 0)), -3)
            2.4
            >>> Q.find_transition_rates(((4, 0, 0), (2, 1, 0)), -3)
            4.0
            >>> Q.find_transition_rates(((0, 0, 0), (1, 0, 0)), ((1, 0, 0), (1, 0, 0)))
            6.0
        """
        if state1 == -1 or state1 == -2 or state1 == -3:
            return 0
        elif state2 == -1:
            if state1[0] == (self.n1, 1, 0):
                return self.r11*self.mu1
            return 0
        elif state2 == -2:
            if state1[1] == (self.n2, 0, 1):
                return self.r22*self.mu2
            return 0
        elif state2 == -3:
            if state1 == ((self.n1, 0, 1), (self.n2+1, 0, 0)):
                return self.r21*self.mu2
            if state1 == ((self.n1+1, 0, 0), (self.n2, 1, 0)):
                return self.r12*self.mu1
            return 0
        else:
            delta = ((state2[0][0]-state1[0][0], state2[0][1]-state1[0][1], state2[0][2]-state1[0][2]), (state2[1][0]-state1[1][0], state2[1][1]-state1[1][1], state2[1][2]-state1[1][2]))
            # Arrival at 1
            if delta == ((1, 0, 0), (0, 0, 0)):
                if state1[0][0]+state1[0][1]+state1[0][2]<self.n1+1:
                    return self.L1
            # Arrival at 2
            if delta == ((0, 0, 0), (1, 0, 0)):
                if state1[1][0]+state1[1][1]+state1[1][2]<self.n2+1:
                    return self.L1

            # Service and transition 1 - 2
            if delta == ((-1, 0, 0), (1, 0, 0)):
                if state1[0][1] == 0 and state1[0][2] == 0:
                    if state1[1][0]+state1[1][1]+state1[1][2]<self.n2+1:
                        return self.r12*self.mu1
            # Service and transition 1 - 2 with an unblock
            if delta == ((0, 0, 0), (1, -1, 0)):
                if state1[0][1] == 0 and state1[0][2] == 0 and state1[0][0]==self.n1+1:
                    if state1[1][0]+state1[1][1]+state1[1][2]<self.n2+1:
                        return self.r12*self.mu1
            # Service want to trans 1-2 but block
            if delta == ((-1, 0, 1), (0, 0, 0)):
                if state1[0][1] == 0 and state1[0][2] == 0:
                    if state1[1][0] + state1[1][1] + state1[1][2] == self.n2+1:
                        return self.r12*self.mu1
            # Service at N1 and exit
            if delta == ((-1, 0, 0), (0, 0, 0)):
                if state1[0][1] == 0 and state1[0][2] == 0:
                    return (1-self.r11-self.r12)*self.mu1
            # Service ar N1 and exit, and unblock
            if delta == ((0, 0, 0), (0, -1, 0)):
                if state1[0][1] == 0 and state1[0][2] == 0 and state1[0][0]==self.n1+1:
                    return (1-self.r11-self.r12)*self.mu1

            # Service and transition 2 - 1
            if delta == ((1, 0, 0), (-1, 0, 0)):
                if state1[1][1] == 0 and state1[1][2] == 0:
                    if state1[0][0]+state1[0][1]+state1[0][2]<self.n1+1:
                        return self.r21*self.mu2
            # Service and transition 2 - 1 with an unblock
            if delta == ((1, -1, 0), (0, 0, 0)):
                if state1[1][1] == 0 and state1[1][2] == 0 and state1[1][0]==self.n1+1:
                    if state1[0][0]+state1[0][1]+state1[0][2]<self.n1+1:
                        return self.r21*self.mu2
            # Service want to trans 2-1 but block
            if delta == ((0, 0, 0), (-1, 0, 1)):
                if state1[1][1] == 0 and state1[1][2] == 0:
                    if state1[0][0] + state1[0][1] + state1[0][2] == self.n1+1:
                        return self.r21*self.mu2
            # Service at N1 and exit
            if delta == ((0, 0, 0), (-1, 0, 0)):
                if state1[1][1] == 0 and state1[1][2] == 0:
                    return (1-self.r21-self.r22)*self.mu2
            # Service ar N1 and exit, and unblock
            if delta == ((0, -1, 0), (0, 0, 0)):
                if state1[1][1] == 0 and state1[1][2] == 0 and state1[1][0]==self.n2+1:
                    return (1-self.r21-self.r22)*self.mu2
        return 0

    def write_transition_matrix(self):
        """
        Writes the transition matrix for the markov chain

            >>> Q = Network(0, 0, 6.0, 7.0, 0.1, 0.1, 0.1, 0.1, 3.0, 4.0)
            >>> Q.write_transition_matrix()
            >>> Q.transition_matrix
            array([[ -6. ,   3. ,   0. ,   3. ,   0. ,   0. ,   0. ,   0. ,   0. ],
                   [  5.6,  -9.3,   0. ,   0.7,   0. ,   3. ,   0. ,   0. ,   0. ],
                   [  0. ,   0. ,  -0.7,   0. ,   0. ,   0. ,   0. ,   0. ,   0.7],
                   [  4.8,   0.6,   0. ,  -8.4,   0. ,   3. ,   0. ,   0. ,   0. ],
                   [  0. ,   0. ,   0. ,   4.8,  -5.4,   0. ,   0. ,   0. ,   0.6],
                   [  0. ,   4.8,   0.6,   5.6,   0. , -11. ,   0. ,   0. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,  -0. ,   0. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,  -0. ,   0. ],
                   [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,  -0. ]])
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
            array([[ 0.45454545,  0.27272727,  0.        ,  0.27272727,  0.        ,
                     0.        ,  0.        ,  0.        ,  0.        ],
                   [ 0.50909091,  0.15454545,  0.        ,  0.06363636,  0.        ,
                     0.27272727,  0.        ,  0.        ,  0.        ],
                   [ 0.        ,  0.        ,  0.93636364,  0.        ,  0.        ,
                     0.        ,  0.        ,  0.        ,  0.06363636],
                   [ 0.43636364,  0.05454545,  0.        ,  0.23636364,  0.        ,
                     0.27272727,  0.        ,  0.        ,  0.        ],
                   [ 0.        ,  0.        ,  0.        ,  0.43636364,  0.50909091,
                     0.        ,  0.        ,  0.        ,  0.05454545],
                   [ 0.        ,  0.43636364,  0.05454545,  0.50909091,  0.        ,
                     0.        ,  0.        ,  0.        ,  0.        ],
                   [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                     0.        ,  1.        ,  0.        ,  0.        ],
                   [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                     0.        ,  0.        ,  1.        ,  0.        ],
                   [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                     0.        ,  0.        ,  0.        ,  1.        ]])
        """
        self.time_step = 1 / max([abs(self.transition_matrix[i][i]) for i in range(len(self.transition_matrix))])
        self.discrete_transition_matrix = self.transition_matrix*self.time_step + np.identity(len(self.transition_matrix))

    def find_mean_time_to_absorbtion(self):
        """
        Finds the mean time to absorbtion

            >>> Q = Network(1, 1, 6.0, 7.0, 0.1, 0.3, 0.3, 0.1, 3.0, 4.0)
            >>> Q.find_mean_time_to_absorbtion()
            >>> Q.mean_steps_to_absorbtion
            {'((0, 0, 0), (1, 0, 0))': 106.88633063796904, '((1, 0, 0), (2, 0, 0))': 86.458216192148484, '((2, 0, 0), (1, 0, 0))': 96.19655504160005, '((1, 0, 1), (2, 0, 0))': 8.4285714285714395, '((1, 0, 0), (0, 0, 0))': 107.47216921155163, '((2, 0, 0), (2, 0, 0))': 77.931931045612842, '((0, 0, 1), (2, 0, 0))': 14.328571428571442, '((1, 0, 0), (1, 0, 0))': 101.85806091072483, '((2, 0, 0), (0, 1, 0))': 91.516388592125281, '((2, 0, 0), (1, 1, 0))': 67.408814472177809, '((0, 0, 0), (2, 0, 0))': 101.06442636219762, '((0, 0, 0), (0, 0, 0))': 110.12924992476033, '((2, 0, 0), (0, 0, 0))': 104.3492838006775}
            >>> Q.mean_time_to_absorbtion
            {'((0, 0, 0), (1, 0, 0))': 6.038775742258137, '((1, 0, 0), (2, 0, 0))': 4.884644982607258, '((2, 0, 0), (1, 0, 0))': 5.434833618169494, '((1, 0, 1), (2, 0, 0))': 0.4761904761904767, '((1, 0, 0), (0, 0, 0))': 6.071873966754328, '((2, 0, 0), (2, 0, 0))': 4.402933957379256, '((0, 0, 1), (2, 0, 0))': 0.8095238095238102, '((1, 0, 0), (1, 0, 0))': 5.754692706820611, '((2, 0, 0), (0, 1, 0))': 5.1704174345833485, '((2, 0, 0), (1, 1, 0))': 3.808407597298181, '((0, 0, 0), (2, 0, 0))': 5.709854596734329, '((0, 0, 0), (0, 0, 0))': 6.2219915211728996, '((2, 0, 0), (0, 0, 0))': 5.89543976275014}
        """
        T = self.discrete_transition_matrix[:-3, :-3]        
        S = np.linalg.inv(np.identity(len(T)) - T)
        steps2absorb = [sum([S[i,j] for j in range(len(S))]) for i in range(len(S))]
        time2absorb = [s*self.time_step for s in steps2absorb]
        self.mean_steps_to_absorbtion = {str(self.State_Space[i]): steps2absorb[i] for i in range(len(steps2absorb))}
        self.mean_time_to_absorbtion = {str(self.State_Space[i]): float(time2absorb[i]) for i in range(len(time2absorb))}


    def find_median_time_to_absorption(self):
        """
        Finds the median time to absorption
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


    def write_results_to_file(self, param):
        """
        Takes the summary statistics and writes them into a .yml file
        """
        results_file = open('%stheoretical_results_%s.yml' % (directory, str(param)), 'w')
        results_file.write(yaml.dump(self.mean_time_to_absorbtion, default_flow_style=False))
        results_file.close()

        # results_file = open('%stheoretical_results_median_%s.yml' % (directory, str(param)), 'w')
        # results_file.write(yaml.dump(self.median_time_to_absorbtion, default_flow_style=False))
        # results_file.close()


if __name__ == '__main__':
    # n2s = [0, 1, 2, 3, 4, 5, 6]
    # for n2 in n2s:
    #     Q = Network(n1, n2, mu1, mu2, r11, r12, r21, r22, L1, L2)
    #     Q.find_mean_time_to_absorbtion()
    #     # Q.find_median_time_to_absorption()
    #     Q.write_results_to_file(n2)
    #     print "Now starting n2 = " + str(n2)
    #     # Q.simulate_markov_chain(1000)
    #     # Q.write_simulation_results_to_file(L1)
    Q = Network(n1, n2, mu1, mu2, r11, r12, r21, r22, L1, L2)