from __future__ import division
import numpy as np
import yaml


# Input parameters here
n1, n2 = 3, 2
mu1, mu2 = 10.0, 8.0
r12, r21 = 0.25, 0.15
L1, L2 = 4.0, 5.0
directory = '/Users/geraintianpalmer/Documents/DetectingDeadlockInQingNetworkSimulation/data_for_graphs/vary_n2/'

class Network:
    """
    A class to hod the queueing network object
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
                else:
                    return 0
            if delta == (0, 1):
                if state1[1]<=self.n2:
                    return self.L2
                else:
                    return 0
            if delta == (-1, 0):
                if state1[1] == self.n2+2:
                    return 0
                else:
                    return (1-self.r12)*self.mu1
            if delta == (0, -1):
                if state1[0] == self.n1+2:
                    return 0
                else:
                    return (1-self.r21)*self.mu2
            if delta == (-1, 1):
                if state1[1] == self.n2+2:
                    return 0
                else:
                    return self.r12*self.mu1
            if delta == (1, -1):
                if state1[0] == self.n1 + 2:
                    return 0
                else:
                    return self.r21*self.mu2
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
                   [  0. ,   0. ,   0. ,   5.4,   0. ,   0. ,   6.3, -13.1,   0. ,
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
                      5.4,   0.6,   6.3, -12.9,   0.6],
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
                   [ 0.   ,  0.   ,  0.   ,  0.27 ,  0.   ,  0.   ,  0.315,  0.345,
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
                     0.   ,  0.27 ,  0.03 ,  0.315,  0.355,  0.03 ],
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
            {'(0, 1)': 5569.3309485569425, '(3, 0)': 5431.1549079579117, '(0, 0)': 5571.622376104965, '(1, 0)': 5568.0109461689945, '(3, 1)': 5232.8255050133139, '(1, 1)': 5562.8594988480072, '(2, 1)': 5547.3770517645999, '(1, 2)': 5552.4262511938905, '(2, 0)': 5561.4483015679771, '(0, 3)': 5462.6681009914637, '(1, 3)': 5218.4671485026602, '(0, 2)': 5565.8065378684014, '(2, 2)': 5518.8950250134531}
            >>> Q.mean_time_to_absorbtion
            {'(0, 1)': 278.46654742784716, '(3, 0)': 271.5577453978956, '(0, 0)': 278.58111880524825, '(1, 0)': 278.40054730844975, '(3, 1)': 261.64127525066573, '(1, 1)': 278.1429749424004, '(2, 1)': 277.36885258823, '(1, 2)': 277.62131255969456, '(2, 0)': 278.0724150783989, '(0, 3)': 273.1334050495732, '(1, 3)': 260.923357425133, '(0, 2)': 278.29032689342006, '(2, 2)': 275.94475125067265}
        """
        T = self.discrete_transition_matrix[:-1, :-1]
        S = np.linalg.inv(np.identity(len(T)) - T)
        steps2absorb = [sum([S[i,j] for j in range(len(S))]) for i in range(len(S))]
        time2absorb = [s*self.time_step for s in steps2absorb]
        self.mean_steps_to_absorbtion = {str(self.State_Space[i]): steps2absorb[i] for i in range(len(steps2absorb))}
        self.mean_time_to_absorbtion = {str(self.State_Space[i]): float(time2absorb[i]) for i in range(len(time2absorb))}

    def write_results_to_file(self):
		"""
		Takes the summary statistics and writes them into a .yml file
		"""
		results_file = open('%stheoretical_results_%s.yml' % (directory, str(n2)), 'w')
		results_file.write(yaml.dump(self.mean_time_to_absorbtion, default_flow_style=False))
		results_file.close()

if __name__ == '__main__':
    n2s = [0, 1, 2, 3, 4, 5, 6]
    for i in range(len(n2s)):
        n2 = n2s[i]
        Q = Network(n1, n2, mu1, mu2, r12, r21, L1, L2)
        Q.find_mean_time_to_absorbtion()
        Q.write_results_to_file()
