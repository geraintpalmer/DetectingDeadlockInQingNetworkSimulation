from __future__ import division
import numpy as np


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
                return (1-self.r12)*self.mu1*min(state1[0], 1)
            if delta == (0, -1):
                return (1-self.r21)*self.mu2*min(state1[1], 1)
            if delta == (-1, 1):
                if state1[1] == self.n2+2:
                    return 0
                else:
                    return self.r12*self.mu1*min(state1[0], 1)
            if delta == (1, -1):
                if state1[0] == self.n1 + 2:
                    return 0
                else:
                    return self.r21*self.mu2*min(state1[1], 1)
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
            {(0, 1): 5569.3309485569425, (1, 2): 5552.4262511938905, (0, 0): 5571.622376104965, (3, 0): 5431.1549079579117, (3, 1): 5232.8255050133139, (2, 1): 5547.3770517645999, (0, 2): 5565.8065378684014, (2, 0): 5561.4483015679771, (1, 3): 5218.4671485026602, (2, 2): 5518.8950250134531, (1, 0): 5568.0109461689945, (0, 3): 5462.6681009914637, (1, 1): 5562.8594988480072}
            >>> Q.mean_time_to_absorbtion
            {(0, 1): 278.46654742784716, (1, 2): 277.62131255969456, (0, 0): 278.58111880524825, (3, 0): 271.55774539789559, (3, 1): 261.64127525066573, (2, 1): 277.36885258823003, (0, 2): 278.29032689342006, (2, 0): 278.07241507839888, (1, 3): 260.92335742513302, (2, 2): 275.94475125067265, (1, 0): 278.40054730844975, (0, 3): 273.13340504957318, (1, 1): 278.14297494240037}
        """
        T = self.discrete_transition_matrix[:-1, :-1]
        S = np.linalg.inv(np.identity(len(T)) - T)
        steps2absorb = [sum(S[i]) for i in range(len(S))]
        time2absorb = [s*self.time_step for s in steps2absorb]
        self.mean_steps_to_absorbtion = {self.State_Space[i]: steps2absorb[i] for i in range(len(steps2absorb))}
        self.mean_time_to_absorbtion = {self.State_Space[i]: time2absorb[i] for i in range(len(time2absorb))}
