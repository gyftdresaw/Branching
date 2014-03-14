
# basic framework for simulating general time dependent galton watson processes
# key features:
#   - general waiting time distributions
#   - general branching distributions
#   - organize and extract simulated process data

import numpy as np

class GWEngine:
    # Galton Watson Engine is created with the arguments
    #   waiting_gen: function that returns waiting time from appropriate distribution
    #   child_gen: function that returns number of children from appropraite distribution
    #   gen_test: number of test RV's drawn from gens to approximate mean/var
    def __init__(self, waiting_gen, child_gen, gen_test = 1000):
        self.waiting_gen = waiting_gen
        self.child_gen = child_gen
        self.test_gens(gen_test)
        
    # sample means/vars for both waiting time and child distributions
    def test_gens(self, gen_test):
        waiting_sample = [self.waiting_gen() for i in xrange(gen_test)]
        child_sample = [self.child_gen() for i in xrange(gen_test)]
        self.waiting_mean = np.mean(waiting_sample)
        self.waiting_std = np.std(waiting_sample)
        self.child_mean = np.mean(child_sample)
        self.child_std = np.std(child_sample)

    # simulate single branching process trial
    #   t: simulation time (all branching stops after time t)
    #   num_std: number of buffer standard deviations for initial size of data matrix
    def simulate_branching(self, t, num_std = 3):
        # hold data in matrix format
        # we'll decide how large this matrix ought to be based on mean/var estimates
        # on average the growth should be something like child^(t/waiting)
        # number of time points ought to be (t/waiting)
        # we'll overestimate this using the num_std buffer
        child_over = self.child_mean + num_std * self.child_std
        waiting_under = self.waiting_mean - num_std * self.waiting_std
        max_branches = int(child_over ** (t / waiting_under))
        max_timepoints = 
        
GW = GWEngine(np.random.exponential,lambda: 2)


