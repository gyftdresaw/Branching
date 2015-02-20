
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
        
    # get samples of both waiting time and child distributions
    def test_gens(self, gen_test):
        self.waiting_sample = sorted([self.waiting_gen() for i in xrange(gen_test)])
        self.child_sample = sorted([self.child_gen() for i in xrange(gen_test)])
        self.test_samples = gen_test

    # simulate single branching process trial
    #   t: simulation time (all branching stops after time t)
    #   prct_thresh: percentile from which process bounds are based off
    #
    # we're going to assume we start with single cell at time 0.0
    def simulate_branching(self, t, prct_thresh = 0.4, blimit = 40000, tlimit = 1000):
        # hold data in matrix format
        # we'll decide how large this matrix ought to be based on distribution samples
        # on average the growth should be something like child^(t/waiting)
        # number of time points ought to be (t/waiting)
        # we'll overestimate this using the percentile_threshold buffer
        child_over = self.child_sample[int(self.test_samples * (1.0 - prct_thresh))]
        waiting_under = self.waiting_sample[int(self.test_samples * prct_thresh)]
        max_branches = int(child_over ** (t / waiting_under))
        max_timepoints = int(t / waiting_under)

        # print max_branches,max_timepoints

        # strict upper limits are set by blimit and tlimit
        max_branches = min(max_branches,blimit)
        max_timepoints = min(max_timepoints,tlimit)
        
        # print max_branches,max_timepoints

        # construct matrix to hold data
        # rows are data for each branch, columns are times of division/birth
        # initialize with infs
        data = np.empty((max_branches,max_timepoints))
        data.fill(np.inf)
        
        # we'll also keep a vector of indices that indicate how many divisions/births
        # have occurred for a branch
        next_index = np.array([0 for i in xrange(max_branches)])

        # initial conditions for single cell
        data[0,0] = 0.0
        next_index[0] = 1

        # we're just gonna go branch by branch and divide until we're past t
        current_branch = 0
        last_branch = 0
        while current_branch <= last_branch:
            # divide current branch until theres a division time past t
            # time of last division on branch
            last_div = data[current_branch,next_index[current_branch]-1]
            while last_div < t:
                # time of next division
                div_time = last_div + self.waiting_gen()
                # record time of next division
                data[current_branch,next_index[current_branch]] = div_time
                next_index[current_branch] += 1

                if next_index[current_branch] == max_timepoints:
                    print 'timepoint max exceeded'
                    # double max timepoints
                    bonus_timepoints = np.empty((max_branches,max_timepoints))
                    bonus_timepoints.fill(np.inf)
                    data = np.append(data,bonus_timepoints,axis=1)
                    max_timepoints *= 2

                # if our division occurs before t, we need to add a new branch
                if div_time < t:
                    # increase number of active branches by 1
                    last_branch += 1
                    if last_branch == max_branches:
                        print 'branch max exceeded'
                        # double max branches
                        bonus_branches = np.empty((max_branches,max_timepoints))
                        bonus_branches.fill(np.inf)
                        data = np.append(data,bonus_branches,axis=0)
                        next_index = np.append(next_index,np.zeros(max_branches))
                        max_branches *= 2
                    # start new branch at division time
                    data[last_branch,next_index[last_branch]] = div_time
                    next_index[last_branch] += 1
                
                # update last division time
                last_div = div_time
            
            # update current branch
            current_branch += 1
        print last_branch
        
        
        
GW = GWEngine(np.random.exponential,lambda: 2)


