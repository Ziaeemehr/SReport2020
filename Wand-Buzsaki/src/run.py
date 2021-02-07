import os
import numpy as np
import networkx as nx
from os import system
from time import time
import lib
try:
    from threading import Thread
    from joblib import Parallel, delayed
except:
    pass
#---------------------------------------------------------#


def run_command(arg):
    command = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}".format(
        *arg)
    system("./prog " + command)
#---------------------------------------------------------#


def Linux():
    arg = []
    for nt in networks:
        for k in K:
            for i in I_app:
                for pi in Phi:
                    for alpha_ in alpha:
                        for beta_ in beta:
                            for E_syn_ in E_syn:
                                arg.append([N,
                                            k * nt[1],
                                            tFinal,
                                            tTransition,
                                            i,
                                            pi,
                                            alpha_,
                                            beta_,
                                            E_syn_,
                                            nt[0],
                                            numEnsembles,
                                            PRINE_SPIKES,
                                            PRINT_VOLTAGES,
                                            PRINT_SYNAPSES,
                                            MEASURE_FREQUENCY])
    Parallel(n_jobs=n_jobs)(
        map(delayed(run_command), arg))


#---------------------------------------------------------#
N = 200
K = [1.25]
Phi = [6.0, 7.0, 8.0, 9.0]
I_app = [0.35]
# Phi = np.arange(4.6, 9.1, 0.1)
# I_app = np.arange(0.1, 3.0, 0.1)
alpha = [12.0]
# tau = np.arange(0.5, 3.5, 0.5)
# beta = 1.0 / tau
beta = [0.1]
E_syn = [-75.0]
tTransition = 1000.0
tFinal = tTransition + 2000.0
numEnsembles = 1
networks = [["SF-DAG", 2]]
PRINE_SPIKES = 1
PRINT_VOLTAGES = 1
PRINT_SYNAPSES = 1
MEASURE_FREQUENCY = 1
n_jobs = 4
#---------------------------------------------------------#


if __name__ == "__main__":

    seed = 125
    start = time()
    # lib.generate_random_graph(N, 1.0, seed=seed, fname="C")
    Linux()
    lib.display_time(time()-start)


# I_mean = [42.16, 50.41]
# I_sigma = [0.0,0.0]
# networks = ['SF-BDG', 'SF-DAG',
# 			'SF-U', 'SW-BDG',
# 			'SW-DAG', 'SW-U']
