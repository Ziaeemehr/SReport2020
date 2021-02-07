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


def runCommand(arg):
    command = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}".format(*arg)
    system("./prog " + command)
#---------------------------------------------------------#

def batchRun():
    arg = []
    for k in couplings:
        for gm_ in gm:
            for gahp_ in gahp:
                for i in I:
                    for nt in networks:
                        arg.append([N,
                                    tSimulation,
                                    tTransition,
                                    k,
                                    gm_,
                                    gahp_,
                                    gSynaptic,
                                    i,
                                    nt,
                                    numEnsemble,
                                    printSpikes,
                                    printVoltages,
                                    printSynapse])
    Parallel(n_jobs=n_jobs)(
        map(delayed(runCommand), arg))


#---------------------------------------------------------#
N = 20
couplings = [1.0] # np.linspace(0,0.5, 21)
I = np.arange(1, 3, 0.1)
gm = [1.0]
gahp = [0.0]
gSynaptic = 0.01
tTransition = 500.0
tSimulation = tTransition + 1000.0
numEnsemble = 1
# networks = ['SF-BDG', 'SF-BDG-T', "SF-DAG", "SF-DAG-T"]
networks = ["complete"]
printSpikes = 1
printVoltages = 1
printSynapse = 1
n_jobs = 4
#---------------------------------------------------------#


if __name__ == "__main__":

    seed = 125
    start = time()
    adj = nx.to_numpy_array(nx.complete_graph(N))
    np.savetxt("dat/complete.dat", adj, fmt="%d")
    batchRun()
    lib.display_time(time()-start)