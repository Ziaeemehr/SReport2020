import numpy as np
import pylab as pl
import networkx as nx
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sys import exit
import os


def generate_random_graph(n, p, seed=125, fname="C"):
    G = nx.erdos_renyi_graph(n, p, seed=seed, directed=False)
    adj = nx.to_numpy_array(G)
    np.savetxt("dat/"+fname+".dat", adj, fmt='%d')


def display_time(time):
    ''' 
    show real time elapsed
    '''

    hour = int(time/3600)
    minute = (int(time % 3600))/60
    second = time-(3600.*hour+60.*minute)
    print "Done in %d hours %d minutes %09.6f seconds" \
        % (hour, minute, second)

#---------------------------------------------------------#


def read_from_file(filename, n):
    spikes = []
    with open(filename, "r") as f:
        for i in range(n):
            data = f.readline().split()
            data = [float(j) for j in data]
            spikes.append(data)
    return spikes
# ------------------------------------------------------------------#


"""

def pairSpikeSynchrony(spikes1, spikes2, nbins=50):
    # REF: Gamma Oscillation by Synaptic Inhibition in a Hippocampal
    # Interneuronal Network Model

    def f(spikes):
        x = np.histogram(spikes, bins=nbins)
        x = np.where(x > 0.5, 1, 0)
        return x

    x = f(spikes1)
    y = f(spikes2)

    numerator = np.sum(x * y)
    denuminator = np.sqrt(np.sum(x) * np.sum(y))

    return numerator / denuminator
# ------------------------------------------------------------------#


def averageSynchrony(spiketrains, nbins):
    
    #spiketrains as list of list 
    #N = len(spiketrains) shows the number of neurons


    N = len(spiketrains)
    synchMat = np.identity(N)

    for i in range(N):
        for j in range(i + 1, N):
            synchMat[i][j] = synchMat[j][i] = pairSpikeSynchrony(
                spiketrains[i], spiketrains[j])

    return np.mean(synchMat)
# ------------------------------------------------------------------#


def calculatePeriod(spiketrains):

    N = len(spiketrains)
    T = 0.0
    counter = 0
    for i in range(N):
        if len(spiketrains[i]) > 3:
            T += spiketrains[i][-1] - spiketrains[i][-2]
            counter += 1
    T /= float(counter)
    tau = np.linspace(1, T, 10, endpoint=True)
    print "Period : ", T

    return T, tau

"""
