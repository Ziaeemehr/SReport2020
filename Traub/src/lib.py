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
    minute = (int(time % 3600))//60
    second = time-(3600.*hour+60.*minute)
    print ("Done in %d hours %d minutes %09.6f seconds" \
        % (hour, minute, second))

#---------------------------------------------------------#


def read_from_file(filename, n):
    spikes = []
    with open(filename, "r") as f:
        for i in range(n):
            data = f.readline().split()
            data = [float(j) for j in data]
            spikes.append(data)
    return spikes
