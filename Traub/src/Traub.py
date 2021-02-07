from scipy.integrate import odeint
import networkx as nx
from time import time
from numpy import exp
import numpy as np
import pylab as pl
import sys


class Traub:
    V_REST = -95.0
    G_SYN = 0.01

    def __init__(self, seed=None):
        if seed:
            np.random.seed(seed)

    def set_params(self, **par):
        self.N = par["N"]
        self.Imin, self.Imax = par["I"]
        self.gm = par['gm']
        self.gahp = par['gahp']
        self.adj = par['adj']
        self.tInjection = par["tInjection"]

        dimension = 7
        v0 = np.random.uniform(-95.0, -55.0, size=self.N).tolist()
        # v0 = np.random.uniform(-95.0, -94.0, size=self.N).tolist()
        self.initial_condition = np.asarray(v0 +
                                            np.random.rand(self.N * (dimension - 1)).tolist())

    def ode_sys(self, x0, t):
        '''
        define Traub Model
        '''
        N = self.N
        v = x0[:N]
        m = x0[N: 2 * N]
        n = x0[2 * N: 3 * N]
        h = x0[3 * N: 4 * N]
        w = x0[4 * N: 5 * N]
        ca = x0[5 * N: 6 * N]
        s = x0[6 * N:]

        dv = np.zeros(N)
        dm = np.zeros_like(dv)
        dn = np.zeros_like(dv)
        dh = np.zeros_like(dv)
        dw = np.zeros_like(dv)
        dca = np.zeros_like(dv)
        ds = np.zeros_like(dv)

        for i in range(N):

            mlInf = 1.0 / (1.0 + exp(-(v[i] + 25.0) / 2.5))
            am = 0.32 * (v[i] + 54.0) / (1.0 - exp(-(v[i] + 54.0) / 4.0))
            bm = 0.28 * (v[i] + 27.0) / (exp((v[i] + 27.0) / 5.0) - 1.0)
            ah = 0.128 * exp(-(v[i] + 50.0) / 18.0)
            bh = 4.0 / (1.0 + exp(-(v[i] + 27.0) / 5.0))
            an = 0.032 * (v[i] + 52.0) / (1.0 - exp(-(v[i] + 52.0) / 5.0))
            bn = 0.5 * exp(-(v[i] + 57.0) / 40.0)
            tw = 100.0 / (3.3 * exp((v[i] + 35.0) /
                                    20.0) + exp(-(v[i] + 35.0) / 20.0))
            wInf = 1.0 / (1.0 + exp(-(v[i] + 35.0) / 10.0))
            ica = mlInf * (v[i] - 120.0)

            I_syn = 0.0
            for j in np.nonzero(self.adj[i, :])[0]:
                I_syn += self.G_SYN * s[j] * (self.V_REST)
                # I_syn += self.G_SYN * s[j] * (v[j])
            I_syn *= - 1.0

            if t > self.tInjection:
                I = self.Imax
                gm = self.gm
                gahp = self.gahp

            else:
                I = self.Imin
                gm = 0.0
                gahp = 0.0

            dv[i] = I - 100.0 * h[i] * m[i] ** 3 * (v[i] - 50.0) - (80.0 * n[i] ** 4 + gm * w[i] +
                                                                    gahp * (ca[i] / (ca[i] + 1.0))) * (v[i] + 100) - 0.2 * (v[i] + 67.0) - ica + I_syn
            dm[i] = am * (1.0 - m[i]) - bm * m[i]
            dn[i] = an * (1.0 - n[i]) - bn * n[i]
            dh[i] = ah * (1.0 - h[i]) - bh * h[i]
            dw[i] = (wInf - w[i]) / tw
            dca[i] = -0.002 * ica - ca[i] / 80.0
            ds[i] = 2.0 * (1.0 - s[i]) / \
                (1.0 + exp(-(v[i] + 10.0) / 10.0)) - 0.1 * s[i]

        return np.concatenate((dv, dm, dn, dh, dw, dca, ds), axis=0)

    def simulate(self, tSimulation, dt):

        t = np.arange(0, tSimulation, dt)
        x0 = self.initial_condition
        sol = odeint(self.ode_sys, x0, t)

        return {"t": t, "x": sol}


if __name__ == "__main__":

    start = time()

    N = 5
    adj = nx.to_numpy_array(nx.complete_graph(N))
    # adj = np.loadtxt("dat/FF.dat")

    params = {
        "N": N,
        "I": [6.0, 6.0],
        "tInjection": 200.0,
        "gm": 1.0,
        "gahp": 0.0,
        "adj": adj,
    }

    traub = Traub(seed=1248)
    traub.set_params(**params)
    sol = traub.simulate(600, 0.02)

    fig, ax = pl.subplots(2, figsize=(10, 4), sharex=True)
    ax[0].plot(sol["t"], sol["x"][:, :N])
    ax[1].plot(sol["t"], np.sum(sol["x"][:, (6 * N):], axis=1))
    ax[1].set_xlabel("Time (ms)", fontsize=20)
    ax[0].set_ylabel(r"$V_i$", fontsize=20)
    ax[1].set_ylabel(r"$\sum s_i$", fontsize=20)
    pl.tight_layout()
    pl.savefig(str("../data/fig/%s.png" % params["I"][1]))
    print "Done in %g seconds. " % (time()-start)
    pl.show()
