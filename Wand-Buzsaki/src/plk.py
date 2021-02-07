import numpy as np
import pylab as pl
from run import *
pl.switch_backend('agg')
from sys import exit
import os 

#---------------------------------------------------------#
os.chdir('../data/')
directories = ['fig', 'npz']
for i in directories:
    if not os.path.exists(i):
        os.makedirs(i)
#---------------------------------------------------------#
nk = len(K)
for cn in networks:
    for p in Phi:
        for ip in I_app:

            fig, ax = pl.subplots(2, figsize=(7,7))
            
            vs, br = np.zeros(nk), np.zeros(nk)
            std_vs, std_br = np.zeros(nk), np.zeros(nk)

            for k in range(nk):
                ifname = str('%s-%.6f-%.6f-%.6f' % (cn, K[k], p, ip))
                C = np.genfromtxt('par-'+ifname+'.txt')
                if num_sim ==1:
                    vs[k] = C[2]
                    br[k] = C[3]
                else:
                    vs[k] = np.mean(C[:,2])
                    std_vs[k] = np.std(C[:,2])

                    br[k] = np.mean(C[:,3])
                    std_br[k] = np.std(C[:,3])
            if num_sim ==1:
                ax[0].plot(K, vs, lw=3)
                ax[1].plot(K, br, lw=3)
            else:
                ax[0].errorbar(K, vs, yerr=std_vs, fmt="--o",
                                ecolor='r', elinewidth=1)
                ax[1].errorbar(K, br, yerr=std_br, fmt="--o",
                                ecolor='r', elinewidth=1)
            fname = str('%s-%.6f-%.6f' % (cn, p, ip))

            ax[0].set_ylim(0, 1)
            ax[1].set_ylim(0, 1)
            
            ax[0].set_ylabel(r"Vol synch")
            ax[1].set_ylabel(r"Burst synch")

            fig.savefig('fig/f-'+fname+'.png')
            np.savez('npz/'+fname, K=K, phi=p, I=ip, 
            vs=vs, br=br, std_vs=std_vs, std_br=std_br)
#---------------------------------------------------------#


                
