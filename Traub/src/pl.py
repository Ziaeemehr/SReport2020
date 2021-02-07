from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import numpy as np
import pylab as pl
import os
import lib
import sys
pl.switch_backend('agg')
#---------------------------------------------------------#
os.chdir('../data/')
directories = ['fig', 'npz']
for i in directories:
    if not os.path.exists(i):
        os.makedirs(i)

#---------------------------------------------------------#


def plot_rasters(K, I_app, Phi):

    T = np.arange(t_trans, tfinal, 0.05)
    nk = len(K)
    for cn in networks:
        for p in Phi:
            for ip in I_app:
                for k in range(nk):
                    fig, ax = pl.subplots(1, figsize=(10, 5))
                    ifname = str('%s-%.6f-%.6f-%.6f' % (cn, K[k], p, ip))
                    spikes = lib.read_from_file('text/spk-'+ifname+'.txt', N)

                    for ii in range(N):
                        ax.plot(spikes[ii], [ii] *
                                len(spikes[ii]), '.', c='royalblue',
                                markersize=2)

                    ax.set_ylabel("Node index")
                    ax.set_xlabel('Time(ms)')
                    ax.set_xlim(t_trans, tfinal)
                    # ax.set_xlim(2500, 2900)
                    ax.set_title(
                        cn + str(', K = %g, phi = %g, I = %g' % (K[k], p, ip)))

                    fig.savefig('fig/r-'+ifname+'.png')
                    pl.close()

# ------------------------------------------------------- #


def plot_phase_space(R, X, Y, name="R", xtickstep=1, ytickstep=1,
                     xlabel=None, ylabel=None, title=None,
                     vmax=None, vmin=None):
    '''
    plot R in 2D plane of X and Y axises
    '''
    print(len(X), len(Y), R.shape)

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    r, c = R.shape
    assert((r > 1) & (c > 1))

    x_step = X[1] - X[0]
    y_step = Y[1] - Y[0]

    f, ax = pl.subplots(1, figsize=(6, 6))
    im = ax.imshow(R, interpolation='nearest', aspect='auto',
                   cmap='afmhot', vmax=vmax, vmin=vmin)
    ax.invert_yaxis()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pl.colorbar(im, cax=cax)

    ax.set_xticks(np.arange(0, len(X), xtickstep))
    ax.set_xticklabels(str("%.1f" % i)for i in X[::xtickstep])
    ax.set_yticks(np.arange(0, len(Y), ytickstep))
    ax.set_yticklabels(str("%.1f" % i)for i in Y[::ytickstep])

    if xlabel:
        ax.set_xlabel(xlabel, fontsize=16)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=16)
    if title:
        ax.set_title(title, fontsize=16)

    pl.savefig("fig/"+name+".png")
    pl.close()
# ------------------------------------------------------- #


def plot_I_couplings(I, couplings):

    nk = len(couplings)
    nI = len(I)

    for cn in networks:
        for i in range(nI):
            v = np.zeros(nk)
            s = np.zeros(nk)
            p = np.zeros(nk)
            fig, ax = pl.subplots(3, figsize=(8, 8), sharex=True)
            for kk in range(len(couplings)):
                ifname = str('%s-%.6f-%.6f' %
                             (cn, couplings[kk], I[i]))
                C = np.genfromtxt('text/par-'+ifname+'.txt')
                if (numEnsemble > 1):
                    v[kk] = np.mean(C[:, 2])
                    s[kk] = np.mean(C[:, 3])
                    p[kk] = np.mean(C[:, 4])
                    stdv[kk] = np.std(C[:, 2])
                    stds[kk] = np.std(C[:, 3])
                    stdp[kk] = np.std(C[:, 4])
                else:
                    v[kk] = C[2]
                    s[kk] = C[3]
                    p[kk] = C[4]

            np.savez('npz/' + cn + str('-%.6f-%.6f' % (couplings[kk], I[i])),
                     couplings=couplings[kk],
                     I=I,
                     voltage=v,
                     spike=s,
                     phase=p)

            if numEnsemble > 1:
                ax[0].errorbar(couplings, v, yerr=stdv)
                ax[1].errorbar(couplings, s, yerr=stds)
                ax[2].errorbar(couplings, p, yerr=stdp)
            else:
                ax[0].plot(couplings, v)
                ax[1].plot(couplings, s)
                ax[2].plot(couplings, p)

            ax[0].set_ylabel("voltage synchrony")
            ax[1].set_ylabel("spike synchrony")
            ax[2].set_ylabel("phase difference")
            pl.tight_layout()
            pl.savefig(str("fig/f-%.6f.png" % I[i]))
            pl.close()


def plotVoltageSpike(couplings,
                     I,
                     gm,
                     nCurves=3,
                     xlim=None):

    nk = len(couplings)
    nI = len(I)
    ngm = len(gm)

    t = np.arange(tTransition + 0.02, tSimulation, 0.02)
    colors = pl.cm.brg(np.linspace(0, 1, nCurves + 1))

    for cn in networks:
        for kk in range(len(couplings)):
            for j in range(ngm):
                for i in range(nI):

                    fig = pl.figure(figsize=(15, 5))
                    gs1 = gridspec.GridSpec(6, 1, hspace=0.)
                    ax = []
                    ax.append(fig.add_subplot(gs1[0:2]))
                    ax.append(fig.add_subplot(gs1[2:4]))
                    ax.append(fig.add_subplot(gs1[4:]))

                    subname = str('%s-%.6f-%.6f' %
                                  (cn, couplings[kk], I[i]))

                    v = np.loadtxt('text/v-' + subname + '.txt')
                    syn = np.loadtxt("text/s-" + subname + ".txt")
                    spks = []
                    with open('text/spk-'+subname+'.txt', 'r') as fin:
                        for line in fin:
                            s = line.split()
                            s = map(float, s)
                            spks.append(s)

                    for ii in range(nCurves):
                        ax[0].plot(t, v[ii, :], lw=1,
                                   label=ii + 1, c=colors[ii])
                    for ii in range(N):
                        ax[2].plot(spks[ii], [ii+1]*len(spks[ii]), 'k.')

                    ax[1].plot(t, np.sum(syn, axis=0), lw=1,
                               label=ii + 1, c=colors[ii])

                    for i1 in range(3):
                        # ax[i1].set_xlim(tfinal-500, tfinal)
                        # ax[i1].set_xlim(39500, 40000)
                        ax[i1].tick_params(labelsize=15)
                        # ax[i1].set_xticks([tfinal-500, tfinal])
                    ax[0].set_ylabel('V (mV)', fontsize=16)
                    ax[1].set_ylabel(r'$\sum s$', fontsize=16)
                    ax[1].set_xlabel('Time (ms)', fontsize=16)
                    # ax[1].set_ylabel('Index', fontsize=16)
                    # ax[0].set_title(cn+'L', fontsize=16)
                    ax[0].set_xticks([])
                    ax[1].set_xticks([])
                    ax[2].set_yticks([])

                    if xlim:
                        for jj in range(3):
                            ax[jj].set_xlim(xlim)
                    else:
                        for jj in range(3):
                            ax[jj].set_xlim(np.min(t), np.max(t))
                    # ax[1].axis('off')
                    # ax[1].set_ylim(0.5,3.5)
                    ax[0].legend()

                    pl.savefig("fig/f-"+subname+".png")
                    pl.close()


if __name__ == "__main__":

    from run import *

    #plot_rasters(K, I_app, Phi)
    # plot_I_couplings(I, couplings)
    plotVoltageSpike(couplings, I, gm, nCurves=int(sys.argv[1]))
