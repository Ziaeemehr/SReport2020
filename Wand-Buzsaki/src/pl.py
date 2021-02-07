from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import numpy as np
import pylab as pl
import os
import lib
from run import *
pl.switch_backend('agg')
#---------------------------------------------------------#
os.chdir('../data/')
directories = ['fig', 'npz']
for i in directories:
    if not os.path.exists(i):
        os.makedirs(i)

#---------------------------------------------------------#


def plot_rasters(K, I_app, Phi):

    T = np.arange(tTransition, tFinal, 0.05)
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
                    ax.set_xlim(tTransition, tFinal)
                    ax.set_xlim(2500, 2900)
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
    print len(X), len(Y), R.shape

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


def plot_synchrony_I_Phi(K, I_app, Phi, beta, xtickstep=1, ytickstep=1):

    nk = len(K)
    nphi = len(Phi)
    nI = len(I_app)

    for cn in networks:
        for kk in range(len(K)):
            for be in beta:
                v, s, p = np.zeros((nphi, nI)), np.zeros(
                    (nphi, nI)), np.zeros((nphi, nI))
                for j in range(nphi):
                    for i in range(nI):
                        ifname = str('%s-%.6f-%.6f-%.6f-%.6f' %
                                     (cn[0], K[kk]*cn[1], Phi[j], I_app[i], be))
                        C = np.genfromtxt('text/par-'+ifname+'.txt')
                        if (numEnsembles > 1):
                            v[j, i] = np.mean(C[:, 2])
                            s[j, i] = np.mean(C[:, 3])
                            p[j, i] = np.mean(C[:, 4])
                        else:
                            try:
                                v[j, i] = C[2]
                                s[j, i] = C[3]
                                p[j, i] = C[4]
                            except:
                                print "error in reading phi=%g, I=%g" % (
                                    Phi[j], I[i])

                np.savez('npz/'+cn[0]+str('-%.3f' % (K[kk]*cn[1])), K=K[kk]*cn[1], phi=Phi, I=I_app,
                         voltage=v, spike=s, phase=p)

                plot_phase_space(v, I_app, Phi, cn[0]+'-VS-' +
                                 str('%.3f-%.3f' % (K[kk]*cn[1], be)), vmax=1, vmin=0,
                                 xtickstep=xtickstep, ytickstep=ytickstep)

                # plot_phase_space(s, I_app, Phi, cn+'-Spike_synchrony-' +
                #                  str('%.3f' % K[kk]), vmax=1, vmin=0,
                #                  xtickstep=2, ytickstep=2)
                # plot_phase_space(p, Phi, I_app, cn+'-Phase_synchrony')


def plot_voltage_spike(K, I_app, Phi, beta, xlim=None):

    nk = len(K)
    nI = len(I_app)
    nphi = len(Phi)
    dt = 0.05
    STEP_PRINT = 5
    t = np.arange(tTransition + STEP_PRINT * dt, tFinal, STEP_PRINT * dt)
    colors = pl.cm.Blues(np.linspace(0, 1, N + 1))

    for cn in networks:
        for kk in range(len(K)):
            for j in range(nphi):
                for i in range(nI):
                    for b in beta:

                        fig = pl.figure(figsize=(10, 5))
                        gs1 = gridspec.GridSpec(6, 1, hspace=0.)
                        ax = []
                        ax.append(fig.add_subplot(gs1[0:3]))
                        ax.append(fig.add_subplot(gs1[3:5]))
                        ax.append(fig.add_subplot(gs1[5:]))

                        subname = str('%s-%.6f-%.6f-%.6f-%.6f' %
                                      (cn[0], K[kk]*cn[1], Phi[j], I_app[i], b))

                        v = np.loadtxt('text/v-'+subname+'.txt')
                        spks = []
                        with open('text/spk-'+subname+'.txt', 'r') as fin:
                            for line in fin:
                                s = line.split()
                                s = map(float, s)
                                spks.append(s)

                        for ii in range(0, N, 15):
                            # print len(t), len(v[ii, :])
                            ax[0].plot(t, v[ii, :], lw=2,
                                       label=ii+1, c=colors[ii])
                        syn = np.loadtxt("text/syn-" + subname + ".txt")
                        syn = np.sum(syn, axis=1)
                        np.savez("npz/"+subname, t=t, syn=syn)
                        ax[2].plot(t, syn, lw=2, c="k")
                        ax[2].set_ylabel(r"$\sum s_i$")

                        for ii in range(N):
                            ax[1].plot(spks[ii], [ii+1]*len(spks[ii]), 'k.')

                        for i1 in range(3):
                            ax[i1].tick_params(labelsize=15)
                        if xlim:
                            for i1 in range(3):
                                ax[i1].set_xlim(xlim)
                            # ax[i1].set_xticks([tfinal-500, tfinal])
                        ax[2].set_xlabel('Time (ms)', fontsize=16)
                        ax[0].set_ylabel('V (mV)', fontsize=16)
                        ax[1].set_ylabel("Node index", fontsize=16)
                        ax[0].set_title(r"$\phi$=%.1f" %
                                        Phi[j], fontsize=16)
                        ax[0].set_xticks([])
                        ax[1].set_xticks([])
                        ax[1].set_yticks([1, 200])
                        ax[1].set_ylim([0, 201])
                        ax[0].set_ylim([-70,15])
                        ax[0].set_yticks([-55, 0])
                        ax[0].axhline(y=-55, lw=2, ls="--", c="gray")
                        ax[2].margins(y=0.1)
                        ax[2].set_xticks([2700,2800,2900, 3000])
                        # ax[1].axis('off')
                        # ax[1].set_ylim(0.5,3.5)
                        # ax[0].legend()
                        pl.savefig("fig/f-"+subname+".eps")
                        pl.close()


def plot_synchrony(K, I_app, Phi,
                   xlabel=None, ylabel=None):
    nk = len(K)
    nx = len(x)
    nI = len(I_app)

    for cn in networks:
        for kk in range(len(K)):
            fig, ax = pl.subplots(2, figsize=(6, 4))
            for i in range(nI):
                Rv = np.zeros(nx)
                Rs = np.zeros(nx)
                Rp = np.zeros(nx)
                Fr = np.zeros(nx)
                for j in range(nx):
                    ifname = str('%s-%.6f-%.6f-%.6f-%.6f' %
                                 (cn, K[kk], x[j], I_app[i], beta[0]))
                    C = np.genfromtxt('text/par-'+ifname+'.txt')
                    if (numEnsembles > 1):
                        Rv[j] = np.mean(C[:, 2])
                        Rs[j] = np.mean(C[:, 3])
                        Rp[j] = np.mean(C[:, 4])
                        Fr[j] = np.mean(C[:, 5])

                    else:
                        try:
                            Rv[j] = C[2]
                            Rs[j] = C[3]
                            Rp[j] = C[4]
                            Fr[j] = C[5]
                        except:
                            print "error in reading x=%g, I=%g" % (
                                x[j], I_app[0])

                subname = str('%s-%.6f-%.6f-%.6f' %
                              (cn, K[kk], I_app[i], beta[0]))
                # np.savez('npz/' + subname, K=K[kk], phi=Phi, I=I_app[i])

                ax[0].plot(x, Rv, lw=2, marker="o", label=I_app[i])
                ax[1].plot(x, Fr, lw=2, marker="s", label=I_app[i])

            ax[1].set_xlabel(xlabel)
            ax[0].set_ylabel("R")
            ax[1].set_ylabel("Frequency [Hz]")
            ax[0].legend()
            ax[1].legend()

            pl.tight_layout()
            pl.savefig(str("fig/R-%.6f.png" % I_app[i]))
            pl.close()


def plot_synchrony_tau(K, I_app, Phi, x,
                       xlabel=None, ylabel=None):
    nk = len(K)
    nI = len(I_app)
    nx = len(x)

    for cn in networks:
        for kk in range(len(K)):
            fig, ax = pl.subplots(2, figsize=(6, 6))
            for i in range(nI):
                Rv = np.zeros(nx)
                Rs = np.zeros(nx)
                Rp = np.zeros(nx)
                Fr = np.zeros(nx)
                for j in range(nx):
                    ifname = str('%s-%.6f-%.6f-%.6f-%.6f' %
                                 (cn, K[kk], Phi[0], I_app[i], x[j]))
                    C = np.genfromtxt('text/par-'+ifname+'.txt')
                    if (numEnsembles > 1):
                        Rv[j] = np.mean(C[:, 2])
                        Rs[j] = np.mean(C[:, 3])
                        Rp[j] = np.mean(C[:, 4])
                        Fr[j] = np.mean(C[:, 5])

                    else:
                        try:
                            Rv[j] = C[2]
                            Rs[j] = C[3]
                            Rp[j] = C[4]
                            Fr[j] = C[5]
                        except:
                            print "error in reading phi=%g, I=%g" % (
                                x[j], I_app[0])

                # subname = str('%s-%.6f-%.6f-%.6f' %
                #               (cn, K[kk], I_app[i], Phi[0]))
                # np.savez('npz/' + subname, K=K[kk], Rs=Rs, Rv=Rv, Rp=Rp, Fr=Fr,
                #          phi=Phi, I=I_app[i], beta=beta)

                ax[0].plot(1.0 / x, Rv, lw=2, marker="o", label=I_app[i])
                ax[1].plot(1.0 / x, Fr, lw=2, marker="s", label=I_app[i])

            ax[1].set_xlabel(r"$\tau_{syn}$")
            ax[0].set_ylabel("R")
            ax[1].set_ylabel("Frequency [Hz]")
            ax[0].legend()
            ax[1].legend()

            pl.tight_layout()
            pl.savefig(str("fig/R-%.6f.png" % I_app[i]))
            pl.close()


# def plotPairSynchrony(net, K, Phi, I_app):

#     subname = str('%s-%.6f-%.6f-%.6f' %
#                   (net, K, Phi, I_app))

#     spiketrains = []
#     with open('text/spk-'+subname+'.txt', 'r') as fin:
#         for line in fin:
#             s = line.split()
#             s = map(float, s)
#             spiketrains.append(s)

#     T, tau = lib.calculatePeriod(spiketrains)

#     synchrony = []
#     for i in tau:
#         nbins = int(round(T / i))
#         synchrony.append(lib.averageSynchrony(spiketrains, nbins))


#     np.savez("../data/npz/synchrony-%s"%subname, zip(tau, synchrony))
#     pl.plot(tau, synchrony)
#     pl.savefig("../data/fig/synchrony-%s.png" % subname)
#     pl.close()
if __name__ == "__main__":

    # plot_rasters(K, I_app, Phi)
    # plot_synchrony_I_Phi(K, I_app, Phi, beta)
    # plot_synchrony_tau(K,
    #                    I_app,
    #                    Phi,
    #                    beta,
    #                    xlabel=r"$\phi$",
    #                    ylabel="synchrony")
    plot_voltage_spike(K, I_app, Phi, beta, [tFinal-300, tFinal])

    # for cn in networks:
    #     for phi in Phi:
    #         for k in K:
    #             for i_app in I_app:
    #                 plotPairSynchrony(cn, k, phi, i_app)
