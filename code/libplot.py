#!/usr/bin/env python
# Contains functions to plot
import numpy as np
import matplotlib.pyplot as plt

try:
    # python2
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle  # noqa: F401

list_color = ["black", "red", "blue", "darkgreen", "orange", "magenta", "darkgoldenrod", "coral",
              "lightskyblue", "darkolivegreen", "wheat", "lightcoral", "blueviolet", "lightseagreen", "gold"]
list_linestyle = ["-", "--", "--", "--", "--", "--", "--"]
list_dash = [(), (1, 1), (2, 2), (3, 3), (3.5, 3.5), (4, 4), (4.5, 4.5)]
list_marker = ['o', 's', 'd', 'v', '^', '<', '>', '*', 'x']

# plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
plt.rcParams["font.size"] = 18
plt.rcParams["legend.fontsize"] = "small"
# rc('text.latex', preamble=r'\usepackage{cmbright}')


def plot_eig_Q(filename, ar_eig, ar_occ, dQ):
    '''
    Plot eigenvalues against Q
    '''
    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
#   np.savetxt(os.path.join(dir_image, "%s.dat" % casename), data)

    ar_dQ = ar_eig[:, 0] * dQ

    for i in range(ar_eig.shape[1]-1):
        all_occ = (ar_occ[:, i+1] >= 0.5).all()
        all_unocc = (ar_occ[:, i+1] < 0.5).all()
        if (all_unocc):
            markerfacecolor = "white"
            fillstyle = "full"
        elif (all_occ):
            markerfacecolor = "black"
            fillstyle = "full"
        else:  # Half filled, should not happen in a well defined defect state
            print("Half filled band: band #%i" % i)
            markerfacecolor = "black"
            fillstyle = "bottom"

        ax.plot(ar_dQ, ar_eig[:, i+1], marker="o", linestyle="-", color="black",
                fillstyle=fillstyle, markerfacecolor=markerfacecolor)

    ax.set_xlabel(r"Q amu^{1/2} Ang")
    ax.set_ylabel(r"eigenvalue (eV)")
#   ax.set_xlim(xmin, xmax)
#   ax.set_xlim(xmin, xmax)
    plt.tight_layout()

    fig.savefig(filename, dpi=150)


def plot_overlap_Q(filename, dic_band_overlap, dic_wif, dQ, dE):
    '''
    Plot <psi|psi> with Q and fitting curve
    '''
    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

# Limit to first 4 points
    nmax = 4

    ymin = 0
    ymax = 0
    for ix, (iband, ar0) in enumerate(dic_band_overlap.items()):
        # Plot points
        ar_Q = ar0[:, 0] * dQ
        ax.plot(ar_Q, ar0[:, 1], "bo", marker="o", color=list_color[ix],
                fillstyle="full", markerfacecolor="white", label="Band %i" % iband)
        ax.plot(ar_Q, dic_wif[iband] / dE * ar_Q, linestyle='-', color=list_color[ix])
        ymin = min(ar0[:nmax, 1].min(), ymin)
        ymax = max(ar0[:nmax, 1].max(), ymax)

    ax.set_xlabel(r"Q amu^{1/2} Ang")
    ax.set_ylabel(r"Wavefunction overlap")
# Limit to first 3 point
    dx = ar_Q[nmax-1] - ar_Q[0]
    ax.set_xlim(ar_Q[0] - dx * 0.1, ar_Q[nmax-1] + dx * 0.1)
    dy = ymax - ymin
    ax.set_ylim(ymin - dy*0.1, ymax + dy * 0.1)
    plt.tight_layout()

    fig.savefig(filename, dpi=150)


def plot_cp_T(filename, list_T_cp):
    '''
    Plot cp-T curve
    '''
    plt.clf()
    list_T = [float(x[0]) for x in list_T_cp]
    list_cp = [x[1] for x in list_T_cp]
    plt.semilogy(list_T, list_cp, color="black")
    plt.ylabel("\\tilde{C}_p (cm^3/s)")
    plt.xlabel("T (K)")
    plt.tight_layout()
    plt.savefig(filename, dpi=150)


def plot_tot_Q(filename, list_data_i, list_data_f, dE_corrected, dQ, enable_plot):
    '''
    Plot total energy v.s. Q

    :return: the data plotted (energy shifted to correct position, minimum set to 0)
    '''
    e0 = min([x["etot"] for x in list_data_f])
    dE = min([x["etot"] for x in list_data_i]) - e0
# Force dE as 1eV to align charge
#   dE_corrected = 1

    ar_ratio_f = np.array([x["ratio"] for x in list_data_f])
    ar_q_f = dQ * ar_ratio_f
    ar_ratio_i = np.array([x["ratio"] for x in list_data_i])
    ar_q_i = dQ * ar_ratio_i

    ar_i = np.array([ar_q_i,
                     [x["etot"]-e0-(dE - dE_corrected) for x in list_data_i]
                     ]).T
    ar_f = np.array([ar_q_f,
                     [x["etot"]-e0 for x in list_data_f]
                     ]).T

    if enable_plot:
        qmax = max(max(ar_q_f), max(ar_q_i))

        xmin = -qmax
        xmax = 2*qmax
        ymin = -dE_corrected/2
        ymax = 2*dE_corrected

        fig, ax = plt.subplots(1, 1, figsize=(6, 4))

        ax.plot(ar_i[:, 0], ar_i[:, 1], marker="o", linestyle="-", color="red", fillstyle="none")
        ax.plot(ar_f[:, 0], ar_f[:, 1], marker="o", linestyle="-", color="blue", fillstyle="none")

    #   ax.plot(ar_q_i, [x["etot"]-e0-(dE - dE_corrected) for x in list_data_i],
    #               marker="o", linestyle="-", color="red", fillstyle = "none")
    #   ax.plot(ar_q_f, [x["etot"]-e0 for x in list_data_f], marker="o", linestyle="-",
    #               color="blue", fillstyle = "none")
        ax.set_xlabel(r"Q amu$^{1/2}$\AA")
        ax.set_ylabel("energy (eV)")
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.plot([0, 0], [ymin, ymax], linestyle="--", color="black")
        ax.plot([qmax, qmax], [ymin, ymax], linestyle="--", color="black")
        ax.plot([xmin, xmax], [0, 0], linestyle="--", color="black")
        ax.plot([xmin, xmax], [dE_corrected, dE_corrected], linestyle="--", color="black")
        fig.savefig(filename, dpi=150)

    return ar_i, ar_f
