#!/usr/bin/env python -u

import numpy as np
import matplotlib.pyplot as plt
import pw2py as pw
import glob
import scipy.constants
import warnings


warnings.warn("This script is not tested and may be wrong!!!")


# constants for conversions
hbar = scipy.constants.physical_constants['reduced Planck constant in eV s'][0]
eV2J = scipy.constants.e
atomic_mass = scipy.constants.atomic_mass


def read_full_dat() -> dict:
    ''' read linear extrapolation data '''
    full_dat = {}
    for lin in sorted(glob.glob('lin*/')):
        full_dat[lin] = []
        for rat in sorted(glob.glob(f'{lin}ratio-*/')):
            ratio = float(rat.split('-')[-1][:-1])
            energy = pw.qeout.final_energy(f'{rat}scf.out')[0]
            full_dat[lin].append([ratio, energy])
        full_dat[lin] = np.array(full_dat[lin])

    # shift minimum to zero
    minE = min([dat.min() for _, dat in full_dat.items()])
    for lin, dat in full_dat.items():
        dat[:, 1] -= minE
        np.savetxt(f'{lin}_barrier.dat', dat)
    full_dat

    return full_dat


def calc_zpl(full_dat: dict, report=False) -> dict:
    ''' calculate position of minimum for each curve '''
    zpl = {}
    for lin, dat in full_dat.items():
        zpl[lin] = (dat[np.argmin(dat[:, 1]), 0], dat[:, 1].min())
        if report:
            print(f'{lin} -> min = {zpl[lin][1]:.4f} at ratio = {zpl[lin][0]:.4f}')
    return zpl


def calc_dQ(full_dat: dict, report=False) -> float:
    ''' calculate dQ from atomic geometry of minimum (units of amu^1/2 Ang) '''
    lins = list(full_dat.keys())
    geo0 = pw.atomgeo.from_file(
        f'{lins[0]}/ratio-{zpl[lins[0]][0]:.4f}/scf.out')
    geo1 = pw.atomgeo.from_file(
        f'{lins[1]}/ratio-{zpl[lins[1]][0]:.4f}/scf.out')
    dQ = geo0.calc_dQ(geo1)
    if report:
        print(f'dQ = {dQ:.4f}')
    return dQ


def calc_polyfit(full_dat: dict) -> dict:
    ''' calculate degree 2 polynomial fit around the minimum '''
    polyfit = {}
    for lin, dat in full_dat.items():
        fit_dat = dat[:6] if zpl[lin][0] < 0.5 else dat[-6:]
        polyfit[lin] = np.polyfit(*fit_dat.T, 2)
    return polyfit


def calc_hw_and_S(polyfit: dict, dQ: float, report=False) -> tuple:
    ''' calculate hbar*omega (hw) in eV and the Huang-Rhys factor (S) '''
    dQ_SI = dQ * (atomic_mass)**(1/2) * 1e-10
    hw, S = {}, {}
    for lin, fit in polyfit.items():
        # hbar omega in eV
        hw[lin] = hbar * (fit[0] * eV2J / (dQ_SI**2))**(1/2)
        # HR factor (unitless)
        S[lin] = hw[lin] * (dQ_SI**2) / hbar**2 / 2 / eV2J
        if report:
            print(f'{lin} -> hw = {hw[lin]:.4f}, S = {S[lin]:.4f}')
    return hw, S


def make_plot(full_dat: dict, polyfit: dict, zpl: dict, dQ: float, filename: str = 'config.png'):
    ''' make configuration plot '''
    plt.figure(dpi=300, figsize=(4, 4))
    xlim = np.array([-0.5, 1.5])
    colors = {k: c for k, c in zip(full_dat, ['red', 'blue'])}

    # plot scatter data
    for lin, dat in full_dat.items():
        qvals = dat[:, 0]*dQ
        yvals = dat[:, 1]
        plt.scatter(qvals, yvals, color=colors[lin])
    # plot fit lines
    for lin, fit in polyfit.items():
        qvals = np.linspace(-0.5, 1.5)
        yvals = fit[0] * qvals**2 + fit[1] * qvals + fit[2]
        qvals *= dQ
        plt.plot(qvals, yvals, color=colors[lin])

    # labels and such
    plt.ylabel('Energy [eV]')
    plt.xlabel(r'$\rm \Delta Q$ [$\sqrt{\rm amu}$ $\rm \AA$]')
    plt.xlim(xlim*dQ)
    # keep default ylim
    ax = plt.gca()
    ylim = ax.get_ylim()
    plt.ylim(ylim)

    # extract zpl for annotating
    z0, z1 = list(zpl.values())
    # lines
    kwargs = {'color': 'grey', 'linestyle': 'dotted', 'zorder': 0}
    plt.hlines((z0[1], z1[1]), *(xlim*dQ), **kwargs)
    plt.vlines((z0[0]*dQ, z1[0]*dQ), *(ylim), **kwargs)

    # save fig
    plt.tight_layout()
    plt.savefig('config.png')
    plt.close()


if __name__ == '__main__':
    full_dat = read_full_dat()
    zpl = calc_zpl(full_dat, report=True)
    dQ = calc_dQ(full_dat, report=True)
    polyfit = calc_polyfit(full_dat)
    hw, S = calc_hw_and_S(polyfit, dQ, report=True)
    make_plot(full_dat, polyfit, zpl, dQ)
