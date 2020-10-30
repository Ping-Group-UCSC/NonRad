#!/usr/bin/env python
import os
import re
import numpy as np
from math import sqrt, exp, pi
from numpy.linalg import norm
from constant import Ry2eV, Electron2Coulomb, Ang2m, AMU2kg, hbar_eVs, hbar_Js, AMU2me, Ang2Bohr, Ha2eV, Kelvin2au, \
    Bohr2m, Bohr2Ang, second2au, GHz2Ha
import glob

from chem_utils import f_Element_Symbol_to_Mass
from libreadqe import read_pos_and_etot_ratio, read_wave, get_ratio_folder, read_eig, get_save_folder
from libmath import overlap_x_quantum_harmonic_num, overlap_x_quantum_harmonic_ladder_hr
from constant import indent


def calc_dQ(vecR, list_pos1, list_pos2):
    '''
    Calculate dQ

    input: vecR in angstrom (lattice vector in columns), list_pos1 in crystal coordinate
    output unit: Amu1/2 * Angstrom

    :return: dQ, Delta Q between two structure
            M,  effective mass
            ar_delta : delta of all atoms (effective phonon mode) , 3*N array (xyz change first, atom index second)
            ar_Qi : normalized phono mode, (1/sqrt(M) * M_i * mu ) for all i-th directions,
                mu is the phonon mode vector, in 3*N array (xyz change first, atom index second)
    '''
    list_delta_pos = [atom1["pos"] - atom2["pos"] for atom1, atom2 in zip(list_pos1, list_pos2)]

    for ix in range(len(list_delta_pos)):
        v = list_delta_pos[ix]
        for i in range(len(v)):
            x = v[i]
            if (abs(x) > 0.5):  # Shift out of boundary ; shift back
                v[i] = x - round(x, 0)
        list_delta_pos[ix] = np.dot(vecR, v)

    list_delta = [{"species": atom["species"], "delta":pos} for atom, pos in zip(list_pos1, list_delta_pos)]

    list_R = [norm(atom["delta"])**2 for atom in list_delta]
    list_mR = [f_Element_Symbol_to_Mass(atom["species"]) * norm(atom["delta"])**2 for atom in list_delta]
    dQ2 = sum(list_mR)
    dR2 = sum(list_R)
    dQ = sqrt(dQ2)
    M = float(dQ2 / dR2)

    ar_delta = np.asarray([atom["delta"] for atom in list_delta]).flatten()
# Normalized
    ar_Qi = np.asarray([atom["delta"] * f_Element_Symbol_to_Mass(atom["species"]) / sqrt(M) / sqrt(dR2)
                        for atom in list_delta]).flatten()

    return dQ, M, ar_delta, ar_Qi


def calc_dE(dir1):
    '''
    Calculate the energy difference between ratio=0 and ratio=1 for given electron configuration
    Always report the absolute value
    Note, the calculation may not contain that information, return 0 fi that is the case
    '''
    list_data1, vecR = read_pos_and_etot_ratio(dir1)
    etot1 = None
    etot2 = None
    for x in list_data1:
        if (x["ratio"] == 0):
            etot1 = x["etot"]
        elif (x["ratio"] == 1):
            etot2 = x["etot"]
    if (etot1 is None or etot2 is None):
        return 0
    else:
        return abs(etot1 - etot2)


def calc_freq(folder, ratio_min, ratio_max, dQ):
    '''
    Fit the effective 1D frequency from scf total energy and Q
    Data are selected from ratio_min to ratio_max

    raises AssertionError if minimum is not at 0 or 1
    '''
    import warnings
    # ignore warning by polyfit
    warnings.filterwarnings("ignore", message="Polyfit may be poorly conditioned")

    tol = 1e-6
    list_data = []
    for filename in glob.glob("%s/ratio-*/scf.out" % folder):
        ratio = float(re.match(".*ratio-(.+)/scf.out", filename).group(1))
        if (ratio > ratio_min - tol and ratio <= ratio_max + tol):
            # Read energy
            # Always read the first converged energy from output (so this script works with relaxation calculation)
            try:
                with open(filename.replace(".in", ".out"), 'r') as f:
                    lines = f.readlines()
                    for line in lines[::-1]:
                        if ("total energy" in line and "is the sum" not in line):
                            tag = line[:2].strip()
                            break

                    for line in lines:
                        if (line.startswith(tag)):
                            etot = float(line.split()[-2]) * Ry2eV
                            break
            except:  # noqa: E722
                raise ValueError("Error reading %s" % (filename.replace(".in", ".out")))
            list_data.append({"ratio": ratio, "etot": etot})

    if (len(list_data) == 0):
        for filename in glob.glob("%s/ratio-*/POSCAR" % folder):
            ratio = float(re.match(".*ratio-(.+)/POSCAR", filename).group(1))
# Read energy
            if (ratio > ratio_min - tol and ratio <= ratio_max + tol):
                with open(filename.replace("POSCAR", "OUTCAR"), 'r') as f:
                    lines = f.readlines()
                for i in range(len(lines)-1, 0, -1):
                    if ("energy  without entropy" in lines[i]):
                        etot = float(lines[i][65:82])
                        break

                list_data.append({"ratio": ratio, "etot": etot})

    print("%sFound: %i points in folder %s in range %.2f ~ %.2f" %
          (indent*2, len(list_data), folder, ratio_min, ratio_max))

    if (len(list_data) < 3):
        print("Not enough points for range %.4f ~ %.4f" % (ratio_min, ratio_max))
        return None, None, None

    ar_ratio = np.asarray([x["ratio"] for x in list_data])
    ar_etot = np.asarray([x["etot"] for x in list_data])

    # ar_ratio_min = ar_ratio[np.argmin(ar_etot)]
    # min_at_0_or_1 = np.any(np.isclose(ar_ratio_min, [0, 1], atol=tol))
    # assert min_at_0_or_1, "Minimum must be at ratio = 0 or 1, true min at: %f" % ar_ratio_min

# Fit with different orders
# Maxmumly one order less than number of points
# Note it is angular frequency
    list_result = []
    for order in range(2, min(5, len(list_data)+1)):
        p = np.polyfit(ar_ratio, ar_etot, deg=order)
        try:
            hfreq = 1/dQ * sqrt(p[order-2] * 2 * Electron2Coulomb / (Ang2m**2 * AMU2kg)) * hbar_eVs
        except ValueError:
            hfreq = 0

        list_result.append({"order": order, "hbarfreq": hfreq})
# Return order 2 as final result and
# estimate error

    freq = list_result[0]["hbarfreq"]
    ar = [x["hbarfreq"] for x in list_result]
    freq_error = max(ar) - min(ar)

    S = dQ**2 * (Ang2m**2 * AMU2kg) * (freq / hbar_eVs) / 2 / hbar_Js

    return freq, freq_error, list_result, S


def calc_wif_standalone(dir1, dir2, ix_defect, ix_bandmin, ix_bandmax, de=None):
    '''
    Compute Wif from a series of SCF calculations; Check Phys. Rev. B 90, 075202(2014) for definition.
    Each SCF calculation must be in ratio-* folder and contains a *.save folder as result

    This function prints everything to stdout

    :param dir1: Folder contains ratio-* for initial state (hole in VB), only the ratio=0 will be used
    :param dir2: Folder contains ratio-* for final state (hole in defect, or the lower energy of two)
    :param defect: Defect band index (1-base)
    :param bandmin: Valence hole band index minimum(1-based)
    :param bandmax: Valence hole band index maximum(1-based)
    :param de: energy difference for two levels; if not specified,
        use energy difference at Q=0 (ratio=0) of the final stat
    '''

    list_data1, vecR = read_pos_and_etot_ratio(dir1)
    list_data2, vecR = read_pos_and_etot_ratio(dir2)

    dE = de
    dic_band_overlap = dict((ix, []) for ix in range(ix_bandmin, ix_bandmax+1))

# Compute dQ
    dQ, M, ar_delta, ar_Qi = calc_dQ(vecR, list_data2[0]["pos"], list_data2[-1]["pos"])
    dQ = dQ / (list_data2[-1]["ratio"] - list_data2[0]["ratio"])
    print("DeltaQ = %.2f" % dQ)

    folder1_q0 = get_save_folder(os.path.join(dir1, get_ratio_folder(0)))
    # ar_eig1_q0 = read_eig(folder1_q0)
    folder2_q0 = get_save_folder(os.path.join(dir2, get_ratio_folder(0)))
    ar_eig2_q0 = read_eig(folder1_q0)
    spin1 = 2
    spin2 = 2

    for iband in range(ix_bandmin, ix_bandmax+1):
        evc1 = read_wave(folder2_q0, spin2, 1, iband)
        if (de is None):
            dE = -ar_eig2_q0[0, spin1-1, iband-1, 0] + ar_eig2_q0[0, spin2-1, ix_defect-1, 0]
        else:
            dE = de
        print("ei-ef : %.2f" % dE)

        # list_ratio = []
        for data in list_data1:
            ratio = data["ratio"]
            print("Ratio %.4f" % ratio)
            folder1 = get_save_folder(os.path.join(dir1, get_ratio_folder(ratio)))
            folder2 = get_save_folder(os.path.join(dir2, get_ratio_folder(ratio)))

            if (folder2 is None):
                print("Skip ratio %.4f as missing .save folder %s" % (ratio, dir1 if folder1 is None else dir2))
                continue

            ar_eig2 = read_eig(folder2)
            if (ar_eig2 is None):
                print("Cannot read %s" % folder1)
                continue
            if (ar_eig2 is None):
                print("Cannot read %s" % folder2)
                continue

            evc2 = read_wave(folder2, spin2, 1, ix_defect)
            s = np.dot(np.conj(evc1), evc2)
# Gamme (real wavefunciton) only contains positive G; double it
            if evc1.dtype == np.float64:
                s *= 2
            print("%.6f %24.16g %24.16g %24.16g" % (
                dE,
                abs(s),
                s.real, s.imag))
            dic_band_overlap[iband].append([ratio, abs(s)])

        dic_band_overlap[iband].sort(key=lambda x: x[0])

    for iband, list_overlap in dic_band_overlap.items():
        # Only fit first several
        ar0 = np.asarray(list_overlap)
        ar_ratio = ar0[:, 0]
        ar_overlap = ar0[:, 1]
        ar_Q = dQ * ar_ratio
        nq = 3
        print("Fitting %i points: Q=[%.3f, %.3f]" % (nq, ar_Q[0], ar_Q[nq-1]))
        p = np.polyfit(ar_Q[:nq], ar_overlap[:nq],  deg=1)
        print("Band %i dS/dQ %.2e Wif %.2e" % (iband, p[0], p[0] * dE))

    with open("overlap.txt", 'w') as f:
        list_band_overlap = list(dic_band_overlap.items())
        list_band_overlap.sort(key=lambda x: x[0])
        for iband, list_overlap in list_band_overlap:
            f.write("# %i\n" % iband)
            np.savetxt(f, np.asarray(list_overlap))
            f.write("\n")


def calc_wif(dir_i, dir_f, ix_defect, ix_bandmin, ix_bandmax, dQ, de=None, spinname="down"):
    '''
    Compute Wif from a series of SCF calculations; Check Phys. Rev. B 90, 075202(2014) for definition.
    Each SCF calculation must be in ratio-* folder and contains a *.save folder as result

    :param dir_i: Folder contains ratio-* for initial state (hole in VB), only the ratio=0 will be used
    :param dir_f: Folder contains ratio-* for final state (hole in defect, or the lower energy of two)
    :param defect: Defect band index (1-base)
    :param bandmin: Valence hole band index minimum(1-based)
    :param bandmax: Valence hole band index maximum(1-based)
    :param de: energy difference for two levels; if not specified,
        use energy difference at Q=0 (ratio=0) of the final stat
    :param spin: compute the overlap of wavefunction of given spin channel, can be either "up" or "down".
        For a defect  all filled initial state and a hole in spin down channel final state (PRB 90, 075202),
        the "down" spin should be investigate. For a defect-defect transition the selected defect
        level spin should be used.

    :return: dic_eig : dic(intial/final) = array of [ratio, eigenvalues for different bands],
        array of [ratio, occupations numbers for different bands]
    '''
    if (ix_bandmin > ix_bandmax):
        raise ValueError("Band max must be equal or larger than band min")

    list_data1, vecR = read_pos_and_etot_ratio(dir_i)
    list_data2, vecR = read_pos_and_etot_ratio(dir_f)

    dE = de
    dic_band_overlap = dict((ix, []) for ix in range(ix_bandmin, ix_bandmax+1))

# Collect all eigenvalues
    ix_plotmin = min(ix_bandmin, ix_defect) - 1
    ix_plotmax = max(ix_bandmax, ix_defect) + 1

# Two state and two spin seperately
    dic_eig_all = {("i", 1): [],
                   ("i", 2): [],
                   ("f", 1): [],
                   ("f", 2): []
                   }

# Read all eigenvalues
    for statename, list_data, dir0 in [
            ("i", list_data1, dir_i),
            ("f", list_data2, dir_f),
    ]:
        for data in list_data:
            ratio = data["ratio"]
            folder = get_save_folder(os.path.join(dir0, get_ratio_folder(ratio)))
            ar_eig1 = read_eig(folder)

            for spin in (1, 2):
                dic_eig_all[(statename, spin)].append(
                    [ratio, ar_eig1[0, spin-1, ix_plotmin - 1:ix_plotmax, 0],
                        ar_eig1[0, spin-1, ix_plotmin-1:ix_plotmax, 1]]
                )

# Organize data
    dic_eig_occ = {}
    for key, val in dic_eig_all.items():
        dic_eig_occ[key] = {"eig": np.asarray([[x[0]] + x[1].tolist() for x in val]),
                            "occ": np.asarray([[x[0]] + x[2].tolist() for x in val])}

# Compute dE
    # The hole is always spin-down as in QE number of electrons is always more in spin up
    spin = {"up": 1, "down": 2}[spinname]

    folder_q0 = get_save_folder(os.path.join(dir_f, get_ratio_folder(0)))
    ar_eig1_q0 = read_eig(folder_q0)

    '''
    loop through bands and compute overlap
    '''
    dic_band_overlap = {}
    print("%sGathering :  band ratio" % (indent*2))
    for iband in range(ix_bandmin, ix_bandmax+1):
        list_overlap = []
#       print(spin, iband, ix_defect, ar_eig1_q0[0,spin-1,iband-1,0], + ar_eig1_q0[0, spin-1, ix_defect-1,0])
        dE = float(-ar_eig1_q0[0, spin-1, iband-1, 0] + ar_eig1_q0[0, spin-1, ix_defect-1, 0])
# Read wavefunction of a perturbed bulk state in final state
# (This and defect wavefunction should be same Hamiltonian to be meaningful)
        evc1 = read_wave(folder_q0, ispin=spin, ik=1, ib=iband)

        # list_ratio = []
        for data in list_data2:
            ratio = data["ratio"]
            folder2 = get_save_folder(os.path.join(dir_f, get_ratio_folder(ratio)))

            if (folder2 is None):
                print("Skip ratio %.4f as missing .save folder %s" % (ratio, os.path.join(dir_f, "ratio-%.4f" % ratio)))
                continue
            print("%sOverlap: %i  %.4f" % (indent*3, iband, ratio))
            # print("%sCompute overlap for band %i ratio %.4f" % (indent*2, iband, ratio))

# Read wavefunction of a defect state of final state
            evc2 = read_wave(folder2, spin, 1, ix_defect)
# Note : if QE is gamma only (wavefunction is real), the evc contains only positive G
# And evc^2 = 0.5
# If QE is not gamma only then all G are included, evc^2 = 1
# Note G=0 is doubled ; but error very small in general
            s = np.dot(np.conj(evc1), evc2)
            if evc1.dtype == np.float64:
                s *= 2
            list_overlap.append([ratio, abs(s)])

        list_overlap.sort(key=lambda x: x[0])

        dic_band_overlap[iband] = np.asarray(list_overlap)

    list_wif = []
    for iband, ar0 in dic_band_overlap.items():
        # Only fit first several
        # print(ar0)
        ar_Q = ar0[:, 0] * dQ
        ar_overlap = ar0[:, 1]
        nq = 3
#       print("Fitting %i points: Q=[%.3f, %.3f]" % (nq, ar_Q[0], ar_Q[nq-1]))
        p = np.polyfit(ar_Q[:nq], ar_overlap[:nq],  deg=1)
#       print("Band %i dS/dQ %.2e Wif %.2e" % (iband, p[0], p[0] * dE))
        list_wif.append((iband, float(p[0] * dE)))

# Find the maxium
    list_wif.sort(key=lambda x: x[1])
    ixband_wifmax, wif = list_wif[-1]
    dic_wif = dict(list_wif)

    return dic_eig_occ, dE, dic_band_overlap, dic_wif, ixband_wifmax, wif


def calc_phonon_part_T0_HR(dE, dQ, freq, order_x):
    '''
    Calculate the phonon part in Eq 22 at T=0 (only ni=0 included)
    And assume freqi=freqf
    '''
    freqi = freq
    freqf = freq

    dQ = dQ * AMU2me ** 0.5 * Ang2Bohr
    dE = dE / Ha2eV
    freqi = freqi / 1000 / Ha2eV
    freqf = freqf / 1000 / Ha2eV
    print("%satomic unit dQ : %24.14g" % (indent*2, dQ))

# Set parameters
    sigma = 0.8 * freqf
    sigmamax = 4 * freqf
    # sf = dE / freqf
    # si = dE / freqi

# Iterate
    s = 0
    n = 0
    m_center = (dE + n * freqi) / freqf
    mmin = max(0, int(m_center - sigmamax / freqf - 1))
    mmax = int(m_center + sigmamax / freqf + 1)
#           print("Maximum final state oscillator quantum number %4i %4i" % (mmin, mmax))
    s1 = 0
    print("%s  n   m   deltacoef      overlap" % (indent*(3+3)))
    for m in range(mmax, mmin, -1):
        dEph = dE + n * freqi - m * freqf
        coef = exp(- dEph ** 2 / (2 * sigma ** 2)) / (sigma * sqrt(2 * pi))
        overlap = overlap_x_quantum_harmonic_ladder_hr(m, dQ, 1, freqi, order_x)

        s1 += overlap**2 * coef
        print("%sContribution: %3i %3i %14.7g %14.7g" % (indent*3, n, m, coef, overlap**2))
    s += s1

    if order_x == 1:
        print(r"%s<\chi_n|Q-Q0|\chi_m> @ T = %d ~=  %24.14g" % (indent*2, 0, s))
    elif order_x == 0:
        print(r"%s<\chi_n|\chi_m> @ T = %d ~=  %24.14g" % (indent*2, 0, s))
    return s


def calc_phonon_part(dE, dQ, freqi, freqf, list_T, order_x):
    '''
    Calculate the phonon part in Eq 22

    :return: list of temperature-phonon part
    '''
    dQ = dQ * AMU2me ** 0.5 * Ang2Bohr
    dE = dE / Ha2eV
    freqi = freqi / 1000 / Ha2eV
    freqf = freqf / 1000 / Ha2eV
#   M = M * AMU2me
#   ar_delta /= Bohr2Ang
#   ar_Qi = ar_Qi / Bohr2Ang * (AMU2me)**0.5
    print("%satomic unit dQ : %24.14g" % (indent*2, dQ))
    print("%satomic unit freqi/f: %24.14g %24.14g" % (indent*2, freqi, freqf))
#   print("M freqi frqef : %24.14g %24.14g %24.14g" % (M, freqi, freqf))
#   print("M*freq i/f : %24.14g %24.14g" % (M*freqi, M*freqf))

    list_T = [x for x in list_T]
# Sort descending to reuse some intergrals (high temperature contains more states and integrals)
    list_T.sort(reverse=True)

# Cached function for overlap
    cache_overlap_x = {}

    def get_overlap_x(ni, mf):
        '''
        Cached overlap
        Note mass is already included in dQ so no mass term here
        '''
        i = (ni, mf)
        if (i not in cache_overlap_x):
            #           print(ni, mf, dQ, freqi, freqf)
            cache_overlap_x[i] = overlap_x_quantum_harmonic_num(
                ni, mf, dQ, 1, freqi, 1, freqf, order_x)

        return cache_overlap_x[i]

# Set parameters
    sigma = 0.8 * freqf
    sigmamax = 4 * freqf
    # sf = dE / freqf
    si = dE / freqi


# Enumerate all temperatures
    list_all = []
    for T in list_T:
        kBT = Kelvin2au * T

# Esimate maximum quantum numbers we must included
        tr = 1e-6

# Bose-Einstein
#       n_be0 = 1/(exp(0.5 * freqi / kBT) - 1)
#       nmax = min(int(si), 5)
#       for i in range(1, nmax+1):
#           n_be1 = 1/(exp((i+0.5) * freqi / kBT) - 1)
#           if (n_be1 / n_be0 < tr):
#               nmax = i - 1
#               break

# Boltzman
        nmax = max(int(si), 5)
        for i in range(1, nmax+1):
            n_b1 = exp(-i * freqi / kBT)
            if (n_b1 < tr):
                nmax = i - 1
                break
        print("\n%sTemperature : %4i K ; Maximum initial state oscillator quantum number %4i" % (indent*2, T, nmax))

# Iterate
        s = 0
        print("%s  n   m   occupation     deltacoef      overlap" % (indent*(3+3)))
        for n in range(nmax, -1, -1):
            n_b = exp(-n * freqi/kBT) * (1 - exp(-freqi / kBT))
            m_center = (dE + n * freqi) / freqf
            mmin = max(0, int(m_center - sigmamax / freqf - 1))
            mmax = int(m_center + sigmamax / freqf + 1)
#           print("Maximum final state oscillator quantum number %4i %4i" % (mmin, mmax))
            s1 = 0
            for m in range(mmax, mmin, -1):
                dEph = dE + n * freqi - m * freqf
                coef = exp(- dEph ** 2 / (2 * sigma ** 2)) / (sigma * sqrt(2 * pi))
                overlap = get_overlap_x(n, m)

                s1 += overlap**2 * coef
                print("%sContribution: %3i %3i %14.7g %14.7g %14.7g" %
                      (indent*3, n, m, n_b, coef, overlap**2))
            s += s1 * n_b

        if order_x == 1:
            print(r"%s<\chi_n|Q-Q0|\chi_m> @ T = %24.14g  %24.14g" % (indent*2, T, s))
        elif order_x == 0:
            print(r"%s<\chi_n|\chi_m> @ T = %24.14g  %24.14g" % (indent*2, T, s))

        list_all.append((T, s))

    return list_all


def calc_cp_T(dim, vol, g, wif, list_phonon_part):
    '''
    Calculate cp from components
    Unit in meV, Ang and AMU

    :return: unit in cm^* s
    '''
    vol = vol / Bohr2Ang**dim
    wif = wif / Ha2eV / (AMU2me**0.5 / Bohr2Ang)

    list_cp = []
    for temperature, osc in list_phonon_part:
        cp = vol * 2 * pi * g * wif ** 2 * osc * (Bohr2m**dim) * 1E2**dim * second2au
        list_cp.append([temperature, cp])

    return list_cp


def calc_lifetime_T(g, wif, list_phonon_part, job):
    '''
    Calculate lifetime (as concentration of calculation)
    '''

    '''
    Some notes on units here!
    if job == 'nonrad':
        cp = 2pi/hbar * wif^2 * xif
        wif [energy/dQ] = [eV/(sqrt(AMU)*Ang)]
        xif [deltaQ^2/energy] = [(me*bohr^2)/Ha]
    elif job == 'isc':
        cp = 4pi*hbar * wif^2 * xif
        NOTE wif is lamba_t in this case
        wif [frequency] = [GHz]
        xif [1/energy] = [1/Ha]

    '''

    # determine the prefactor
    if job == 'nonrad':
        # wif [eV/(sqrt(AMU)*Ang)] -> [Ha/(sqrt(me)*bohr)]
        wif = wif / Ha2eV / (AMU2me**0.5 / Bohr2Ang)
        prefactor = 2 * pi * g * wif ** 2

    elif job == 'isc':
        # convert GHz to Ha (note this multiplying hbar^2 in given units)
        wif *= GHz2Ha
        prefactor = 4 * pi * g * wif ** 2

    list_lifetime = []
    for temperature, osc in list_phonon_part:
        # note second2au is Ha -> Hz (1/hbar in given units)
        cp = prefactor * osc * second2au

        list_lifetime.append([temperature, 1/cp])

    return list_lifetime
