#!/usr/bin/env python
# Connect all script to calculate
import sys
import os
import yaml
import numpy as np
from numpy.linalg import norm

np.seterr(all="log")

from libreadqe import read_pos_and_etot_ratio
from libcalc import (
    calc_dQ,
    calc_freq,
    calc_wif,
    calc_phonon_part,
    calc_cp_T,
    calc_phonon_part_T0_HR,
    calc_lifetime_T,
    calc_dE,
)

# from libplot import plot_tot_Q, plot_eig_Q, plot_overlap_Q, plot_cp_T

from io_package import read_cell_and_pos_auto


def saveinput():
    with open(file_input, "w") as f:
        yaml.dump(dinput, f, default_flow_style=False, width=200)


def save():
    with open(file_store, "w") as f:
        yaml.dump(data, f, default_flow_style=False, width=200)


def save_dic_ndarray(filename, dic):
    """
    Same a dictionary with named ndarray
    """
    with open(filename, "w") as f:
        for name, ar in dic.items():
            f.write("# %s\n" % name)
            np.savetxt(f, ar)
            f.write("\n")


def savedata_single_array(filename, database, term, ar):
    """
    Save an array to a file and record the filename in the database
    """
    np.savetxt(filename, ar)
    database[term] = filename
    # Call save all after
    save()


def main():
    file_input = "input.yaml"
    file_store = "nonrad.yaml"

    if not os.path.exists(file_input):
        dinput = yaml.load(
            """
    title: Nonradiative
    dimension : The periodical dimension of the system, can be "xyz" (3D) or "xy" (2D)
    g : The degeneracy factor of sites, mostly 1
    dE : The energy difference between two states (charge transition level) as initial - final
    folder_init_state : The folder name of scf calculation (ratio*) with initial state charge
    folder_final_state : The folder name of scf calculation (ratio*) with final state charge
    ratio_init : The number to represent the initial state equalibrium geometry in ratio*, should be 0 or 1
    ratio_final : The number to represent the final state equalibrium geometry in ratio*, should be 0 or 1
    bulkband_index : The bulk band indicies in format [min,max], start from 1
    defectband_index : The defect band index
    """
        )
    else:
        with open(file_input, "r") as f:
            dinput = yaml.load(f.read())

    if "The " in dinput["folder_init_state"]:
        print("Please modify %s according to instruction inside" % file_input)
        saveinput()
        sys.exit(1)

    if not os.path.exists(file_store):
        data = {}
    else:
        with open(file_store, "r") as f:
            data = yaml.load(f.read())

    # Set default input values
    if "defectband_spin" not in dinput:
        dinput["defectband_spin"] = "down"  # This is for calculation as PRB 90,075202

    if not "unit" in data:
        data["unit"] = {
            "dQ": "aMU^1/2 Ang",
            "M": "aMU",
            "hbarFreqi": "meV",
            "hbarFreqi_error": "meV",
            "hbarFreqf": "meV",
            "hbarFreqf_error": "meV",
        }

    # Compute dQ
    if not "dQ" in data and "folder_init_state" in dinput:
        print("Compute dQ ...")
        (vecR, list_pos_f), _ = read_cell_and_pos_auto(
            dinput["folder_final_state"] + "/ratio-0.0000/scf"
        )
        (vecR, list_pos_i), _ = read_cell_and_pos_auto(
            dinput["folder_final_state"] + "/ratio-1.0000/scf"
        )
        dQ, M, _, _ = calc_dQ(vecR, list_pos_f, list_pos_i)
        data["dQ"] = dQ
        data["M"] = M

    # Plot Q-etot curve
    if not "filename_Q_etot" in data:
        filename = "tot_Q.png"
        data["filename_Q_etot"] = filename
        print("Plot Q-etot ...")
        list_data_i, vecR = read_pos_and_etot_ratio(dinput["folder_init_state"])
        list_data_f, vecR = read_pos_and_etot_ratio(dinput["folder_final_state"])
        # plot_tot_Q(filename, list_data_i, list_data_f, dinput["dE"], data["dQ"])

    if not "dErelf" in data:
        dEf = calc_dE(dinput["folder_final_state"])
        dEi = calc_dE(dinput["folder_init_state"])
        data["dErelf"] = dEf
        data["dEreli"] = dEi

    # Compute frequency
    if not "hbarFreqi" in data:
        # From 0-0.15 and 0-0.20
        r1 = dinput["ratio_init"]
        r2 = abs(r1 - 0.15)
        if r1 > r2:
            r1, r2 = (r2, r1)
        freq, freq_error, list_result, S = calc_freq(
            dinput["folder_init_state"], ratio_min=r1, ratio_max=r2, dQ=data["dQ"]
        )
        data["hbarFreqi"] = freq * 1000
        data["hbarFreqi_error"] = freq_error * 1000
        data["hbarFreqi_order"] = list_result
        data["S_i"] = S

        r1 = dinput["ratio_final"]
        r2 = abs(r1 - 0.15)
        if r1 > r2:
            r1, r2 = (r2, r1)

        freq, freq_error, list_result, S = calc_freq(
            dinput["folder_final_state"], ratio_min=r1, ratio_max=r2, dQ=data["dQ"]
        )
        data["hbarFreqf"] = freq * 1000
        data["hbarFreqf_error"] = freq_error * 1000
        data["hbarFreqf_order"] = list_result
        data["S_f"] = S

    if not "hbarFreqi_to02_order" in data:
        r1 = dinput["ratio_init"]
        r2 = abs(r1 - 0.20)
        if r1 > r2:
            r1, r2 = (r2, r1)
        freq, freq_error, list_result, S = calc_freq(
            dinput["folder_init_state"], ratio_min=r1, ratio_max=r2, dQ=data["dQ"]
        )
        data["hbarFreqi_to02_order"] = list_result

        r1 = dinput["ratio_final"]
        r2 = abs(r1 - 0.20)
        if r1 > r2:
            r1, r2 = (r2, r1)
        freq, freq_error, list_result, S = calc_freq(
            dinput["folder_final_state"], ratio_min=r1, ratio_max=r2, dQ=data["dQ"]
        )
        data["hbarFreqf_to02_order"] = list_result

    if not isinstance(dinput["bulkband_index"], list) or not isinstance(
        dinput["defectband_index"], int
    ):
        print("Please modify %s according to instruction inside" % file_input)
        saveinput()
        sys.exit(1)

    if "Wif" not in data or "WifBand" not in data:
        # Calculate Wif
        # Plot and store data
        dic_eig_occ, diffEigQ0, dic_band_overlap, dic_wif, ixband_wifmax, wif = calc_wif(
            dinput["folder_init_state"],
            dinput["folder_final_state"],
            dinput["defectband_index"],
            dinput["bulkband_index"][0],
            dinput["bulkband_index"][1],
            data["dQ"],
            de=None,
            spinname=dinput["defectband_spin"],
        )

        for (statename, spin), val in dic_eig_occ.items():
            casename = "%s-%s" % (statename, "up" if spin == 1 else "down")
            filename = "eig-%s.txt" % casename
            np.savetxt(filename, val["eig"])
            filename = "occ-%s.txt" % casename
            np.savetxt(filename, val["occ"])
            data["filename_eigvals_%s" % casename] = filename
            # plot_eig_Q("eig-%s.png" % casename, val["eig"], val["occ"], data["dQ"])

        filename = "overlap.txt"
        data["filename_overlap"] = filename
        save_dic_ndarray(filename, dic_band_overlap)

        # plot_overlap_Q("overlap.png", dic_band_overlap, dic_wif, data["dQ"], diffEigQ0)

        data["diffEigQ0"] = diffEigQ0
        data["Wif"] = wif
        data["WifBand"] = ixband_wifmax
        data["Wif_Allband"] = dic_wif

    if not isinstance(dinput["dE"], float) or "The " in dinput["dimension"]:
        print("Please modify %s according to instruction inside" % file_input)
        saveinput()
        sys.exit(1)

    if "filename_phonon_part" not in data:
        list_phonon_part = calc_phonon_part(
            dinput["dE"],
            data["dQ"],
            data["hbarFreqi"],
            data["hbarFreqf"],
            dinput["temperature"],
        )

        print("phonon part done ... saving to phonon_part.txt")
        savedata_single_array(
            "phonon_part.txt", data, "filename_phonon_part", list_phonon_part
        )
    else:
        print("phonon part read from phonon_part.txt")
        list_phonon_part = np.loadtxt("phonon_part.txt")

    # Calculate the apporximate phonon part at T=0
    if "phonon_part_T0_HR" not in data:
        print("Calculating the approximate phonon part at T=0")
        data["phonon_part_T0_HR"] = calc_phonon_part_T0_HR(
            dinput["dE"], data["dQ"], data["hbarFreqf"]
        )

    if "filename_Cp_tilde" in data:
        list_cp = np.loadtxt("cp_T.txt")
    else:
        # Read and compute volume
        (vecR, list_pos_f), package = read_cell_and_pos_auto(
            dinput["folder_final_state"] + "/ratio-0.0000/scf"
        )
        if dinput["dimension"] == "xyz":
            vol = np.dot(np.cross(vecR[:, 0], vecR[:, 1]), vecR[:, 2])
            dim = 3
        elif dinput["dimension"] == "xy":
            vol = norm(np.cross(vecR[:, 0], vecR[:, 1]))
            dim = 2
        else:
            print("Unsupported dimension '%s'" % dinput["dimension"])
        list_cp = calc_cp_T(dim, vol, dinput["g"], data["Wif"], list_phonon_part)
        savedata_single_array("cp_T.txt", data, "filename_Cp_tilde", list_cp)

        filename = "cp_T.png"
        # plot_cp_T(filename, list_cp)

    if "filename_lifetime" in data:
        pass
    else:
        list_lifetime = calc_lifetime_T(dinput["g"], data["Wif"], list_phonon_part)
        savedata_single_array(
            "lifetime_T.txt", data, "filename_lifetime", list_lifetime
        )

    save()


if __name__ == "__main__":
    main()
