#!/usr/bin/env python

from __future__ import print_function

from io_package import read_cell_and_pos_auto
from libplot import plot_tot_Q, plot_eig_Q, plot_overlap_Q, plot_cp_T
from libcalc import calc_dQ, calc_freq, calc_wif, calc_phonon_part, calc_cp_T, \
    calc_phonon_part_T0_HR, calc_lifetime_T, calc_dE
from libreadqe import read_pos_and_etot_ratio
import sys
import os
import yaml
import numpy as np
from numpy.linalg import norm
from datetime import datetime

from constant import indent

np.seterr(all="log")


def main():
    '''
    Main entrance to nonradiative code
    '''

    def get_datetime():
        return datetime.now().strftime("%m/%d/%Y at %H:%M:%S")

    def intro():
        print("Beginning nonradiative code on {}\n".format(get_datetime()))

    def outro():
        print("Ending nonradiative code on {}\n".format(get_datetime()))

    def display_input():
        print('Input parsed as:')
        for key in dinput:
            print('{}{} : {}'.format(indent, key, dinput[key]))
        print()

    # subfunctions in main
    def saveinput():
        with open(file_input, 'w') as f:
            yaml.dump(dinput, f, default_flow_style=False, width=200)

    def save():
        with open(file_store, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, width=200)

    def save_dic_ndarray(filename, dic):
        '''
        Same a dictionary with named ndarray
        '''
        with open(filename, 'w') as f:
            for name, ar in dic.items():
                f.write("# %s\n" % name)
                np.savetxt(f, ar)
                f.write("\n")

    def savedata_single_array(filename, database, term, ar):
        '''
        Save an array to a file and record the filename in the database
        '''
        np.savetxt(filename, ar)
        database[term] = filename
        # Call save all after
        save()

    def error_with_input():
        print("Please modify %s according to instruction inside" % file_input)
        saveinput()
        sys.exit(1)

    # begin main code
    intro()
    file_input = "input.yaml"

    # read input file
    if (not os.path.exists(file_input)):
        raise FileNotFoundError('Missing input file: input.yaml')
    else:
        with open(file_input, 'r') as f:
            dinput = yaml.load(f.read(), Loader=yaml.FullLoader)

    # report parse input
    display_input()

    # check to see what job is running
    try:
        dinput['job'] = dinput['job'].lower()
    except KeyError:
        raise ValueError('job type must be specified in input file: nonrad or isc')
    if dinput['job'] == 'nonrad':
        print("Calculating nonradiative recombination rate")
    elif dinput['job'] == 'isc':
        print("Calculating intersystem crossing rate")
    else:
        raise ValueError('job type not recognized: {}'.format(dinput['job']))

    # name of file where data will be written to
    file_store = "{}.yaml".format(dinput['job'])

    # check to see if data file already exists
    if (not os.path.exists(file_store)):
        data = {}
    else:
        with open(file_store, 'r') as f:
            data = yaml.load(f.read(), Loader=yaml.FullLoader)

    # specific to nonrad checks
    if dinput['job'] == 'nonrad':
        # check band index provided is of an okay format
        if (not isinstance(dinput['bulkband_index'], list) or not isinstance(dinput['defectband_index'], int)):
            print("Error: Encountered when checking provided bulkband_index (must be a list of integers)")
            error_with_input()
        # check that dimension was provided
        if "The " in dinput['dimension']:
            print("Error: Encountered when checking provided dimension")
            error_with_input()

    # check that type of dE is correct
    if not isinstance(dinput['dE'], float):
        print("Error: Encountered when checking provided dE (must be a float)")
        error_with_input()

    # Set default input values
    if ("enable_plot" not in dinput):
        dinput['enable_plot'] = False

    if ("defectband_spin" not in dinput):
        dinput['defectband_spin'] = "down"  # This is for calculation as PRB 90,075202

    if ("unit" not in data):
        data['unit'] = {
            "dQ": "aMU^1/2 Ang",
            "M": "aMU",
            "hbarFreqi": "meV",
            "hbarFreqi_error": "meV",
            "hbarFreqf": "meV",
            "hbarFreqf_error": "meV",
        }

    # check for lambda_t(ransverse)
    if (dinput['job'] == 'isc'):
        if ('lambda_t' not in dinput):
            data['lambda_t'] = 1.0
        else:
            data['lambda_t'] = dinput['lambda_t']

    # Compute dQ
    if ("dQ" not in data and "folder_init_state" in dinput):
        print("\n{}Compute dQ".format(indent))
        (vecR, list_pos_f), package = read_cell_and_pos_auto(dinput['folder_final_state'] + "/ratio-0.0000/scf")
        (vecR, list_pos_i), package = read_cell_and_pos_auto(dinput['folder_final_state'] + "/ratio-1.0000/scf")
        dQ, M, ar_delta, ar_Qi = calc_dQ(vecR, list_pos_f, list_pos_i)
        data['dQ'] = dQ
        data['M'] = M
        print("{}dQ = {}".format(indent*2, dQ))
        print("{}M = {}".format(indent*2, M))
    else:
        print("\n{}Read dQ from {}".format(indent, file_store))

    # Plot Q-etot curve
    if ("filename_Q_etot" not in data):
        print("\n{}Generating Q-etot curves".format(indent))
        filename = "tot_Q.png"
        data['filename_Q_etot'] = filename
        list_data_i, vecR = read_pos_and_etot_ratio(dinput['folder_init_state'])
        list_data_f, vecR = read_pos_and_etot_ratio(dinput['folder_final_state'])
        ar_i, ar_f = plot_tot_Q(filename, list_data_i, list_data_f, dinput['dE'], data['dQ'], dinput['enable_plot'])

        # Save Q-etot data
        print("{}Saving Q-etot data to {}".format(indent*2, "tot_Q_i.txt"))
        savedata_single_array("tot_Q_i.txt", data, "filename_Q_etot_i", ar_i)
        print("{}Saving Q-etot data to {}".format(indent*2, "tot_Q_f.txt"))
        savedata_single_array("tot_Q_f.txt", data, "filename_Q_etot_f", ar_f)
    else:
        print("\n{}Read Q-etot from {}".format(indent, file_store))

    if ("dErelf" not in data):
        print("\n{}Calculating dE".format(indent))
        dEf = calc_dE(dinput['folder_final_state'])
        dEi = calc_dE(dinput['folder_init_state'])
        data['dErelf'] = dEf
        data['dEreli'] = dEi
        print("{}dErelf = {}".format(indent*2, data['dErelf']))
        print("{}dEreli = {}".format(indent*2, data['dEreli']))
        save()
    else:
        print("\n{}Read dE from {}".format(indent, file_store))

    # Compute frequency
    if ("hbarFreqi" not in data):
        print("\n{}Computing frequency".format(indent))
        # From 0-0.15 and 0-0.20

        # computing initial frequency
        r1 = dinput['ratio_init_min']
        r2 = dinput['ratio_init_max']
        freq, freq_error, list_result, S = calc_freq(
            dinput['folder_init_state'], ratio_min=r1, ratio_max=r2, dQ=data['dQ'])
        data['hbarFreqi'] = freq * 1000
        data['hbarFreqi_error'] = freq_error * 1000
        data['hbarFreqi_order'] = list_result
        data['S_i'] = S
        print("{}hbarFreqi = {}".format(indent*3, data['hbarFreqi']))
        print("{}S_i = {}".format(indent*3, data['S_i']))

        # computing final frequency
        r1 = dinput['ratio_final_min']
        r2 = dinput['ratio_final_max']
        freq, freq_error, list_result, S = calc_freq(
            dinput['folder_final_state'], ratio_min=r1, ratio_max=r2, dQ=data['dQ'])
        data['hbarFreqf'] = freq * 1000
        data['hbarFreqf_error'] = freq_error * 1000
        data['hbarFreqf_order'] = list_result
        data['S_f'] = S
        print("{}hbarFreqf = {}".format(indent*3, data['hbarFreqf']))
        print("{}S_f = {}".format(indent*3, data['S_f']))
        save()
    else:
        print("\n{}Read frequencies from {}".format(indent, file_store))

    # if ("hbarFreqi_to02_order" not in data):
    #     print("\n{}Computing frequency to higher orders".format(indent))
    #     r1 = dinput['ratio_init']
    #     r2 = abs(r1 - 0.20)
    #     if (r1 > r2):
    #         r1, r2 = (r2, r1)
    #     freq, freq_error, list_result, S = calc_freq(
    #         dinput['folder_init_state'], ratio_min=r1, ratio_max=r2, dQ=data['dQ'])
    #     data['hbarFreqi_to02_order'] = list_result

    #     r1 = dinput['ratio_final']
    #     r2 = abs(r1 - 0.20)
    #     if (r1 > r2):
    #         r1, r2 = (r2, r1)
    #     freq, freq_error, list_result, S = calc_freq(
    #         dinput['folder_final_state'], ratio_min=r1, ratio_max=r2, dQ=data['dQ'])
    #     data['hbarFreqf_to02_order'] = list_result
    #     # print("{}hbarFreq_to02_order (i,f) = ({}, {})".format(
    #     #     indent*2, data['hbarFreqi_to02_order'], data['hbarFreqf_to02_order']))
    #     save()

    if dinput['job'] == 'nonrad':
        if ("Wif" not in data or "WifBand" not in data):
            print("\n{}Calculating electron wavefunction overlap Wif".format(indent))
            # Calculate Wif
            # Plot and store data
            dic_eig_occ, diffEigQ0, dic_band_overlap, dic_wif, ixband_wifmax, wif = calc_wif(
                dinput['folder_init_state'], dinput['folder_final_state'], dinput['defectband_index'],
                dinput['bulkband_index'][0],
                dinput['bulkband_index'][1],
                data['dQ'],
                de=None,
                spinname=dinput['defectband_spin'])

            for (statename, spin), val in dic_eig_occ.items():
                casename = "%s-%s" % (statename, "up" if spin == 1 else "down")
                filename = "eig-%s.txt" % casename
                np.savetxt(filename, val["eig"])
                filename = "occ-%s.txt" % casename
                np.savetxt(filename, val["occ"])
                data['filename_eigvals_%s' % casename] = filename
                if dinput['enable_plot']:
                    plot_eig_Q("eig-%s.png" % casename, val["eig"], val["occ"], data['dQ'])

            filename = "overlap.txt"
            print("\n{}Saving overlap data to {}".format(indent*2, filename))
            data['filename_overlap'] = filename
            save_dic_ndarray(filename, dic_band_overlap)

            if dinput['enable_plot']:
                plot_overlap_Q("overlap.png", dic_band_overlap, dic_wif, data['dQ'], diffEigQ0)

            data['diffEigQ0'] = diffEigQ0
            data['Wif'] = wif
            data['WifBand'] = ixband_wifmax
            data['Wif_Allband'] = dic_wif
            print("{}diffEigQ0 = {}".format(indent*3, data['diffEigQ0']))
            print("{}Wif = {}".format(indent*3, data['Wif']))
            print("{}WifBand = {}".format(indent*3, data['WifBand']))
            print("{}Wif_Allband = {}".format(indent*3, data['Wif_Allband']))
            save()
        else:
            print("\n{}Read Wif from {}".format(indent, file_store))

    if ("filename_phonon_part" not in data):
        if dinput['job'] == 'nonrad':
            data['order_x'] = 1
        elif dinput['job'] == 'isc':
            data['order_x'] = 0
        print("\n{}Calculating phonon wavefunction overlap Xif".format(indent))
        list_phonon_part = calc_phonon_part(
            dinput['dE'],
            data['dQ'],
            data['hbarFreqi'],
            data['hbarFreqf'],
            dinput['temperature'],
            data['order_x']
        )

        savedata_single_array("phonon_part.txt", data, "filename_phonon_part", list_phonon_part)
    else:
        print("\n{}Read Xif from {}".format(indent, file_store))
        list_phonon_part = np.loadtxt("phonon_part.txt")

    # Calculate the approximate phonon part at T=0 based on HR theory
    if ("phonon_part_T0_HR" not in data):
        print("\n{}Calculating HR phonon part at T=0".format(indent))
        data['phonon_part_T0_HR'] = calc_phonon_part_T0_HR(
            dinput['dE'],
            data['dQ'],
            data['hbarFreqf'],
            data['order_x']
        )

    if dinput['job'] == 'nonrad':
        # Calculate Cp
        if ("filename_Cp_tilde" not in data):
            print("\n{}Calculating capture coefficient Cp".format(indent))
            # Read and compute volume
            (vecR, list_pos_f), package = read_cell_and_pos_auto(dinput['folder_final_state'] + "/ratio-0.0000/scf")
            if (dinput['dimension'] == "xyz"):
                vol = np.dot(np.cross(vecR[:, 0], vecR[:, 1]), vecR[:, 2])
                dim = 3
            elif (dinput['dimension'] == "xy"):
                vol = norm(np.cross(vecR[:, 0], vecR[:, 1]))
                dim = 2
            else:
                print("Unsupported dimension '%s'" % dinput['dimension'])

            # if dinput['job'] == 'nonrad':
            list_cp = calc_cp_T(dim, vol, dinput['g'], data['Wif'], list_phonon_part)
            # elif dinput['job'] == 'isc':
            #     list_cp = calc_cp_T(dim, vol, dinput['g'], data['lambda_t'], list_phonon_part)

            filename = "cp_T.txt"
            print('{}Saving cp_T to {}'.format(indent*2, filename))
            savedata_single_array(filename, data, "filename_Cp_tilde", list_cp)

            if dinput['enable_plot']:
                filename = "cp_T.png"
                plot_cp_T(filename, list_cp)
        else:
            filename = "cp_T.txt"
            print("\n{}Read Cp from {}".format(indent, filename))
            list_cp = np.loadtxt(filename)

    if ("filename_lifetime" not in data):
        print("\n{}Calculating lifetime".format(indent))

        if dinput['job'] == 'nonrad':
            list_lifetime = calc_lifetime_T(dinput['g'], data['Wif'], list_phonon_part, dinput['job'])
        elif dinput['job'] == 'isc':
            list_lifetime = calc_lifetime_T(dinput['g'], data['lambda_t'], list_phonon_part, dinput['job'])

        print('%s%s   %s' % (indent*2, 'T(K)', 'lifetime(s)'))
        for T, lif in list_lifetime:
            print('%s%4d  %.12e' % (indent*2, T, lif))
        # print('\n', list_lifetime, '\n')
        filename = "lifetime_T.txt"
        print('{}Saving lifetime_T to {}'.format(indent*2, filename))
        savedata_single_array(filename, data, "filename_lifetime", list_lifetime)
    else:
        print("\n{}Umm ... nothing to do? (delete {} if you want to run a calculation from scratch)".format(
            indent, file_store))

    save()

    print('\nDone! :)\n')
    outro()


if __name__ == "__main__":
    main()
