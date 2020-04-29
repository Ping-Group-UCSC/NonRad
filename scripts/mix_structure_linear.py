#!/usr/bin/env python

import os
import numpy as np
from string import digits, ascii_letters
import re
import sys
import shutil

sys.path.append("../code")

from io_package import read_cell_and_pos_auto, write_cell_and_pos_auto


def mix_structure(vecR, list_pos1, list_pos2, ratio, head, tail, filename_out):
    '''
    Mix two structures
    If ratio is 0, as pos1, 1 as pos2
    Otherwise use linear extrapolation

    Also check atom names (except number suffix, which can be different as polaron positions)

    the starting magnetization will also be changed according to ratio, but will never be zero

    there must be one polaron in each list_pos, and they must not be the same atom
    After mixing the polaron in list_pos1 will be named as "element, suffix_polaron+1" and in list_pos2 will be named as "element, suffix_polaron+2"

    :param suffix_polaron: the digit suffix to indicate where is the polaron
    :param ix_polaron_species: the 1-based index of species of the first polaron atom
    '''
    if (len(list_pos1) != len(list_pos2)):
        raise ValueError("Different number of atoms")


    l3 = []
    suffix1 = None
    suffix2 = None
    for p1, p2 in zip(list_pos1, list_pos2):
        element1 = p1["species"].translate(None, digits)
        element2 = p2["species"].translate(None, digits)
        if (element1 != element2):
            raise ValueError("Different atom found : %s %s" % (p1["species"], p2["species"]))
#Match atom positions to the nearest position (0.01 and 0.99 => 0.01 and -0.01)
        vec2 = p2["pos"].copy()
        for i in range(len(vec2)):
            if (abs(vec2[i] - p1["pos"][i]) > 0.5):
                vec2[i] -= round(vec2[i] - p1["pos"][i], 0)

        pos3 = p1["pos"] * (1-ratio) + vec2 * ratio
        name = p1["species"]
        l3.append({"species": name, "pos" : pos3})

    with open(filename_out, 'w') as f:
        for line in head:
            f.write(line)
        for atom in l3:
            f.write("%-4s %.10f %.10f %.10f\n" % (atom["species"], atom["pos"][0], atom["pos"][1], atom["pos"][2]))
        for line in tail:
            f.write(line)

    return

filename_scf = "template.in" 
with open(filename_scf, 'r') as f:
   head = f.readlines() 

filename_occ = "occ.in" 
if (os.path.exists(filename_occ)):
    with open(filename_occ, 'r') as f:
       tail = f.readlines() 
else:
    tail = []

(vecR, list_pos1), prog = read_cell_and_pos_auto("../relax-gs/relax")
(vecR2, list_pos2), prog = read_cell_and_pos_auto("../relax-cdftup1/relax")
if ((vecR != vecR2).any()):
    print(vecR, vecR2)
    raise ValueError("Inconsistent lattice vectors!")
for ratio in (0, 0.05, 0.10, 0.15, 0.2, 0.25, 0.5, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00):
    dirname = "ratio-%.4f" % ratio
    if (not os.path.exists(dirname)):
        os.makedirs(dirname)
    if (os.path.exists(os.path.join(dirname, "scf.in"))):
        continue
    mix_structure(vecR, list_pos1, list_pos2, ratio, head, tail, os.path.join(dirname, "scf.in"))

