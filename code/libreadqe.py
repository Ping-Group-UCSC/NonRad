#!/usr/bin/env python
import os
import glob
import re
import xml.etree.ElementTree as ET
import mmap
import numpy as np

from io_package import read_cell_and_pos_auto

Ha2eV = 27.2113834
Ry2eV = Ha2eV / 2


def get_ratio_folder(ratio):
    '''
    Get folder name of given ratio
    '''
    return "ratio-%.4f" % ratio


def get_save_folder(folder):
    '''
    Get the *.save folder under a given directory, assuming only one .save
    If not found ,return None
    '''
    filename = None
    for filename in glob.glob("%s/*.save" % folder):
        break
    if filename is None:
        # also check temp which is a common outdir name
        for filename in glob.glob("%s/temp/*.save" % folder):
            break
    return filename


def read_pos_and_etot_ratio(folder=None):
    '''
    Read ratio-***/scf.out from given directory

    Return a list of ratio/atoms/etot, and vecR (which should be the same for all structures)
    '''
    if (folder is None):
        folder = os.getcwd()

    # check for folders/files
    if not glob.glob('%s/ratio-*' % folder):
        raise FileNotFoundError("folders %s/ratio-* do not exist" % folder)
    elif not glob.glob("%s/ratio-*/scf.out" % folder):
        raise FileNotFoundError("files %s/ratio-*/scf.out do not exist" % folder)

    list_data = []
    for filename in glob.glob("%s/ratio-*/scf.out" % folder):
        ratio = float(re.match(".*ratio-(.+)/scf.out", filename).group(1))
# Read structure
        (vecR, list_pos), prog = read_cell_and_pos_auto(filename.replace(".out", ""))
# Remove numbers in species and convert to cartesian
        for atom in list_pos:
            atom["species"] = filter(lambda x: x.isalpha(), atom["species"])
            atom["posxyz"] = np.dot(vecR, atom["pos"])

# Read energy
# Always read first \!\! for hybrid or \! for non-yhybrid
        with open(filename.replace(".in", ".out"), 'r') as f:
            lines = f.readlines()
            for line in lines[::-1]:
                if ("total energy" in line and "is the sum" not in line):
                    tag = line[:2].strip()
                    break

            if (tag != "!" and tag != "!!"):
                print("Cannot recognize total energy in %s" % filename.replace(".in", ".out"))
                continue

            for line in lines:
                if (line.startswith(tag)):
                    etot = float(line.split()[-2]) * Ry2eV
                    break
        list_data.append({"ratio": ratio, "pos": list_pos, "etot": etot})

    list_data.sort(key=lambda x: x["ratio"])
    return list_data, vecR


def read_eig_61(folder):
    '''
    Read eigenvalues from save folder in qe-6.1

    :return: Array in ik, ispin, ib, 1/2 (for eigenvalues and occupations numbers)
    '''
    t1 = ET.parse(os.path.join(folder, "data-file.xml")).getroot()
    nk = int(t1.find("BRILLOUIN_ZONE/NUMBER_OF_K-POINTS").text)
    list_spin = [".1", ".2"] if t1.find(
        "SPIN/LSDA").text.strip() == "T" else [""]
    # list_data = []
    ar = None
    for ik in range(1, nk+1):
        for ispin, stspin in enumerate(list_spin):
            filename = t1.find("EIGENVALUES/K-POINT.%i/DATAFILE%s" %
                               (ik, stspin)).attrib["iotk_link"]
            t2 = ET.parse(os.path.join(folder, filename)).getroot()
            data = np.asarray([float(x)
                               for x in t2.find("EIGENVALUES").text.split()])
            occ = np.asarray([float(x)
                              for x in t2.find("OCCUPATIONS").text.split()])
            unit = t2.find("UNITS_FOR_ENERGIES").attrib["UNITS"]
            if (unit == "Hartree"):
                data = data * Ha2eV
            else:
                raise ValueError("Unkown unit")
            if (ar is None):
                ar = np.zeros((nk, len(list_spin), len(data), 2),
                              dtype=np.float64)
            ar[ik-1, ispin, :, 0] = data
            ar[ik-1, ispin, :, 1] = occ

    return ar


def parse_attribute(st):
    try:
        # st is bytes object but this works in python2
        l1 = [x.split("=") for x in st.split()[1:]]
    except TypeError:
        # raises a TypeError in python3, need to decode st
        l1 = [x.split("=") for x in st.decode().split()[1:]]

    for x in l1:
        x[1] = x[1][1:-1]
        if (x[1].isdigit()):
            x[1] = int(x[1])
    return dict(l1)


def read_wave_61(prefix, ispin, ik, ib):
    '''
    Read wavefunction from qe-6.1
    '''
    folder = os.path.join(prefix, "K%05i" % ik)
    if (not os.path.exists(folder)):
        raise ValueError("Cannot find folder %s" % folder)
    filename = os.path.join(folder, "evc.dat")
    if (not os.path.exists(filename)):
        filename = os.path.join(folder, "evc%i.dat" % ispin)
        if (not os.path.exists(filename)):
            raise ValueError("Cannot find file %s" % filename)

    evc = None
    with open(filename, 'rb') as f:
        mm = mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_READ)
        ix = mm.find(b"gamma_only")
        mm.seek(ix + 12)
        gamma_only = mm.read(1) == 'T'
        ix = mm.find(b"<evc.%i" % ib)
        if (ix == -1):
            mm.close()
            print("Reading file %s" % filename)
            raise ValueError("Cannot find band %i" % ib)
        # ix2 = mm.find(b">", ix)
        # ix3 = mm.find(b"</evc.%i" % ib, ix2)
        # print("Block: %i %i" % (ix2, ix3))

        mm.seek(ix)
        st = mm.readline().strip()[1:-1]
        dic_attr = parse_attribute(st)
        if (dic_attr["type"] == "complex"):
            sizev = dic_attr["kind"] * 2
        # Skip trailing marker of previous record
        # Leading marker of this record
        # Dummy 0 for iotk not raw (which is true for evc*.dat)
        header_record = np.fromstring(mm.read(4*3), dtype=np.int32)
#       print(header_record[1])
        if (header_record[2] != 0):
            print("Dummy not zero, wrong format")
#       print(mm.tell())
        if (header_record[1] != dic_attr["size"] * sizev + 4):
            print("Record size not consistent with data length")

        data = mm.read(dic_attr["size"] * sizev)
        data = np.fromstring(data, dtype='<f%s' % dic_attr["kind"])
#       print(mm.tell())

#       print(data[:5])
        trailing_record = np.fromstring(mm.read(4), dtype=np.int32)
        if (trailing_record[0] != header_record[1]):
            print("Leading and trailing marker not consistent")

        mm.close()

        if (dic_attr["type"] == "complex" and dic_attr["kind"] == 8):
            dt = np.complex128
        elif (dic_attr["type"] == "real" and dic_attr["kind"] == 8):
            dt = np.float64
        else:
            raise ValueError("Unknown data type %s" % dic_attr)
        if (gamma_only):  # Real for gamma-only
            dt = np.float64

        evc = data.view(dtype=dt)

    return evc


def read_wave_66(prefix, ispin, ik, ib):
    '''
    Read wavefunction from qe-6.6 (dat format, not hdf5)

    WARNING!!! This should only be used for reading kp=(0, 0, 0), others are not tested

    WARNING!!! gamma_only is not implemented (should read real wfc instead of complex)

    WARNING!!! npol != 1 is not implemented (e.g. noncollinear case)
    '''
    sspin = 'up' if ispin == 1 else 'dw'
    filename = os.path.join(prefix, f'wfc{sspin}{ik}.dat')
    with open(filename, 'rb') as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buffer:
            # ------------------------------------------------------------
            buffer.seek(56)
            ngw, igwx_, npol, nbnd_ = \
                np.frombuffer(buffer.read(16), dtype=np.int32)
            # ------------------------------------------------------------
            gksize = 12 * igwx_  # 3*int32 = 12 bytes
            seek = 56 + 16 + 8 + 80 + gksize + 8
            wfcsize = 16 * igwx_  # complex128 = 16 bytes
            seek += wfcsize * ib
            buffer.seek(seek)
            # ------------------------------------------------------------
            return np.frombuffer(buffer.read(wfcsize), dtype=np.complex128)


def read_eig_66(folder):
    '''
    Read eigenvalues from save folder in qe-6.6

    :return: Array in ik, ispin, ib, 1/2 (for eigenvalues and occupations numbers)
    '''
    root = ET.parse(os.path.join(folder, "data-file-schema.xml")).getroot()
    band_child = root.find('output').find('band_structure')
    lsda = band_child.find('lsda').text == 'true'
    nspin = 2 if lsda else 1
    bnd_str = 'nbnd_up' if lsda else 'nbnd'
    nbnd = int(band_child.find(bnd_str).text)
    nk = int(band_child.find('nks').text)
    eig, occ = [], []
    for ks_child in band_child.findall('ks_energies'):
        keig = np.fromstring(ks_child.find('eigenvalues').text, sep=' ') * Ha2eV
        kocc = np.fromstring(ks_child.find('occupations').text, sep=' ')
        eig.append(keig.reshape(nspin, nbnd))
        occ.append(kocc.reshape(nspin, nbnd))
    ar = np.zeros((nk, nspin, nbnd, 2), dtype=np.float64)
    ar[:, :, :, 0] = eig
    ar[:, :, :, 1] = occ
    return ar


def read_eig(folder):
    '''
    check for 6.1 or 6.6 style folder than call read_eig_61 or read_eig_66
    see read_eig_61 and read_eig_66 for more details
    '''
    file61 = os.path.join(folder, 'data-file.xml')
    file66 = os.path.join(folder, 'data-file-schema.xml')
    if os.path.exists(file61):
        return read_eig_61(folder)
    elif os.path.exists(file66):
        return read_eig_66(folder)
    else:
        raise ValueError(
            f"Unable to locate an xml data file in folder: {folder}")


def read_wave(prefix, ispin, ik, ib):
    '''
    check for 6.1 or 6.6 style folder than call read_wave_61 or read_wave_66
    see read_wave_61 and read_wave_66 for more details
    '''
    file61_nospin = os.path.join(prefix, "K%05i" % ik, "evc.dat")
    file61 = os.path.join(prefix, "K%05i" % ik, "evc%i.dat" % ispin)
    file66_nospin = os.path.join(prefix, f'wfc{ik}.dat')
    sspin = 'up' if ispin == 1 else 'dw'
    file66 = os.path.join(prefix, f'wfc{sspin}{ik}.dat')
    if os.path.exists(file61_nospin) or os.path.exists(file61):
        return read_wave_61(prefix, ispin, ik, ib)
    elif os.path.exists(file66_nospin) or os.path.exists(file66):
        return read_wave_66(prefix, ispin, ik, ib)
    else:
        raise ValueError(
            f"Unable to locate a wfc file in folder: {prefix}")
