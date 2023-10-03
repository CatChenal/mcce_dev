#!/usr/bin/env python

"""
Charge microstate analysis tools
This library offers the following functions
 * Read microstate file and convert to a list of charge microstates
 * Group microstates by residue charge state
 * Group/rank microstates by population
 * Bin microstates
 * Find microstates within an energy band
 * Distance score of two microstate groups
"""


import sys
from dataclasses import dataclass
import numpy as np
import math
import tracemalloc
import zlib
import matplotlib.pyplot as plt


# the constants should be in their own module
ph2Kcal = 1.364
Kcal2kT = 1.688


class Microstate:
    """Class of a microstate.
    """
    def __init__(self, state, E, count):
        self._stateid = zlib.compress(" ".join([str(x) for x in state]).encode())
        self.E = E
        self.count = count

    def state(self):
        return [int(i) for i in zlib.decompress(self.stateid).decode().split()]


@dataclass
class Conformer:
    """Dataclass of a conformer.
    """
    iconf : int = 0
    ires : int = 0
    confid : str = ""
    resid : str = ""
    occ : float = 0.0
    crg : float = 0.0

    def load_from_head3_line(self, line):
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3]+self.confid[5:11]
        self.crg = float(fields[4])


@dataclass
class Conformer:
    """Dataclass of a conformer.
    """
    iconf : int = 0
    confid : str = ""
    resid : str = ""
    crg : float = 0.0

    def load(self, line):
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3]+self.confid[5:11]
        self.crg = float(fields[4])


# UNUSED
""" class Free_Res:
    def __init__(self):
        self.resid = ""
        self.charges = []     # a list of charge choice
 """
#UNUSED
""" def readheadlst(fname):
    conformers = []
    lines = open(fname).readlines()

    for line in lines[1:]:
        if len(line) > 80:
            conf = Conformer()
            conf.load(line)
            conformers.append(conf)

    return conformers
 """

class Charge_Microstate:
    def __init__(self, crg_stateid, total_E, count):
        self.crg_stateid = crg_stateid
        self.E = self.average_E = 0
        self.total_E = total_E
        self.count = count

    def state(self):
        return [int(i) for i in zlib.decompress(self.crg_stateid).decode().split()]

class Subset_Microstate:
    def __init__(self, subset_stateid, total_E, count):
        self.subset_stateid = subset_stateid
        self.E = self.average_E = 0
        self.total_E = total_E
        self.count = count

    def state(self):
        return [int(i) for i in zlib.decompress(self.subset_stateid).decode().split()]
    


class MSout:
    def __init__(self, fname):
        self.T = 273.15
        self.pH = 7.0
        self.Eh = 0.0
        self.N_ms = 0
        self.N_uniq = 0
        self.lowest_E = 0.0
        self.highest_E = 0.0
        self.average_E = 0.0
        self.fixed_iconfs = list()
        self.fixed_crg = 0.0
        self.fixed_ne = 0.0
        self.fixed_nh = 0.0
        self.free_residues = list()   # free residues, referred by conformer indices
        self.iconf2ires = dict()      # from conformer index to free residue index
        self.microstates = dict()
        self.conformers = list()

        self.load_msout(fname)


    @staticmethod
    def check_next_lines(lines: list):
        """Continue reading lines until non-comment line found.
           Return a data line and the remainding lines.
        """
        while True:
            line = lines.pop(0).strip()
            if line and line[0] != "#":
                break
        return line, lines


    def load_msout(self, fname: str):
        """Populate MSout variables from a microstate file.
        """
        lines = open(fname).readlines()

        # Get a valid line
        line, lines = self.check_next_lines(lines)

        fields = line.split(",")
        for field in fields:
            key, value = field.split(":")
            key = key.strip().upper()
            value = float(value)
            if key == "T":
                self.T = value
            elif key == "PH":
                self.pH = value
            elif key == "EH":
                self.Eh = value

        # second line, confirm this is from Monte Carlo sampling
        line, lines = self.check_next_lines(lines)

        key, value = line.split(":")
        if key.strip() != "METHOD" or value.strip() != "MONTERUNS":
            raise ValueError(f"This file {fname} is not a valid microstate file.")

        # Third line, fixed conformer indicies
        line, lines = self.check_next_lines(lines)

        self.fixed_iconfs = [int(i) for i in line.split(":")[1].split()]

        # 4th line, free residues
        line, lines = self.check_next_lines(lines)

        residues = line.split(":")[1].split(";")

        for r in residues:
            if r.strip():
                self.free_residues.append([int(i) for i in r.split()])

        for r, res in enumerate(self.free_residues):
            for iconf in res:
                self.iconf2ires[iconf] = r

        # find the next MC record
        found_mc = False
        new_mc = False

        for line in lines:
            if line.find("MC:") == 0: # ms starts
                found_mc = True
                new_mc = True
                continue
            elif new_mc:
                # current_state is a list of res indices:
                current_state = [int(c) for c in line.split(":")[1].split()]
                new_mc = False
                continue
            elif found_mc:
                fields = line.split(",")
                if len(fields) < 3:
                    continue  # or break?

                flipped = [int(c) for c in fields[2].split()]
                for ic in flipped:
                    current_state[self.iconf2ires[ic]] = ic

                ms = Microstate(current_state,
                                float(fields[0]),
                                int(fields[1]),
                                )
                # ms.state == current_state
                key = ",".join([f"{i}" for i in ms.state])
                if key in self.microstates:
                    self.microstates[key].count += ms.count
                else:
                    self.microstates[key] = ms

        # find N_ms, lowest, highest, average E
        E_sum = 0.0
        msvalues = self.microstates.values()
        self.N_uniq = len(msvalues)

        # Preset lowest/highest E to first one
        first_E = next(iter(msvalues)).E
        self.lowest_E = first_E
        self.highest_E = first_E

        for ms in msvalues:
            self.N_ms += ms.count
            E_sum += ms.E * ms.count

            if self.lowest_E > ms.E:
                self.lowest_E = ms.E

            if self.highest_E < ms.E:
                self.highest_E = ms.E

        self.average_E = E_sum / self.N_ms


# TODO:
# Define fn `ms_group_by` 
# with switch =["E_bounds","Conf_index","Conf_ID"]

def group_ms_by_energy(microstates:list, energy_bounds:list) -> list:
    """Divide the microstates into N=len(energy_bounds) bands using each energy in energy_bounds
    as lower boundary. 
 
    Args:
    microstates (list(MSout.microstates.values())): list of microstates objects.
    energy_bounds (list): list of energies to use as lower boundaries.
 
    Return:
    The list of N lists of energy bands sorted ascendingly.
    """

    N = len(energy_bounds)
    energy_bounds.sort()
    energy_bounds.append(1.0e100)    # add a big number as the righmost (?) boundary

    resulted_bands = [[] for i in range(N)]

    for ms in microstates:
        it = -1
        for n in range(N):
            if energy_bounds[n] <= ms.E < energy_bounds[n+1]:
                it = n
                break
        if it > -1:
            resulted_bands[it].append(ms)

    return resulted_bands


def group_ms_by_iconf(microstates:list, iconfs:list) -> tuple:
    """
    Split the `microstates` into two groups: the first contains one of the conformers indices
    in `iconfs`,  the second one contains none of the listed conformers.
    Args:
      microstates (list(MSout.microstates.values())): list of microstates objects.
      iconfs (list): list of conformer indices.

    Return:
    A tuple: group 1, group2.

    """
    ingroup = []
    outgroup = []
    for ms in microstates:
        contain = False
        for ic in iconfs:
            if ic in ms.state:
                ingroup.append(ms)
                contain = True
                break
        if not contain:
            outgroup.append(ms)

    return ingroup, outgroup

#REUME HERE: ADD conformers
def group_ms_by_confid(microstates:list, confids:list) -> tuple:
    """
    Group conformers in `microstates` by the conformer IDs in `confids`.
    An ID is considered a match as long as it is a substring of the conformer name.
    The first group contains all the matched microstates, the second group contains 
    all the unmatched ones.
    
    Args:
      microstates (list(MSout.microstates.values())): list of microstates objects.
      iconfs (list): list of conformer indices.

    Return:
      A tuple: group 1, group2.

    """
    ingroup = []
    outgroup = []
    for ms in microstates:
        contain = True
        # where is `conformers` comming from???
        names = [conformers[ic].confid for ic in ms.state]
        for confid in confids:
            in_names = False
            for name in names:
                if confid in name:
                    in_names = True
                    break
            contain = contain and in_names
        if contain:
            ingroup.append(ms)
        else:
            outgroup.append(ms)

    return ingroup, outgroup


def ms_energy_stat(microstates:list) -> tuple:
    """
    Return the lowest, average, and highest energy of the given macrostates collection.
    Args:
      microstates (list(MSout.microstates.values())).
    Return:
      A 3-tuple: lowest_E, average_E, highest_E.
    """

    lowest_E = highest_E = next(iter(microstates)).E
    N_ms = 0
    total_E = 0.0
    for ms in microstates:
        if lowest_E > ms.E:
            lowest_E = ms.E
        elif highest_E < ms.E:
            highest_E = ms.E
        N_ms += ms.count
        total_E += ms.E*ms.count

    average_E = total_E/N_ms

    return lowest_E, average_E, highest_E


def ms_convert2occ(microstates:list) -> dict:
    """
    Convert the conformers that appear at least once in the microstates to
    conformer occupancy.
    Args:
      microstates (list(MSout.microstates.values())).
    Return
      occ (dict): the conformer id is the key.
    """
    occurance = dict()  # occurance of conformer
    occ = dict()
    N_ms = 0
    for ms in microstates:
        N_ms += ms.count
        for ic in ms.state:
            if ic in occurance:
                occurance[ic] += ms.count
            else:
                occurance[ic] = ms.count

    for key in occurance.keys():
        occ[key] = occurance[key]/N_ms

    return occ


def ms_counts(microstates:list) -> int:
    """Return the sum total of microstates.count.
    Args:
      microstates (list(MSout.microstates.values())).
    Return:
      The total count.
    """
    return np.sum([ms.count for ms in microstates])


def ms_charge(ms:Microstate, conformers:list) -> float:
    """Sum the conformers' charge of the microstate ms.
    Args:
      ms (Microstate): a microsstate object.
      conformers (list): list of Conformer objects populated from head3.lst.
                         (Use read_conformers() to create the list).
    Return:
      The total charge (float).
    """
    #crg = 0.0
    #for ic in ms.state:
    #    crg += conformers[ic].crg
    return np.sum(conformers[ic].crg for ic in ms.state)


def ms_convert2sumcrg(microstates:list, free_res:list, conformers:list) -> list:
    """
    Given a list of microstates, convert to net charge of each free residue.
    Args:
      microstates (list(MSout.microstates.values())).
      free_res (list): MSout.free_residues list.
      conformers (list): list of Conformer objects populated from head3.lst.
                         (Use read_conformers() to create the list).
    Return:
      List of charges.
    """
    iconf2ires = dict()
    #for i_res in range(len(free_res)):
    #    for iconf in free_res[i_res]:
    #        iconf2ires[iconf] = i_res
    for r, fres in enumerate(free_res):
        for iconf in fres:
            iconf2ires[iconf] = r

    charges_total = np.zeros(len(free_res))
    N_ms = 0
    for ms in microstates:
        N_ms += ms.count
        for ic in ms.state:
            charges_total[iconf2ires[ic]] += conformers[ic].crg * ms.count

    return (charges_total/N_ms).tolist()


def read_conformers() -> list:
    """Return a list of Conformer objects popuated from the head3.lst file.
    Assumed: "head3.lst" file in working directory.
    """
    conformers = []
    lines = open("head3.lst").readlines()
    lines.pop(0)
    for line in lines:
        conf = Conformer()
        conf.load_from_head3_line(line)
        conformers.append(conf)

    return conformers


def e2occ(energies:list) -> float:
    """Return the Boltzmann distributed occupancy.
    Args:
      energies (list): list of energy values in Kcal/mol.
    Return:
      occupancy (float).  
    """
    e = np.array(energies)
    Pi_raw = np.exp(-Kcal2kT * (e - min(e)))

    return Pi_raw/np.sum(Pi_raw)


def bhata_distance(prob1:list, prob2:list) -> float:
    """
    Return the Bhattacharyya distance[1] of the distributions prob1 and prob2.
    Args:
      prob1, prob2 (lists):
    Return:
      The Bhattacharyya distance (float) if prob1 and prob2 have the same length,
      else 10,000.0.
    
    [1]: https://en.wikipedia.org/wiki/Bhattacharyya_distance
    """

    d_max = 10_000.0   # Max possible value set to this
    if len(prob1) != len(prob2):
        return d_max

    p1 = np.array(prob1) / np.sum(prob1)
    p2 = np.array(prob2) / np.sum(prob2)

    bc = np.sum(np.sqrt(p1 * p2))
    if bc <= np.exp(-d_max):
        d = d_max
    else:
        d = -np.log(bc)

    return d


def changed_res(msgroup1: list, msgroup2:list, free_res:list) -> list:
    """Return a list of Bhattacharyya distances of free residues occupancy in each group.
    Args:
      msgroup1 (list): list(MSout.microstates.values()).
      msgroup2 (list): list(MSout.microstates.values()).
      free_res (list): MSout.free_residues list.
    Return:
      List of floats.
    """
    occ1 = ms_convert2occ(msgroup1)
    occ2 = ms_convert2occ(msgroup2)

    bhd = list()
    for res in free_res:
        p1, p2 = list(), list()
        for ic in res:
            p1.append(occ1[ic] if ic in occ1 else 0.0)
            p2.append(occ2[ic] if ic in occ2 else 0.0)

        bhd.append(bhata_distance(p1, p2))

    return bhd


def changed_conf(msgroup1: list, msgroup2:list):
    """Given two group of microstates, calculate what changed at conformer level.

    Args:
      msgroup1 (list): list(MSout.microstates.values()).
      msgroup2 (list): list(MSout.microstates.values()).

    """
    occ1 = ms_convert2occ(msgroup1)
    occ2 = ms_convert2occ(msgroup2)

    all_keys = set(occ1.keys())
    all_keys |= set(occ2.keys())
    all_keys = list(all_keys)
    all_keys.sort()

    diff_occ = dict()
    for key in all_keys:
        p1 = occ1[key] if key in occ1 else 0.0
        p2 = occ2[key] if key in occ2 else 0.0

        diff_occ[key] = p2 - p1

    return diff_occ


def main(msout_filepath:str=None):
    """Main microstate analysis function available with module import.
    """

    conformers = read_conformers()

    file_used = "ms_out/pH4eH0ms.txt"
    if msout_filepath is not None:
        file_used = msout_filepath
        
    print(f"Microstate analysis of file: {file_used}\n\n")
    msout = MSout(file_used)

    # e_step = (msout.highest_E - msout.lowest_E)/20
    # ticks = [msout.lowest_E + e_step*(i) for i in range(20)]
    # ms_in_bands = group_ms_by_energy(msout.microstates.values(), ticks)
    # print([len(band) for band in ms_in_bands])
    # netural, charged = group_ms_by_iconf(msout.microstates.values(), [12, 13, 14, 15])
    # l_E, a_E, h_E = ms_energy_stat(msout.microstates.values())
    # print(l_E, a_E, h_E)

    # charge over energy bands

    # e_step = (msout.highest_E - msout.lowest_E) / 20
    # ticks = [msout.lowest_E + e_step*(i+1) for i in range(19)]
    # ms_in_bands = group_ms_by_energy(msout.microstates.values(), ticks)
    # for band in ms_in_bands:
    #     band_total_crg = 0.0
    #     for ms in band:
    #         band_total_crg += ms_charge(ms)
    #     print(band_total_crg/ms_counts(band))

    # netural, charged = group_ms_by_iconf(msout.microstates.values(), [12, 13, 14, 15])
    # diff_occ = changed_conf(netural, charged)
    # for key in diff_occ.keys():
    #     print("%3d, %s: %6.3f" % (key, conformers[key].confid, diff_occ[key]))

    # diff_bhd = changed_res(netural, charged, msout.free_residues)
    # for ir in range(len(msout.free_residues)):
    #     print("%s: %6.4f" % (conformers[msout.free_residues[ir][0]].resid, diff_bhd[ir]))
    # charges = ms_convert2sumcrg(msout.microstates.values(), msout.free_residues)
    # for ir in range(len(msout.free_residues)):
    #     print("%s: %6.4f" % (conformers[msout.free_residues[ir][0]].resid, charges[ir]))

    microstates = list(msout.microstates.values())

    Res_OI = ["GLU-1A0035"]  # residues of interest
    charged_res, _ = group_ms_by_confid(microstates, Res_OI)
    print(f"Number of microstates: {len(microstates)}\n",
          f"Number of confomers for res {Res_OI}: {len(charged_res)}\n")


if __name__ == "__main__":

    main()
