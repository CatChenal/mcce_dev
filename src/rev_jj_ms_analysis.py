#!/usr/bin/env python

"""
Charge microstate analysis tools
This library offers the following functions:

 * Read microstate file and convert to a list of charge microstates
 * Group microstates by residue charge state
 * Group/rank microstates by population
 * Bin microstates
 * Find microstates within an energy band
 * Distance score of two microstate groups
"""

import sys
import zlib
import math
import tracemalloc
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt


# the constants should be in their own module
ph2Kcal = 1.364
Kcal2kT = 1.688

@dataclass
class Conformer:
    """Class of a conformer.
    """
    iconf : int = 0
    confid : str = ""
    resid : str = ""
    crg : float = 0.0

   # load_from_head3_line
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
# UNUSED
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

class Microstate:
    """Base class for building a microstate object.
    """
    def __init__(self, state: list, E: float, count: int):
        self._stateid = zlib.compress(" ".join([str(x) for x in state]).encode())
        self.E = E
        self.count = count
        self.average_E = 0.0
        self.total_E = 0.0
        
    def state(self):
        return [int(i) for i in zlib.decompress(self._stateid).decode().split()]


# Not needed:
class Charge_Microstate(Microstate):
    def __init__(self, crg_stateid, total_E, count):
        super().__init__(crg_stateid, total_E, count)
        
        self._crg_stateid = super()._stateid # rename
        self.E = 0.0
        self.total_E = total_E
        

class Subset_Microstate(Microstate):
    def __init__(self, subset_stateid, total_E, count):
        super().__init__(subset_stateid, total_E, count)

        self._subset_stateid = super()._stateid # rename
        self.E = 0.0
        self.total_E = total_E


class MC:
    """
    Class for a microstate.

    Vars:
    Methods:
    """
    def __init__(self):
        self.T = 298.15
        self.pH = 7.0
        self.Eh = 0.0

        self.N_ms = 0
        self.N_uniq = 0

        self.lowest_E = 0.0
        self.highest_E = 0.0
        self.average_E = 0.0
        
        self.fixed_crg = 0.0
        self.fixed_ne = 0.0
        self.fixed_nh = 0.0
        self.fixed_confs = list()      # fixed conformers

        self.free_residues = list()       # list of conformer groups constituents of free residues
        self.free_residue_names = list()  # free residue names
        self.iconf2ires = dict()          # from conformer index to free residue index
        self.ires_by_iconf = dict()       # index of free residue by index of conf
        self.iconf_by_confname = dict()
    
        self.microstates = list()        # a list of microstates
        self.microstates_by_id = dict()
        self.conformers = list()

        self.method = ""
        self.counts = 0

        #TODO: trap open error
        # populate self.conformers and self.iconf_by_confname:
        lines = open("head3.lst").readlines()
        iconf = 0
        for line in lines[1:]:
            if len(line) > 80:
                conf = Conformer()
                conf.load(line)
                self.conformers.append(conf)
                self.iconf_by_confname[conf.confid] = iconf
                iconf += 1


    @staticmethod
    def fields_from_line(line: str):
        """Return the line and the fields extracted from it.
        Used by MC.read_ms.
        """
        line = line.split("#")[0]  #remove comment
        fields = line.split(",")

        return line, fields
    

    def read_ms(self, fname: str, mc_numbers: list=None):
        """MC class method. Read a file in ms_out folder (a.k.a "a ticker tape" file), 
        for example, "/ms_out/pH7eH0.txt". All runs of Monte Carlo sampling (currently 6) 
        will be included with default `mc_numbers` (None), else the one or more runs of 
        interest must be passed as a list.  These numbers correspond to the `MC:n` line 
        in the microstate file.

        Args:
          fname (str): the name of microstate file as found in mc_out folder.
          mc_numbers (list(int)): list of MC runs to consider.
        """
        
        if mc_numbers:
            MC_Segments = mc_numbers
        else:
            MC_Segments = [0,1,2,3,4,5]
            # Alternative: 
            # Set selection to single mc run with a warning message.
            # Purpose: 'entice' users to use a subset:
            #
            # MC_Segments = [0,]
            # msg = "You have not selected any specific MC run (`mc_number` arg is None),\n" 
            # msg = msg + "so only the first one will be processed.\n"
            # msg = msg + "The call to process all MC runs is `.readms(msfile, "
            # msg = msg + "mc_numbers=[0,1,2,3,4,5])\n"
            # msg = msg + "Enter Y/y if that's ok to continue (press anything to redo): "
            # Ask for redo:
            # reply = input(msg).lower()[0]
            # if reply != 'y':
            #    sys.exit(0)
            #

        # Only start variables assignment after a possible exit point (if possible).
        self.microstates = list()   # reset
        self.counts = 0             # reset
        fields = list()
        n_lines = 0

        with open(fname) as f:
            while len(fields) != 3 and n_lines < 10:
                _, fields = self.fields_from_line(f.readline())
                n_lines += 1

            if n_lines >= 10:
                msg = "Not found: A condition line like 'T:298.15,pH:5.00,eH:0.00' \
                    in the first 10 lines."
                raise ValueError(msg)

            for field in fields:
                key, value = field.split(":")
                key = key.strip()
                value = float(value.strip())
                if key.upper() == "T":
                    self.T = value
                elif key.upper() == "PH":
                    self.pH = value
                elif key.upper() == "EH":
                    self.Eh = value

            # method record
            fields = list()  # reset
            while len(fields) != 2:
                _, fields = self.fields_from_line(f.readline())

            if fields[0].strip() == "METHOD":
                self.method = fields[1].strip()
            else:
                raise ValueError("Not found: A line of METHOD record.")

            # fixed conformer, skip
            fields = list()  # reset
            while len(fields) != 2:
                _, fields = self.fields_from_line(f.readline())

            self.fixedconfs = [int(x) for x in fields[1].strip(" \n").split()]
                
            # free residues
            fields = list()
            while len(fields) != 2:
                line, fields = self.fields_from_line(f.readline())

            n_res = int(fields[0])
            self.free_residues = [[int(xx) for xx in x.strip().split()] 
                                  for x in fields[1].strip(" ;\n").split(";")]
            self.free_residue_names = [self.conformers[x[0]].resid for x in self.free_residues]

            if len(self.free_residues) != n_res:
                msg = "The number of free residues don't match."
                msg = msg + f"\tProblem line:\n\t{line}"
                raise ValueError(msg)

            for ires, res in enumerate(self.free_residues):
                for iconf in res:
                    self.ires_by_iconf[iconf] = ires

            # read MC microstates
            newmc = False
            found_mc = False
            self.microstates_by_id.clear()

            while True:
                line = f.readline()
                if not line:
                    break

                if line.find("MC:") == 0:   # found a MC start
                    newmc = True
                    found_mc = True
                
                    fields = line.split(":")
                    which_mc = int(fields[1].strip())
                    if which_mc in MC_Segments:
                        print(f"Reading: {line.strip()}")
                    continue
                
                elif newmc:
                    #f1, f2 = line.split(":")
                    f2 = line.split(":")[1]
                    current_state = [int(c) for c in f2.split()]
                    newmc = False
                    continue

                elif found_mc:
                    if which_mc in MC_Segments:
                        fields = line.split(",")
                        if len(fields) >= 3:
                            state_e = float(fields[0])
                            count = int(fields[1])
                            flipped = [int(c) for c in fields[2].split()]

                            for ic in flipped:
                                current_state[self.ires_by_iconf[ic]] = ic

                            #print(flipped, current_state)
                            ms = Microstate(current_state, state_e, count)
                            if ms._stateid in self.microstates_by_id:
                                self.microstates_by_id[ms._stateid].count += ms.count
                            else:
                                self.microstates_by_id[ms._stateid] = ms
                            self.counts += ms.count

        # convert microstates to a list
        #self.microstates = [item[1] for item in self.microstates_by_id.items()]
        self.microstates = list(self.microstates_by_id.values())


    def get_occ(self, microstates: list) -> list:
        conf_occ = np.zeros(len(self.conformers))
        total_counts = 0

        for ms in microstates:
            total_counts += ms.count
            for iconf in ms.state():
                conf_occ[iconf] += ms.count

        return (conf_occ/total_counts).tolist()


    def confnames_by_iconfs(self, iconfs):
        # Is iconfs arg a str?
        confnames = [self.conformers[ic].confid for ic in list(iconfs)]
        return confnames
    

    def select_by_conformer(self, microstates: list, conformer_in: list=None) -> tuple:
        """Select microstate if confomer is in the list.
          Return all if the list is empty.
        Args:
          
        Return:
          
        """
        if conformer_in is None:
            return [], microstates

        iconf_in = set([self.iconf_by_confname[confid] for confid in conformer_in])
     
        selected = []
        unselected = []
        for ms in microstates:
            state = set(ms.state())
            if state & iconf_in:
                selected.append(ms)
            else:
                unselected.append(ms)

        return selected, unselected


    def select_by_energy(self, microstates, energy_in: list=None):
        """Select microstate if energy is in the list AND energy is in the range. 
        Return all if the list is empty.
        Args:
          
        Return:
          
        """

        if energy_in is None:
            return [], microstates
        
        selected = []
        unselected = []
        energy_in.sort()

        for ms in microstates:
            if energy_in[0] <= ms.E < energy_in[1]:
                selected.append(ms)
            else:
                unselected.append(ms)

        return selected, unselected


    def convert_to_charge_ms(self) -> list:
        """
        Return a list of charged microstates.
        """

        charge_microstates = []
        charge_ms_by_id = dict()

        for ms in self.microstates:
            current_crg_state = [round(self.conformers[ic].crg) for ic in ms.state()]
            crg_ms = Charge_Microstate(current_crg_state, ms.E*ms.count, ms.count)
            if crg_ms._crg_stateid in charge_ms_by_id:
                charge_ms_by_id[crg_ms._crg_stateid].count += crg_ms.count
                charge_ms_by_id[crg_ms._crg_stateid].total_E += crg_ms.total_E
            else:
                charge_ms_by_id[crg_ms._crg_stateid] = crg_ms
        
        for stateid in charge_ms_by_id:
            crg_ms = charge_ms_by_id[stateid]
            crg_ms.E = crg_ms.average_E = crg_ms.total_E / crg_ms.count
            charge_microstates.append(crg_ms)

        del charge_ms_by_id

        return charge_microstates


    def convert_to_subset_ms(self, res_of_interest: list):
        """
        
        """
        
        iconfs_of_interest = []
        for res in res_of_interest:
            if res in self.free_residue_names:
                ires = self.free_residue_names.index(res)
                conf_select = self.free_residues[ires]
            else:         # this reside is fixed on one or more conformers
                i_fixed = self.fixed_resnames.index(res)
                conf_select = list(self.fixedconfs[i_fixed])

            iconfs_of_interest.append(conf_select)   # This is a list of list
    
        # prepare a list of free residues for grouping microstates
        i_free_res_of_interest = []
        for iconfs in iconfs_of_interest:
            if len(iconfs) > 1:
                i_free_res_of_interest.append(self.ires_by_iconf[iconfs[0]])

        subset_ms_by_id = dict()
        for ms in self.microstates:
            current_sub_state = [ms.state()[i] for i in i_free_res_of_interest]
            #sub_stateid = zlib.compress(" ".join([str(x) for x in current_sub_state]).encode())
            sub_ms = Subset_Microstate(current_sub_state, ms.E*ms.count, ms.count)
                                        
            if sub_ms._subset_stateid in subset_ms_by_id:
                subset_ms_by_id[sub_ms._sub_stateid].count += sub_ms.count
                subset_ms_by_id[sub_ms._sub_stateid].total_E += sub_ms.total_E
            else:
                subset_ms_by_id[sub_ms._sub_stateid] = sub_ms  # new k,v pair to store this ministate

        subset_microstates = []
        for stateid in subset_ms_by_id:
            sub_ms = subset_ms_by_id[stateid]
            sub_ms.E = sub_ms.average_E = sub_ms.total_E / sub_ms.count
            subset_microstates.append(sub_ms)

        del subset_ms_by_id
        
        return subset_microstates
    
    def __repr__(self):
         return "{}({})".format(self.__name__, **vars(self))

    
def get_erange(microstates: list) -> tuple:
    """Return the energy range of the microstates.

    """
    emin = emax = microstates[0].E
    for ms in microstates[1:]:
        if emin > ms.E:
            emin = ms.E
        if emax < ms.E:
            emax = ms.E
    return emin, emax


def bin_mscounts_total(microstates: list, nbins: int=100, erange: list=None) -> tuple:
    """Return two lists, one as energy range and one as the total counts 
    of the in-group microstates.

    """
    if erange:   # use the range if range arrary is provided
        erange.sort()
        n_rng = len(erange)
        counts = np.zeros(n_rng)

        for ms in microstates:
            ibin = -1
            for ie in range(n_rng-1, -1, -1):
                if ms.E > erange[ie]:
                    ibin = ie
                    break
            if ibin >= 0:
                counts[ibin] += ms.count
    else:
        lowest_E, highest_E = get_erange(microstates)
        E_range = highest_E - lowest_E + 0.01

        bin_size = E_range / nbins
        counts = np.zeros(nbins)
        
        for ms in microstates:
            ibin = int((ms.E - lowest_E) / bin_size)
            counts[ibin] += ms.count
        erange = [lowest_E + i*bin_size for i in range(nbins)]

    return erange, counts.tolist()


def bin_mscounts_unique(microstates: list, nbins: int=100, erange: list=None) -> tuple:
    """Return two lists, one as energy range and one as counts of to counts.
    """

    if erange:   # use the range if range arrary is provided
        erange.sort()
        n_rng = len(erange)
        counts = np.zeros(n_rng)

        for ms in microstates:
            ibin = -1
            for ie in range(n_rng-1, -1, -1):
                if ms.E > erange[ie]:
                    ibin = ie
                    break
                
            if ibin >= 0:
                counts[ibin] += 1
    else:
        lowest_E, highest_E = get_erange(microstates)
        E_range = highest_E - lowest_E + 0.01

        bin_size = E_range / nbins
        counts =  np.zeros(nbins)

        # sort by ms energy
        microstates.sort(key=lambda x: x.E)
        for ms in microstates:
            ibin = int((ms.E - lowest_E) / bin_size)
            counts[ibin] += 1

        erange = [lowest_E + i*bin_size for i in range(nbins)]

    return erange, counts.tolist()


def get_count(microstates: list) -> int:
    """Return the sum total of the microstates.count."""
    count = 0
    for ms in microstates:
        count += ms.count
    return count


def average_e(microstates: list) -> float:
    """Return the average energy of the microstates."""
    t = 0.0
    c = 0
    for ms in microstates:
        t += ms.E * ms.count
        c += ms.count
    return t/c


def bhata_distance(prob1: list, prob2: list) -> float:
    """
    Return the Bhattacharyya distance[1] of the distributions prob1 and prob2.
    Args:
      prob1, prob2 (lists):
    Return:
      The Bhattacharyya distance (float) if prob1 and prob2 have the same length,
      else 10,000.0.
    
    [1]: https://en.wikipedia.org/wiki/Bhattacharyya_distance
    """
    d_max = 100_000.0   # Max possible value

    if len(prob1) != len(prob2):
        return d_max
    
    p1 = np.array(prob1)
    p2 = np.array(prob2)

    bc = np.sum(np.sqrt(p1 * p2))
    if bc <= np.exp(-d_max):
        d = d_max
    else:
        d = -np.log(bc)

    return d


def plot_prob(p1:list, p2:list, d_in_title:float):
    """Bar plot of two lists holding conformer occupancies.
    Args:
      p1, p2
      d_in_title : Bhattacharyya dist???
    """

    def label_bars(rects):
        """Annotate the bars (rects) with their heights."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate(f"{height:.3f}",
                        xy=(rect.get_x() + rect.get_width()/2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')
    
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 5))

    x = np.arange(len(p1))

    rects1 = ax.bar(x - width/2, p1, width, label="group1")
    label_bars(rects1)

    rects2 = ax.bar(x + width/2, p2, width, label="group2")
    label_bars(rects2)
               
    ax.set_ylabel("occ")
    ax.set_title(f"d={d_in_title:.3f}")

    ax.legend()
    plt.show()


if __name__ == "__main__":

    msfile = "ms_out/pH5eH0ms.txt"
    #tracemalloc.start()
    mc = MC()

    mc.readms(msfile, mc_numbers=[1,2])
    #print("Loaded ms", tracemalloc.get_traced_memory())


    # Example 1: Bin microstates based on energy
    erange, total_counts = bin_mscounts_total(mc.microstates)
    _, uniq_counts = bin_mscounts_unique(mc.microstates)
    for i in range(len(erange)):
        print("%8.3f %6d %6d" % (erange[i], total_counts[i], uniq_counts[i]))
    #print("bin ms", tracemalloc.get_traced_memory())


    # Example 2: When GLU35 is ionized, what residues change conformation?
    glu_charged_confs = ["GLU-1A0035_011", "GLU-1A0035_012", "GLU-1A0035_013", "GLU-1A0035_011"]
    glu_charged_ms, glu_neutral_ms = mc.select_by_conformer(mc.microstates, conformer_in=glu_charged_confs)
    conf_occ_glu_charged = mc.get_occ(glu_charged_ms)
    conf_occ_glu_neutral = mc.get_occ(glu_neutral_ms)
    for res in mc.free_residues:
        resid = mc.conformers[res[0]].resid
        prob1 = [conf_occ_glu_neutral[ic] for ic in res]
        prob2 = [conf_occ_glu_charged[ic] for ic in res]
        d = bhata_distance(prob1, prob2)
        print("%s, d= %.3f" % (resid, d))
        for ic in res:
            print("%s %6.3f %6.3f" % (mc.conformers[ic].confid, conf_occ_glu_neutral[ic], conf_occ_glu_charged[ic]))
        print()
    #print("Grouping ms", tracemalloc.get_traced_memory())

    # Example 3: Which charge microstate is the most dominant?
    charge_microstates = mc.convert_to_charge_ms()
    charge_microstates.sort(key=lambda x: x.count)
    count = 0
    for crg_ms in charge_microstates:
        count += crg_ms.count
        print(crg_ms.state(), crg_ms.count, crg_ms.average_E)
    print("%d charge microstates" % len(charge_microstates))
    print("%d total microstates" % count)

    #print("charge microstates", tracemalloc.get_traced_memory())


    #tracemalloc.stop()

#===========================================================================
