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

import zlib
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt


# the constants should be in their own module
ph2Kcal = 1.364
Kcal2kT = 1.688


@dataclass
class Conformer:
    """Class of a conformer."""

    iconf: int = 0
    confid: str = ""
    resid: str = ""
    crg: float = 0.0

    def load(self, line):
        """Populate self.[iconf,confid,resid,crg] from head3 line."""
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3] + self.confid[5:11]
        self.crg = float(fields[4])


class Microstate:
    """Base class for building a microstate object."""

    def __init__(self, state: list, E: float, count: int):
        self.stateid = zlib.compress(" ".join([str(x) for x in state]).encode())
        self.E = E
        self.count = count
        self.average_E = 0.0
        self.total_E = 0.0

    def state(self):
        """Decode the ._stateid binary variable."""
        return [int(i) for i in zlib.decompress(self.stateid).decode().split()]


# Not needed:
class Charge_Microstate(Microstate):
    """Class for a charge subset."""

    def __init__(self, crg_stateid, total_E, count):
        Microstate.__init__(self, crg_stateid, total_E, count)

        self.crg_stateid = self.stateid  # rename
        self.E = 0.0
        self.total_E = total_E


class Subset_Microstate(Microstate):
    """Class for a subset."""

    def __init__(self, subset_stateid, total_E, count):
        Microstate.__init__(self, subset_stateid, total_E, count)

        self.subset_stateid = self.stateid  # rename
        self.E = 0.0
        self.total_E = total_E


class MC:
    """
    Class for a microstate.
    Methods:
    """

    def __init__(self, T: float = 298.15, pH: float = 7.0, Eh: float = 0.0):
        """T, pH, Eh will determine which ms_out file to read."""
        self.T = T
        self.pH = pH
        self.Eh = Eh

        self.N_ms = 0
        self.N_uniq = 0

        self.lowest_E = 0.0
        self.highest_E = 0.0
        self.average_E = 0.0

        self.fixed_crg = 0.0
        self.fixed_ne = 0.0
        self.fixed_nh = 0.0
        self.fixed_confs = []  # fixed conformers
        self.fixed_residue_names = []

        self.free_residues = (
            []
        )  # list of conformer groups constituents of free residues
        self.free_residue_names = []  # free residue names
        self.iconf2ires = {}  # dict, from conformer index to free residue index
        self.ires_by_iconf = {}  # dict, index of free residue by index of conf
        self.iconf_by_confname = {}  # dict

        self.microstates = []  # a list of microstates
        self.microstates_by_id = {}  # dict
        self.conformers = []

        self.method = ""
        self.counts = 0

        # populate self.conformers and self.iconf_by_confname:
        with open("head3.lst") as f:
            lines = f.readlines()

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
        Used by MC.readms.
        """
        line = line.split("#")[0]  # remove comment
        fields = line.split(",")

        return line, fields

    def readms(self, fname: str, mc_numbers: list = None):
        """MC class method. Read a microstate (ms) file in ms_out folder (a.k.a "a ticker
        tape" file), for example, "/ms_out/pH7eH0ms.txt". All runs of Monte Carlo sampling
        (currently 6) will be included with default `mc_numbers` (None), else the one or more
        runs of interest must be passed as a list.  These numbers correspond to the `MC:n` line
        in the microstate file.

        Args:
          fname (str): the name of microstate file as found in mc_out folder.
          mc_numbers (list(int)): list of MC runs to consider.
        """

        if mc_numbers:
            MC_Segments = mc_numbers
        else:
            msg = (
                "You have not selected any specific MC run (`mc_number` arg is None),\n"
            )
            msg = msg + "so all (6) will be processed.\n"
            msg = (
                msg
                + "The call to process selected MC runs is, e.g. `MC.readms(msfile, "
            )
            msg = msg + "mc_numbers=[2,3])\n"
            msg = msg + "Enter Y/y if you want all runs or anything else to redo: "
            # Ask for redo:
            reply = input(msg).lower()[0]
            if reply != "y":
                return

            MC_Segments = [0, 1, 2, 3, 4, 5]

        # Only start variables assignment after a possible exit point (if possible).
        self.microstates = []  # reset
        self.counts = 0  # reset

        fields = []
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
                parts = field.split(":")
                key = parts[0].upper().strip()

                try:
                    value = float(parts[1])
                except ValueError:
                    print(
                        f"Unrecognized experimental value \
                          (number expected), found: '{parts[1]}')."
                    )

                if key == "T":
                    self.T = value
                elif key == "PH":
                    self.pH = value
                elif key == "EH":
                    self.Eh = value
                else:
                    raise ValueError(f"Unrecognized experimental condition part: {key}")

            # method record
            fields = []  # reset
            while len(fields) != 2:
                _, fields = self.fields_from_line(f.readline())

            if fields[0].strip() == "METHOD":
                self.method = fields[1].strip()
            else:
                raise ValueError("Not found: A line for METHOD record.")

            # fixed conformer, skip
            fields = []  # reset
            while len(fields) != 2:
                _, fields = self.fields_from_line(f.readline())

            self.fixed_confs = [int(x) for x in fields[1].strip(" \n").split()]

            # free residues
            fields = []
            while len(fields) != 2:
                line, fields = self.fields_from_line(f.readline())

            n_res = int(fields[0])
            self.free_residues = [
                [int(xx) for xx in x.strip().split()]
                for x in fields[1].strip(" ;\n").split(";")
            ]
            self.free_residue_names = [
                self.conformers[x[0]].resid for x in self.free_residues
            ]

            if len(self.free_residues) != n_res:
                msg = "Mismatch between the number of free residues \
                    indicator and the number of residues listed on the same line."
                msg = msg + f"\tProblem line:\n\t{line}"
                raise ValueError(msg)

            for ires, res in enumerate(self.free_residues):
                for iconf in res:
                    self.ires_by_iconf[iconf] = ires

            # read MC microstates
            new_mc = False
            found_mc = False
            mc_in_segments = False

            self.microstates_by_id.clear()

            while True:
                line = f.readline()
                if not line:
                    break

                if line.find("MC:") == 0:  # found a MC start
                    new_mc = True
                    found_mc = True

                    fields = line.split(":")
                    mc_in_segments = int(fields[1]) in MC_Segments
                    if mc_in_segments:
                        print(f"Reading: {line}")

                    continue

                if new_mc:
                    current_state = [int(c) for c in line.split(":")[1].split()]
                    if not current_state:
                        msg = "The current ms state line cannot be empty.\n"
                        msg = msg + f"\tProblem line:\n\t{line}"
                        raise ValueError(msg)

                    new_mc = False
                    continue

                if found_mc:
                    if mc_in_segments:
                        fields = line.split(",")
                        if len(fields) >= 3:
                            state_e = float(fields[0])
                            count = int(fields[1])
                            flipped = [int(c) for c in fields[2].split()]

                            for ic in flipped:
                                current_state[self.ires_by_iconf[ic]] = ic

                            # print(flipped, current_state)
                            ms = Microstate(current_state, state_e, count)

                            if ms.stateid in self.microstates_by_id:
                                self.microstates_by_id[ms.stateid].count += ms.count
                            else:
                                self.microstates_by_id[ms.stateid] = ms
                            self.counts += ms.count

        # convert microstates to a list
        self.microstates = list(self.microstates_by_id.values())

    def get_occ(self, microstates: list) -> list:
        """Return the average occupancy in each microstate."""

        conf_occ = np.zeros(len(self.conformers))
        total_counts = 0

        for ms in microstates:
            total_counts += ms.count
            for iconf in ms.state():
                conf_occ[iconf] += ms.count

        return (conf_occ / total_counts).tolist()

    def confnames_by_iconfs(self, iconfs):
        """Return the conformers id given their numbers."""

        confnames = [self.conformers[ic].confid for ic in list(iconfs)]
        return confnames

    def select_by_conformer(
        self, microstates: list, conformer_in: list = None
    ) -> tuple:
        """Select microstate if confomer is in the list.
        Return all if the list is empty.
        Args:
        Return:
        """

        if conformer_in is None:
            return [], microstates

        selected = []
        unselected = []

        # set comprehension (no extra list created):
        iconf_in = {self.iconf_by_confname[confid] for confid in conformer_in}

        for ms in microstates:
            state = set(ms.state())
            if state & iconf_in:
                selected.append(ms)
            else:
                unselected.append(ms)

        return selected, unselected

    def select_by_energy(self, microstates, energy_in: list = None):
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
        charge_ms_by_id = {}  # dict

        for ms in self.microstates:
            current_crg_state = [round(self.conformers[ic].crg) for ic in ms.state()]
            crg_ms = Charge_Microstate(current_crg_state, ms.E * ms.count, ms.count)

            if crg_ms.crg_stateid in charge_ms_by_id:
                charge_ms_by_id[crg_ms.crg_stateid].count += crg_ms.count
                charge_ms_by_id[crg_ms.crg_stateid].total_E += crg_ms.total_E
            else:
                charge_ms_by_id[crg_ms.crg_stateid] = crg_ms

        for stateid in charge_ms_by_id.keys():
            crg_ms = charge_ms_by_id[stateid]
            crg_ms.E = crg_ms.average_E = crg_ms.total_E / crg_ms.count
            charge_microstates.append(crg_ms)

        del charge_ms_by_id

        return charge_microstates

    def convert_to_subset_ms(self, res_of_interest: list):
        """Reduce the ms space to that of the residues in `res_of_interest`."""

        iconfs_of_interest = []
        for res in res_of_interest:
            if res in self.free_residue_names:
                ires = self.free_residue_names.index(res)
                conf_select = self.free_residues[ires]
            else:  # this residue is fixed on one or more conformers
                i_fixed = self.fixed_residue_names.index(res)
                conf_select = list(self.fixed_confs[i_fixed])

            iconfs_of_interest.append(conf_select)  # This is a list of list

        # prepare a list of free residues for grouping microstates
        i_free_res_of_interest = []
        for iconfs in iconfs_of_interest:
            if len(iconfs) > 1:
                i_free_res_of_interest.append(self.ires_by_iconf[iconfs[0]])

        subset_ms_by_id = {}  # dict
        for ms in self.microstates:
            current_sub_state = [ms.state()[i] for i in i_free_res_of_interest]
            sub_ms = Subset_Microstate(current_sub_state, ms.E * ms.count, ms.count)

            if sub_ms.subset_stateid in subset_ms_by_id:
                subset_ms_by_id[sub_ms.subset_stateid].count += sub_ms.count
                subset_ms_by_id[sub_ms.subset_stateid].total_E += sub_ms.total_E
            else:
                # new k,v pair to store this ministate
                subset_ms_by_id[sub_ms.subset_stateid] = sub_ms

        subset_microstates = []
        for stateid in subset_ms_by_id.items():
            sub_ms = subset_ms_by_id[stateid]
            sub_ms.E = sub_ms.average_E = sub_ms.total_E / sub_ms.count
            subset_microstates.append(sub_ms)

        del subset_ms_by_id

        return subset_microstates

    def __repr__(self):
        return "{}({})".format(type(self).__name__, **vars(self))


def get_erange(microstates: list) -> tuple:
    """Return the energy range of the microstates."""

    emin = emax = microstates[0].E
    for ms in microstates[1:]:
        if emin > ms.E:
            emin = ms.E
        if emax < ms.E:
            emax = ms.E
    return emin, emax


def bin_mscounts_total(
    microstates: list, nbins: int = 100, erange: list = None
) -> tuple:
    """Return two lists, one as energy range and one as the total counts
    of the in-group microstates.
    """

    if erange:  # use the provided range list
        erange.sort()
        n_rng = len(erange)
        counts = np.zeros(n_rng)

        for ms in microstates:
            ibin = -1
            for ie in range(n_rng - 1, -1, -1):
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
        erange = [lowest_E + i * bin_size for i in range(nbins)]

    return erange, counts.tolist()


def bin_mscounts_unique(
    microstates: list, nbins: int = 100, erange: list = None
) -> tuple:
    """Return two lists, one as energy range and one as counts of unique counts."""

    if erange:  # use the range if range arrary is provided
        erange.sort()
        n_rng = len(erange)
        counts = np.zeros(n_rng)

        for ms in microstates:
            ibin = -1
            for ie in range(n_rng - 1, -1, -1):
                if ms.E > erange[ie]:
                    ibin = ie
                    break

            if ibin >= 0:
                counts[ibin] += 1
    else:
        lowest_E, highest_E = get_erange(microstates)
        E_range = highest_E - lowest_E + 0.01

        bin_size = E_range / nbins
        counts = np.zeros(nbins)

        # sort by ms energy
        microstates.sort(key=lambda x: x.E)
        for ms in microstates:
            ibin = int((ms.E - lowest_E) / bin_size)
            counts[ibin] += 1

        erange = [lowest_E + i * bin_size for i in range(nbins)]

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
    return t / c


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
    d_max = 100_000.0  # Max possible value

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


def plot_occ_pair(p1: list, p2: list, d_in_title: float = None):
    """Bar plot of two lists holding conformer occupancies.
    Args:
      p1, p2
      d_in_title : Bhattacharyya dist???
    """

    def label_bars(rects):
        """Annotate the bars (rects) with their heights."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate(
                f"{height:.3f}",
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha="center",
                va="bottom",
            )

    width = 0.35

    _, ax = plt.subplots(figsize=(8, 5))

    x = np.arange(len(p1))

    rects1 = ax.bar(x - width / 2, p1, width, label="group1")
    label_bars(rects1)

    rects2 = ax.bar(x + width / 2, p2, width, label="group2")
    label_bars(rects2)

    ax.set_ylabel("occupancy")

    if d_in_title is None:
        ax.set_title("d=<Bhattacharyya dist.?>")
    else:
        ax.set_title(f"d={d_in_title:.3f}")

    ax.legend()
    plt.show()


def example_1(mc: MC):
    """Example 1: Bin microstates based on energy."""
    print(__doc__)

    erange, total_counts = bin_mscounts_total(mc.microstates)
    _, uniq_counts = bin_mscounts_unique(mc.microstates)
    for i, rng in enumerate(erange):
        msg = "{:,.3f}\t{:,d}t{:,d}"
        print(msg.format(rng, total_counts[i], uniq_counts[i]))


def example_2(mc: MC):
    """Example 2: When GLU35 is ionized, what residues change conformation?"""
    print(__doc__)

    glu_charged_confs = [
        "GLU-1A0035_011",
        "GLU-1A0035_012",
        "GLU-1A0035_013",
        "GLU-1A0035_011",
    ]

    glu_charged_ms, glu_neutral_ms = mc.select_by_conformer(
        mc.microstates, conformer_in=glu_charged_confs
    )
    conf_occ_glu_charged = mc.get_occ(glu_charged_ms)
    conf_occ_glu_neutral = mc.get_occ(glu_neutral_ms)

    for res in mc.free_residues:
        resid = mc.conformers[res[0]].resid
        prob1 = [conf_occ_glu_neutral[ic] for ic in res]
        prob2 = [conf_occ_glu_charged[ic] for ic in res]
        d = bhata_distance(prob1, prob2)
        print(f"{resid}, d= {d:.3f}")

        for ic in res:
            print(
                f"{mc.conformers[ic].confid}\t",
                f"{conf_occ_glu_neutral[ic]:.3f}\t",
                f"{conf_occ_glu_charged[ic]:.3f}",
            )


def example_3(mc: MC):
    """Example 3: Which charge microstate is the most dominant?"""
    print(__doc__)

    charge_microstates = mc.convert_to_charge_ms()
    charge_microstates.sort(key=lambda x: x.count)

    count = 0
    for crg_ms in charge_microstates:
        count += crg_ms.count
        print(crg_ms.state(), crg_ms.count, crg_ms.average_E)
    print(
        f"Number of charge microstates: {len(charge_microstates)}\n",
        f"Total number of microstates: {count:,}",
    )


def main():
    """The main function ia automatically run when the module is
    run at the command line, else it is available when the module
    is imported.
    """
    import tracemalloc

    ms_file = "ms_out/pH5eH0ms.txt"
    print(f"Running examples using {ms_file}")

    tracemalloc.start()
    mc = MC()

    mc.readms(ms_file, mc_numbers=[1, 2])
    print("MEM", "Loaded ms", tracemalloc.get_traced_memory())

    example_1(mc)
    print("MEM", "Binning ms", tracemalloc.get_traced_memory())

    example_2(mc)
    print("MEM", "Grouping ms", tracemalloc.get_traced_memory())

    example_3(mc)
    print("MEM", "Charge microstates", tracemalloc.get_traced_memory())

    tracemalloc.stop()


if __name__ == "__main__":
    main()
