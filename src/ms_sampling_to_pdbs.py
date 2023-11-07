#!/usr/bin/env python

__doc__ = """
MODULE `ms_sampling_to_pdbs.py` is a command-line executable module that writes a
collection of pdb files from sampled MCCE microstates. The microstates sample size defines
the size of the collection.

PRE-REQUISITES:
---------------
 - MCCE Step4 was run with the `--ms` flag to enable the persistence of the microstate folder, `ms_out`;
 - The folder provided at the command-line is a MCCE simulation output folder that contains the following
   folders or files (required):
     MCCE output:        Info extracted:
     ...................................
    o run.prm.record     STEP$, MONTE_T and MONTE_RUNS
    o name.txt           Atom and cofactors renaming rules
    o step2_out.pdb      Coordinates
    o head3.lst          Conformer index and id, e.g. 00017, GLU-1A0007_005
    o ms_out/            Microstates data from "msout files", e.g. ms_out/pH5eH0ms.txt

NOTE:
-----
The MCCE executable is not required for the purposes of this module.

"""
USAGE = """
> ms_sampling_to_pdbs.py [<arguments>]

Minimal number of arguments: --mcce_dir, --pH, --Eh, --sample_size
> ms_sampling_to_pdbs.py --mcce_dir a/path/to/mcce/output --pH 7 --Eh 0 --sample_size 99

"""

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datetime import datetime
from functools import lru_cache
import numpy as np
import operator
from pathlib import Path
import subprocess
import sys
from typing import Union
import zlib


# ... default params ......................................................
ROOMT = 298.15  # default simulation temp
MONTE_RUNS = 6
VALID_MC_METHODS = [
    "MONTERUNS",
]


def get_msout_filename(
    mcce_output_path: str, pH: Union[float, int], Eh: Union[float, int]
) -> Path:
    """Return the ms_out filename from path, pH and Eh values."""
    if not Path(mcce_output_path).exists():
        raise FileNotFoundError(f"Folder not found: {mcce_output_path}")

    msout_dir = Path(mcce_output_path).joinpath("ms_out")
    if not msout_dir.exists():
        raise FileNotFoundError(f"Folder 'ms_out' not found in {mcce_output_path}")
    if not msout_dir.is_dir():
        raise TypeError(f"'ms_out' must be a directory in {mcce_output_path})")

    prec_ph, prec_eh = 0, 0
    if isinstance(pH, float):
        if not pH.is_integer():
            prec_ph = 1
    if isinstance(Eh, float):
        if not Eh.is_integer():
            prec_eh = 1
    ms_file = f"pH{pH:.{prec_ph}f}eH{Eh:.{prec_eh}f}ms.txt"
    fname = msout_dir.joinpath(ms_file)
    if not fname.exists():
        raise FileNotFoundError(f"File {fname} not found in {msout_dir}")

    return fname


def mkdir_from_msout_file(msout_file: Path) -> tuple:
    """Create a directory with the same name as msout_file w/o the extension.
    Returns:
        tuple: path, created (bool).
    """

    msout_file_dir = msout_file.parent.joinpath(msout_file.stem)
    exists = msout_file_dir.exists()
    if not exists:
        Path.mkdir(msout_file_dir)
        print(f"Created `msout_file_dir` = {msout_file_dir}")

    return msout_file_dir, not exists


def clear_folder(dir_path: str):
    """Delete all files in folder."""

    if not check_dir_noerr(dir_path):
        # no folder, nothing to clear
        return
    p = Path(dir_path)
    for f in p.iterdir():
        if not f.is_dir():
            f.unlink()
    return


def list_folder(dir_path: str):
    """List the files in folder non-recursively."""

    if not check_dir_noerr(dir_path):
        return
    p = Path(dir_path)
    for f in p.iterdir():
        if not f.is_dir():
            print("\t", f)
    return


def check_dir_noerr(dir_path) -> bool:
    """Return True if `dir_path exists and is a folder, else False.
    Usage: Use when it is not critical to raise FileNotFoundError.
    """

    p = Path(dir_path)
    if not p.exists():
        print(f"Path not found: {p}")
        return False
    if not p.is_dir():
        print(f"Path is not a folder: {p}")
        return False

    return True


def check_path(filepath: str, raise_err=True) -> None:
    """
    Returns:
        None if 'all clear', else error if either the parent folder
        or the file itself does not exists.
    """

    p = Path(filepath)
    if not p.exists():
        if raise_err:
            raise FileNotFoundError(f"Path not found: {p}")
        else:
            return False
    if raise_err:
        return
    else:
        return True


def check_msout_split(msout_file_dir: Path, runprm_mcruns: int) -> bool:
    """Return True if the header and the MCi files exist."""
    b = msout_file_dir.joinpath("header").exists()
    if not b:
        # print(f"Split file (header) not found.")
        return False
    for i in range(runprm_mcruns):
        b = b and msout_file_dir.joinpath(f"MC{i}").exists()
        if not b:
            print(f"Split file MC{i} not found.")
            return False
    return True


def read_conformers(head_3_path: str) -> tuple:
    """Return a list of conformers found in head3.lst file and
    a dictionnary, `iconf_by_confname` :: confid -> index.
    Uses Conformer class.
    Used by MS class.
    """

    conformers = []
    iconf_by_confname = {}
    with open(head_3_path) as h3:
        for nl, line in enumerate(h3):
            if nl == 0:
                continue
            if len(line) <= 80:
                continue

            conf = Conformer()
            conf.load_from_head3lst(line)
            conformers.append(conf)
            iconf_by_confname[conf.confid] = nl

    return conformers, iconf_by_confname


def split_msout_file(msout_fname: str, MC_RUNS: int, overwrite: bool = False):
    """Split the msout file into a "header" portion (preceeding MC:0 line) and MCi files
    for the MC records in a folder created with the name of the msout_file as per arguments.
    Note: Each file created is free of comments or blank lines.
    """

    msout_file_dir, created = mkdir_from_msout_file(msout_fname)

    if check_msout_split(msout_file_dir, MC_RUNS) and not overwrite:
        print(
            "The ms_out file is already split into header and MCi files. Set `overwrite` to True to replace them."
        )
        return

    steps_done = {"exper": False, "method": False, "fixed": False, "free": False}
    MC_done = dict((i, False) for i in range(MC_RUNS))

    header_file = msout_file_dir.joinpath("header")
    header = open(header_file, "w")

    # open files for each MC records:: MCi
    for k in MC_done.keys():
        p = msout_file_dir.joinpath(f"MC{k}")
        globals()[f"MC{k}"] = open(p, "w")

    with open(msout_fname) as fh:
        for nl, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            if line[0] == "#":
                continue

            if not steps_done["exper"]:
                if line[0] != "T":
                    raise ValueError(
                        f"The first data line (experiemental variables) must start with T.\n{line}"
                    )
                fields = line.split(",")
                for field in fields:
                    parts = field.split(":")
                    try:
                        value = float(parts[1])
                    except ValueError:
                        print(
                            f"Unrecognized experimental value \
                              (number expected), found: '{parts[1]}')."
                        )
                        raise

                header.write(line + "\n")
                steps_done["exper"] = True
                continue

            if not steps_done["method"]:
                key, value = line.split(":")
                if key.strip() != "METHOD":
                    raise ValueError(
                        f"""This file: {fname} is not a valid Monte Carlo microstate file: the method line
                        starts with {key} instead of METHOD."""
                    )
                if value.strip() not in VALID_MC_METHODS:
                    raise ValueError(
                        f"""This file: {fname} provides an unknown method: {value}. Supported methods are:
                        {VALID_MC_METHODS}."""
                    )
                header.write(line + "\n")
                steps_done["method"] = True
                continue

            if not steps_done["fixed"]:
                header.write(line + "\n")
                steps_done["fixed"] = True
                continue

            if not steps_done["free"]:
                header.write(line + "\n")
                steps_done["free"] = True
                continue

            if line.startswith("MC:0"):
                header.close()
                continue

            if not MC_done[0]:
                MC0.write(line + "\n")
                if line.startswith("MC:1"):
                    MC_done[0] = True
                    MC0.close()
                continue

            if not MC_done[1]:
                MC1.write(line + "\n")
                if line.startswith("MC:2"):
                    MC_done[1] = True
                    MC1.close()
                continue

            if not MC_done[2]:
                MC2.write(line + "\n")
                if line.startswith("MC:3"):
                    MC_done[2] = True
                    MC2.close()
                continue

            if not MC_done[3]:
                MC3.write(line + "\n")
                if line.startswith("MC:4"):
                    MC_done[3] = True
                    MC3.close()
                continue

            if not MC_done[4]:
                MC4.write(line + "\n")
                if line.startswith("MC:5"):
                    MC_done[4] = True
                    MC4.close()
                continue

            if not MC_done[5]:
                MC5.write(line + "\n")

    MC5.close()
    MC_done[5] = True

    return


# ... classes .................................................................
class Conformer:
    """Minimal Conformer class for use in microstate analysis.
    Attributes: iconf, confid, resid, crg.
    """

    def __init__(self):
        self.iconf = 0
        self.confid = ""
        self.resid = ""
        self.crg = 0.0

    def load_from_head3lst(self, line):
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3] + self.confid[5:11]
        self.crg = float(fields[4])


class Microstate:
    """Sortable class for microstates."""

    def __init__(self, state, E, count):
        self.stateid = zlib.compress(" ".join([str(x) for x in state]).encode())
        self.E = E
        self.count = count

    def state(self):
        return [int(i) for i in zlib.decompress(self.stateid).decode().split()]

    def _check_operand(self, other):
        """Fails on missing attribute."""

        if not (
            hasattr(other, "stateid")
            and hasattr(other, "E")
            and hasattr(other, "count")
        ):
            return NotImplemented("Comparison with non Microstate object.")
        return

    def __eq__(self, other):
        self._check_operand(other)
        return (self.stateid, self.E, self.count) == (
            other.stateid,
            other.E,
            other.count,
        )

    def __lt__(self, other):
        self._check_operand(other)
        return (self.stateid, self.E, self.count) < (
            other.stateid,
            other.E,
            other.count,
        )


class MS:
    """The Microstates main class.
    Uses split ms_out files, so processes one MC runs at a time.
    For a given instance, loading of MC data is implemented by `get_mc_data` method,
    which is cached.
    Example:
        ms = MS("some/folder", 7,0,) # selected_MC=0 by default"""

    # do something

    def __init__(
        self,
        mcce_output_path: str,
        pH: Union[int, float],
        Eh: Union[int, float],
        selected_MC: int = 0,
        overwrite_split_files: bool = False,
    ):
        """MS.init

        Parameters:
            mcce_output_path (str): A MCCE simulation output folder.
            pH (int or float): A pH point.
            Eh (int or float): A Eh point.
            selected_MC (int, 0): The index of an MC run; one of `range(constants.MONTERUNS)`.
            overwrite_split_files (bool, False): whether to redo the splitting of msout_file.
        """

        self.mcce_out = Path(mcce_output_path)
        self.pH = pH
        self.Eh = Eh
        self._selected_MC = selected_MC
        self.overwrite_split_files = overwrite_split_files

        self.T = None
        self.MC_RUNS = None
        self.MC_NITER = None
        self._get_runprm_data()

        self.conformers = []
        self.iconf_by_confname = {}
        self._get_conformer_data()

        self.fname = get_msout_filename(self.mcce_out, self.pH, self.Eh)
        self.msout_fpath = self.mcce_out.joinpath(self.fname)
        print("msout_fpath:", self.msout_fpath)
        self.msout_file_dir, created = mkdir_from_msout_file(self.fname)
        print("msout_file_dir, created:", created, self.msout_file_dir)
        if created:
            split_msout_file(self.msout_fpath, self.MC_RUNS)
        elif self.overwrite_split_files:
            clear_folder(self.msout_file_dir)
            split_msout_file(self.msout_fpath, self.MC_RUNS, overwrite=True)

        self.method = ""
        self.fixed_iconfs = []
        self.fixed_residue_names = []
        self.fixed_crg = 0.0
        self.fixed_ne = 0.0
        self.fixed_nh = 0.0
        self.free_residues = []  # free residues, referred by conformer indices
        self.free_residue_names = []
        self.ires_by_iconf = {}  # index of free residue by index of conf
        self._get_header_data()

        self.microstates = []  # list of microstates
        self.counts = 0
        self.get_mc_data(self.selected_MC)

        check_msout_split(self.msout_file_dir, self.MC_RUNS)

    def __repr__(self):
        return f"""{type(self).__name__}("{self.mcce_out}", {self.pH}, {self.Eh}, selected_MC={self.selected_MC}, overwrite_split_files={self.overwrite_split_files})"""

    @property
    def selected_MC(self):
        return self._selected_MC

    @selected_MC.setter
    def selected_MC(self, value):
        if self.T is None:
            # init
            what = MONTE_RUNS
        else:
            what = self.MC_RUNS
        if value >= what:
            print(
                f"This value is beyond the range of MONTE_RUNS({RUNS}) used in this simulation: {value}."
            )
        else:
            self._selected_MC = value

    def __eq__(self, other):
        return self.selected_MC == other.selected_MC

    def __hash__(self):
        return hash(self.selected_MC)

    def _get_runprm_data(self):
        """Populate class vars: T, MC_RUNS, MC_NITER."""
        print("Getting runprm data.")
        runprm = self.mcce_out.joinpath("run.prm.record")
        check_path(runprm)

        prmdata = (
            subprocess.check_output(
                f"grep -E 'MONTE_T)|MONTE_RUNS|MONTE_NITER' {runprm} | sed -e 's/(//g; s/MONTE_//g; s/)//g'",
                stderr=subprocess.STDOUT,
                shell=True,
            )
            .decode()
            .splitlines()
        )

        self.T = float(prmdata[0].split()[0])
        self.MC_RUNS = int(prmdata[1].split()[0])
        self.MC_NITER = int(prmdata[2].split()[0])

        return

    def _get_conformer_data(self):
        """Populate class vars: conformers, iconf_by_confname."""
        print("Getting conformers data.")
        head3_path = self.mcce_out.joinpath("head3.lst")
        check_path(head3_path)
        self.conformers, self.iconf_by_confname = read_conformers(head3_path)

        return

    def _get_header_data(self):
        """Populate class vars: T, pH, Eh, method, fixed_iconfs, free_residues,
        free_residue_names, and ires_by_iconf from the header of the split msout file.
        """
        print("Getting header data.")
        steps_done = {"exper": False, "method": False, "fixed": False, "free": False}

        header_file = self.msout_file_dir.joinpath("header")
        with open(header_file) as fh:
            for nl, line in enumerate(fh):
                if not steps_done["exper"]:
                    fields = line.split(",")
                    for field in fields:
                        parts = field.split(":")
                        key = parts[0].upper().strip()
                        value = float(parts[1])
                        if key == "T":
                            self.T = value
                        elif key == "PH":
                            self.pH = value
                        elif key == "EH":
                            self.Eh = value
                        else:
                            raise ValueError(
                                f"Unrecognized experimental condition part: {key}"
                            )
                    steps_done["exper"] = True
                    continue

                if not steps_done["method"]:
                    _, value = line.split(":")
                    self.method = value.strip()
                    steps_done["method"] = True
                    continue

                if not steps_done["fixed"]:
                    _, fields = line.split(":")
                    self.fixed_iconfs = [int(x) for x in fields.strip("\n").split()]
                    self.fixed_residue_names = [
                        self.conformers[fc].resid for fc in self.fixed_iconfs
                    ]
                    steps_done["fixed"] = True
                    continue

                if not steps_done["free"]:
                    n_res, fres = line.split(":")
                    self.free_residues = [
                        [int(n) for n in grp.strip().split()]
                        for grp in fres.strip(" ;\n").split(";")
                    ]
                    if len(self.free_residues) != int(n_res):
                        msg = "Mismatch between the number of free residues indicator"
                        msg = (
                            msg
                            + " and the number of residues listed on the same line.\n"
                        )
                        raise ValueError(msg + f"\t{line}")

                    self.free_residue_names = [
                        self.conformers[g[0]].resid for g in self.free_residues
                    ]
                    for ires, res in enumerate(self.free_residues):
                        for iconf in res:
                            self.ires_by_iconf[iconf] = ires
                    steps_done["fee"] = True

        return

    @lru_cache(maxsize=None)
    def get_mc_data(self, MC: int):
        """Populate class vars microstates with the data from a MC file identified
        in `self.selected_MC`.
        """
        print("Getting MC data for:", MC)
        if MC != self.selected_MC:
            self.selected_MC = MC

        microstates_by_id = {}

        MC_file = self.msout_file_dir.joinpath(f"MC{self.selected_MC}")
        with open(MC_file) as fh:
            for nl, line in enumerate(fh):
                line = line.strip()

                if nl == 0:
                    _, confs = line.split(":")
                    current_state = [int(c) for c in confs.split()]
                    if not current_state:
                        msg = "The current ms state line cannot be empty.\n"
                        msg = msg + f"\tProblem line in {MC_file}: {nl}"
                        raise ValueError(msg)
                    continue

                fields = line.split(",")
                if len(fields) >= 3:
                    state_e = float(fields[0])
                    count = int(fields[1])
                    # flipped confs:
                    for ic in [int(c) for c in fields[2].split()]:
                        current_state[self.ires_by_iconf[ic]] = ic

                    ms = Microstate(current_state, state_e, count)
                    if ms.stateid in microstates_by_id.keys():
                        microstates_by_id[ms.stateid].count += count
                    else:
                        microstates_by_id[ms.stateid] = ms
                    self.counts += count

        self.microstates = list(microstates_by_id.values())

        return

    def sort_microstates(self, by: str = "E", reverse: bool = False):
        """Return a sorted copy of MS.microstates."""
        by = by.lower()
        if by not in ["energy", "count"]:
            raise ValueError(f"Values for `by` are 'energy' or 'count'; Given: {by}")
        if by[0] == "e":
            by = "E"
        return sorted(self.microstates, key=operator.attrgetter(by), reverse=reverse)

    def sample_microstates(
        self,
        size: int,
        kind: str = "deterministic",
        sort_by: str = "energy",
        reverse: bool = False,
    ) -> tuple:
        """
        Implement a sampling of MS.microstates depending on `kind`.
        Args:
            size (int): sample size
            kind (str, 'deterministic'): Sampling kind: one of ['deterministic', 'random'].
                If 'deterministic', the microstates in ms_list are sorted then sampled at
                regular intervals otherwise, the sampling is random. Case insensitive.
            sort_by (str, "energy"): Only applies if kind is "deterministic".
            reverse (bool, False): Only applies if kind is "deterministic".
        Returns:
            A 3-tuple: cumsum of MC.count in ms.microstates,
                       array of indices for selection,
                       ms_list?
        """

        kind = kind.lower()
        if kind not in ["deterministic", "random"]:
            raise ValueError(
                f"Values for `kind` are 'deterministic' or 'random'; Given: {kind}"
            )

        if kind == "deterministic":
            if (sort_by is None) or (sort_by == "energy"):
                sort_by = "E"
            else:
                sort_by = sort_by.lower()
                if sort_by not in ["energy", "count"]:
                    raise ValueError(
                        f"Values for `sort_by` are 'energy' or 'count'; Given: {by}"
                    )

            ms_list = self.sort_microstates(by=sort_by, reverse=reverse)
            count_selection = np.arange(size, self.counts - size, self.counts / size)
        else:
            ms_list = self.microstates
            rng = np.random.default_rng()
            count_selection = rng.integers(low=1, high=n_counts + 1, size=size)

        ms_cumsum = np.cumsum([mc.count for mc in ms_list])

        return ms_cumsum, count_selection, ms_list

    def get_selected_confs(self, selected_ms):
        """Return the list of conformer ids for selected_ms.
        Args:
            selected_ms (list): List of Microstates objects.
        """

        return [
            conf.confid
            for conf in self.conformers
            if conf.iconf in selected_ms.state() or conf.iconf in self.fixed_iconfs
        ]

    def verify_counts(self):
        """Verify that the microstates counts iteratively added in MS.counts
        equal NITER * len(ms.conformers).
        """
        ms_calc = self.MC_NITER * (len(self.ires_by_iconf) + 1)
        if ms_calc != self.counts:
            print(
                f"Unexpected microstates count: {self.counts = :,}, should be: {ms_calc:,}"
            )

        return

    def get_pdb_remark(ms_index: int):
        """Return a REMARK 250 string with data from the MS.microstates[ms_index]] to prepend in a pdb.

        > REMARK 250 is mandatory if other than X-ray, NMR, neutron, or electron study.
        [Ref]: https://www.wwpdb.org/documentation/file-format-content/format33/remarks1.html
        """

        R250 = """
    REMARK 250
    REMARK 250 EXPERIMENTAL DETAILS
    REMARK 250   EXPERIMENT TYPE               : MCCE simulation
    REMARK 250   DATE OF DATA COLLECTION       : {DATE}
    REMARK 250   REMARK: Date of data collection is the date this pdb was created.
    REMARK 250 EXPERIMENTAL CONDITIONS
    REMARK 250   TEMPERATURE                   : {T:.2f} (K)
    REMARK 250   PH                            : {PH:.2f}
    REMARK 250   EH                            : {EH:.2f}
    REMARK 250   METHOD                        : {METHOD}
    REMARK 250   SELECTED MONTERUN             : {MC}
    REMARK 250   SELECTED MICROSTATE INDEX     : {MS:,}
    REMARK 250   SELECTED MICROSTATE ENERGY    : {E:.2f} (kcal/mol)
    REMARK 250
    """
        dte = datetime.today()
        remark = R250.format(
            DATE=dte.strftime("%d-%b-%y"),
            T=self.T,
            PH=self.pH,
            EH=self.Eh,
            METHOD=self.method,
            MC=self.selected_MC,
            MS=ms_index,
            E=self.microstates[ms_index].E,
        )

        return remark


def ms_to_pdb(
    selected_confs: list,
    ms_index: int,
    mc_run: int,
    remark_data: str,
    step2_path: str,
    output_folder: str,
    output_pdb_format: str,
) -> None:
    """Create a new pdb file in `output_folder` from the `selected_confs`
    Args:
        selected_confs (list): List of selected Microstates objects.
        ms_index (int): Index of selected ms, part of output pdb filename.
        mc_run (int): Index of MC record used, part of output pdb filename.
        remark_data (str): Used to create pdb REMARK 250 section to provide
                           experimental data (T, PH, EH, METHOD) that can be
                           prepended into a pdb.
        step2_path (str): path to step2_out.pdb.
        output_folder (str): path to folder for pdb created from selected_ms.
        output_pdb_format (str, ): PDB format to use when writing the output pdb;
            one of ["standard", "gromacs", "amber"].

    Returns:
        None: The created file names format is f"mc{mc_run}_ms{ms_index}.pdb".

    Pre-requisites:
        The user must acertain that the `selected_confs` come from the same MCCE output files directory
        as the one provided in `step2_path`.
    """

    sep = "_"
    file_name = Path(output_folder).joinpath(f"mc{mc_run}_ms{ms_index}.pdb")
    if file_name.exists():
        print(
            f"File already exists: {file_name}.",
            "Set `clear_pdbs_folder` to True and rerun if you meant to replace it.",
        )
        return

    # step2_out.pdb format:
    # ATOM      1  CA  NTR A0001_001   2.696   5.785  12.711   2.000       0.001      01O000M000 "
    # 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    #         10        20        30        40        50        60        70        80        90
    # len = 91
    with open(step2_path) as pdb:
        with open(file_name, "w") as output_pdb:
            output_pdb.write(remark_data)
            for line in pdb:
                confID = line[17:20] + line[80:82] + line[21:26] + sep + line[27:30]
                if confID[3:5] == "BK":
                    output_pdb.write(line)
                    continue

                if confID in selected_confs:
                    output_pdb.write(line)
    return


def pdbs_from_ms_samples(
    ms: MS,
    mcce_dir: str,
    n_sample_size: int,
    sample_kind: str = "deterministic",
    ms_sort_by: str = "energy",
    ms_sort_reverse: bool = False,
    output_dir: str = None,
    clear_pdbs_folder: bool = False,
    list_files: bool = False,
) -> None:
    """Create `n_sample_size` MCCE_PDB files in `output_dir/pdbs_from_ms`.

    Args:
        ms (MS): A microstate class instance.
        mcce_dir (str): MCCE simulation output folder.
        n_sample_size (int): How many microstates, hence pdb files to create.
        sample_kind (str, "deterministic"): Kind of sampling; either regularly spaced
            or random.
        ms_sort_by (str): Either 'energy' or 'count' (case insensitive).
            Only applies if `sample_kind` is "deterministic".
        ms_sort_reverse (bool, False): If True, descending order.
            Only applies if `sample_kind` is "deterministic".
        output_dir (str, None): Output folder path. Folder "output_dir/pdbs_from_ms"
            will be created if necessary.
        clear_pdbs_folder (bool, False): Whether to delete existing pdb files.
        list_files (bool, False): Whether to list output folder contents.
    """
    ms_sort_by = ms_sort_by.lower()
    if ms_sort_by not in ["energy", "count"]:
        raise ValueError(
            f"Values for `ms_sort_by` are 'energy' or 'count'; Given: {ms_sort_by}"
        )
    if ms_sort_by == "energy":
        ms_sort_by = "E"

    mcce_dir = Path(mcce_dir)
    step2_path = mcce_dir.joinpath("step2_out.pdb")
    check_path(step2_path)

    if not output_dir or output_dir is None:
        output_dir = ms.msout_file_dir

    pdb_out_folder = Path(output_dir).joinpath("pdbs_from_ms")
    if not pdb_out_folder.exists():
        Path.mkdir(pdb_out_folder)
    elif clear_pdbs_folder:
        clear_folder(pdb_out_folder)

    ms_cumsum, count_selection, ms_sampled = ms.sample_microstates(
        n_sample_size, kind=sample_kind, sort_by=ms_sort_by, reverse=ms_sort_reverse
    )

    # Summarize what's being done:
    print(
        f"Creating n={n_sample_size:,} MCCE_PDB files in {output_dir} from (n) microstates sorted by '{ms_sort_by}'.\n",
        "NOTE: the output pdb will be free of any water molecules in step2_out.pdb.",
    )

    mc_run = ms.selected_MC  # part of pdb name
    for c in count_selection:
        ms_index = np.where((ms_cumsum - c) > 0)[0][0]
        ms_selection = ms_sampled[ms_index]

        confs_for_pdb = ms.get_selected_confs(ms_selection)

        # gather initial data for REMARK section of pdb:
        remark_data = get_pdb_remark(ms, ms_index)
        # write the pdb in the folder
        ms_to_pdb(
            confs_for_pdb, ms_index, mc_run, remark_data, step2_path, pdb_out_folder
        )
        # pdb names: = Path(pdb_out_folder).joinpath(f"mc{mc_run}_ms{ms_index}.pdb")

    print("PDB files creation over.")
    if list_files:
        print(f"Files in {pdb_out_folder}:\n")
        list_folder(pdb_out_folder)

    return


# ... cli ....................................................................
def generate_parser():
    def arg_int_or_float(x):
        """Replaces typing with Union[int, float] does not work in argparse."""
        x = float(x)
        if x.is_integer():
            return int(x)
        else:
            return x

    p = ArgumentParser(
        prog=__name__,
        description=__doc__,
        usage=USAGE,
        formatter_class=RawDescriptionHelpFormatter,
        epilog=">>> END of %(prog)s.",
    )

    p.add_argument(
        "--mcce_dir",
        type=str,
        required=True,
        help="The folder with files from a MCCE simulation; required.",
    )
    p.add_argument(
        "--pH",
        type=arg_int_or_float,
        required=True,
        help="The pH point; part of experiemntal variables defining a microstate; required.",
    )
    p.add_argument(
        "--Eh",
        type=arg_int_or_float,
        required=True,
        help="The Eh point; part of experiemntal variables defining a microstate; required.",
    )
    p.add_argument(
        "-MC",
        type=list,
        default=[0],
        help="List for (zero-based) indices of the MONTERUNS to use; default: %(default)s.",
    )
    p.add_argument(
        "-overwrite_split_files",
        default=False,
        action="store_false",
        # help='the bar to %(prog)s (default: %(default)s)'
        help="The MS class uses a split msout file (header and MCi files); if True the file will be split anew; default: %(default)s.",
    )
    p.add_argument(
        "--sample_size",
        type=int,
        required=True,
        help="The size of the microstates sample, hence the number of pdb files to write; required",
    )
    p.add_argument(
        "-sample_kind",
        type=str,
        choices=["d", "deterministic", "r", "random"],
        default="d",
        help="""The sampling kind: 'deterministic': regularly spaced samples,
        'random': random indices over the microstates space chosen; default: %(default)s.""",
    )
    p.add_argument(
        "-sort_by",
        type=str,
        choices=["energy", "count"],
        default="energy",
        help="""The name referencing the Conformer variable to use when sorting:
        'energy'-> Conf.E, 'count': Conf.count; default: %(default)s.""",
    )
    p.add_argument(
        "-reverse_sort",
        default=False,
        action="store_false",
        help="The sort order: False : ascendingly, True : descendingly; default: %(default)s.",
    )
    p.add_argument(
        "-output_dir",
        type=str,
        default=None,
        help="""The path of the folder receiving the created pdb files. If not provided, the path
        defaults to mcce_dir/ms_out/msout_file_dir/pdbs_from_ms, otherwise the actual output folder
        will be: output_dir/pdbs_from_ms; default: %(default)s.""",
    )
    p.add_argument(
        "-clear_pdbs_folder",
        default=False,
        action="store_false",
        help="Whether to clear an existing pdbs folder; default: %(default)s.",
    )
    p.add_argument(
        "-list_files",
        default=False,
        action="store_false",
        help="Whether to list the pdb files created; default: %(default)s.",
    )

    return p


def check_mcce_dir(mcce_dir):
    """Check the integrity of the mcce output folder:
    0. mcce_dir/run.prm.record file exists
    1. Step4 was run as indicated in run.prm.record
    2. These also exist: head3.lst, step2_out.pdb, and ms_out folder
    """

    for f in ["run.prm.record", "name.txt", "head3.lst", "step2_out.pdb", "ms_out"]:
        fp = mcce_dir.joinpath(f)
        if not check_path(fp, raise_err=False):
            raise EnvironmentError("MCCE output missing: {fp}.")

    step4_done = None
    runprm = mcce_dir.joinpath("run.prm.record")
    try:
        step4_done = subprocess.check_output(
            f"grep 'STEP4' {runprm}", stderr=subprocess.STDOUT, shell=True
        ).decode()
    except CalledProcessError:
        pass

    if step4_done is None:
        raise EnvironmentError(
            "MCCE Step4 not done (missing in run.prm.record), but required."
        )

    return


def cli_ms2pdb(argv=None):
    """"""

    cli_parser = generate_parser()
    args = cli_parser.parse_args(argv)

    if args is None:
        cli_parser.print_help()
        return

    mcce_dir = Path(args.mcce_dir)
    check_mcce_dir(mcce_dir)

    kind = "deterministic"
    if args.sample_kind[0] == "r":
        kind = "random"

    mc = args.MC.pop(0)
    # use 1st MC to populate all attributes:
    ms = MS(
        mcce_dir,
        args.pH,
        args.Eh,
        selected_MC=mc,
        overwrite_split_files=args.overwrite_split_files,
    )
    # write 1st batch of sampled pdbs:
    pdbs_from_ms_samples(
        ms,
        mcce_dir,
        args.sample_size,
        sample_kind=kind,
        ms_sort_by=args.sort_by,
        ms_sort_reverse=args.reverse_sort,
        output_dir=args.output_dir,
        clear_pdbs_folder=args.clear_pdbs_folder,
        list_files=args.list_files,
    )

    while args.MC:
        mc = args.MC.pop(0)
        # update the class attribute for MC selection:
        ms.selected_MC = mc

        # write another batch of sampled pdbs:
        pdbs_from_ms_samples(
            ms,
            mcce_dir,
            args.sample_size,
            sample_kind=kind,
            ms_sort_by=args.sort_by,
            ms_sort_reverse=args.reverse_sort,
            output_dir=args.output_dir,
            clear_pdbs_folder=args.clear_pdbs_folder,
            list_files=args.list_files,
        )
    return


if __name__ == "__main__":
    sys.exit(cli_ms2pdb(sys.argv[1:]))