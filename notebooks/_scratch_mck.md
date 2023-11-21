---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python [conda env:mce]
    language: python
    name: conda-env-mce-py
---

```python jupyter={"source_hidden": true}
import sys
from pathlib import Path
import time
import numpy as np
from pprint import pprint as pp
import matplotlib as mpl
from matplotlib import pyplot as plt
plt.ion()
#plt.style.use('seaborn-v0_8-muted')

from IPython.display import Markdown
# To get multiple outputs into 1 cell w/o using print:
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

# autoreload extension
from IPython import get_ipython

ipython = get_ipython()
if 'autoreload' not in ipython.extension_manager.loaded:
    %load_ext autoreload
%autoreload 2

# -----------------------------------------
# TWO USEFUL FUNCTIONS:

def add_to_sys_path(this_path, up=False):
    """
    Prepend this_path to sys.path.
    If up=True, path refers to parent folder (1 level up).
    """
    if up:
        newp = Path(this_path).parent
    else:
        newp = Path(this_path)
    src = newp.joinpath("src")
    if src.exists():
        newp = str(src)
    else:
        newp = str(newp)
    if newp not in sys.path:
        sys.path.insert(1, newp)
        print('Path added to sys.path: {}'.format(newp))


# Filtered dir() for method discovery:
def fdir(obj, start_with_str='_', exclude=True):
    return [d for d in dir(obj) if not d.startswith(start_with_str) == exclude]

```

```python jupyter={"source_hidden": true}
# Insert current src dir into sys.path so that modules in ../src can be imported:
# CHANGE THIS IF NEEDED:

add_to_sys_path(Path.cwd(), up=True)
```

---
# Test executable module

```python
from functools import partial
import subprocess

import ms_sampling_to_pdbs as sampling

load_npz = partial(np.load, allow_pickle=True)
save_npz = partial(np.savez, allow_pickle=True)
```

```python
HERE = Path.cwd()  # do not change

test_folder = "granepura_GunnerLab/mcce_data" # you can change this to any other mcce folder
mcce_dir = HERE.parent.parent.joinpath(test_folder)  # do not change
mcce_dir, mcce_dir.exists()
```

<!-- #raw jupyter={"source_hidden": true} -->
mcce_dir = DATA.joinpath("4lzt")
mcce_dir
mcce_dir = Path("../data/4lzt").resolve()
mcce_dir

nametxt = mcce_dir.joinpath("name.txt")
runprm = mcce_dir.joinpath("run.prm.record")
runprm2 = runprm.relative_to(Path.cwd().parent)  # use when running subprocesses from this nbk

msout_dir = mcce_dir.joinpath("ms_out")
!ls -l {msout_dir}
<!-- #endraw -->

<!-- #raw -->
Markdown(filename="../README.md")
<!-- #endraw -->

<!-- #raw jupyter={"source_hidden": true} -->
TEST_DIR = Path.cwd().parent.joinpath("tests")
TEST_DIR

test_runprm = TEST_DIR.joinpath("run.prm.record")
<!-- #endraw -->

```python
pH = 7.0
Eh= 0.0
overwrite = False

ms = sampling.MS(mcce_dir,pH, Eh, overwrite_split_files = overwrite)
ms.T, ms.MC_RUNS
```

```python
ms.get_mc_data(0)
print(f"{ms.counts = :,}; {ms.MC_RUNS = :,}")

print(f"{len(ms.microstates) = :,}")
```

```python
#which = msout_dir.joinpath("pH6eH0ms")
!ls -l {ms.msout_file_dir}
```

```python
print(ms)
ms.mcce_out
ms.msout_file_dir
ms.msout_fpath
ms.msout_fpath.name
```

```python

```

# Check output of  sampling

<!-- #raw jupyter={"source_hidden": true} -->
#ok

mc_run = ms.selected_MC  # part of pdb name
mc_run

n_sample_size = 4
sample_kind = "deterministic"
ms_sort_by = "energy"
output_pdb_format = "mcce"  # not implemented, no effect

if ms.sampled_cumsum is None:
    ms.get_sampling_params(n_sample_size,
                           kind=sample_kind,
                           sort_by=ms_sort_by
                           )

if ms.sampled_ms is None:
    ms_list = ms.microstates
else:
    ms_list = ms.sampled_ms

#...................................................
step2_path = ms.mcce_out.joinpath("step2_out.pdb")
pdb_out_folder = Path.cwd()
if pdb_out_folder.name != "pdbs_from_ms":
    pdb_out_folder = pdb_out_folder.joinpath("pdbs_from_ms")

if not pdb_out_folder.exists():
    Path.mkdir(pdb_out_folder)
#...................................................

for c in ms.sampled_count_selection:
    ms_index = np.where((ms.sampled_cumsum - c) > 0)[0][0]
    ms_selection = ms_list[ms_index]

    confs_for_pdb = ms.get_selected_confs(ms_selection)

    # gather initial data for REMARK section of pdb:
    remark_data = ms.get_pdb_remark(ms_index)
    # write the pdb in the folder
    sampling.ms_to_pdb(confs_for_pdb, ms_index, mc_run,
                       remark_data, step2_path,
                       pdb_out_folder, output_pdb_format)

confs_for_pdb[:6]
<!-- #endraw -->

<!-- #raw jupyter={"source_hidden": true} -->
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
                    for i, res in enumerate(self.free_residues):
                        for iconf in res:
                            self.ires_by_iconf[iconf] = i
                    steps_done["fee"] = True
<!-- #endraw -->

<!-- #raw jupyter={"source_hidden": true} -->
# from jj_ms_analysis.py

def convert_to_subset_ms(self, res_of_interest):
    """What does this do?
    """
    iconfs_of_interest = []

    for ires, res in enumerate(res_of_interest):
        if res in self.free_residue_names:
            conf_select = self.free_residues[ires]
        else:  # this reside is fixed on one or more conformers
            i_fixed = fixed_resnames.index(res)
            conf_select = [self.fixedconfs[i_fixed]]

        iconfs_of_interest.append(conf_select)  # a list of list

    # prepare a list of free residues for grouping microstates
    i_free_res_of_interest = [self.ires_by_iconf[iconfs[0]]
                             for iconfs in iconfs_of_interest
                             if len(iconfs) > 1
                             ]
    subset_ms_by_id = {}
    for ms in self.microstates:
        current_sub_state = [ms.state()[i] for i in i_free_res_of_interest]
        sub_stateid = zlib.compress(
            " ".join([str(x) for x in current_sub_state]).encode()
        )
        sub_ms = Subset_Microstate(sub_stateid, ms.E * ms.count, ms.count)

        if sub_stateid in subset_ms_by_id:
            subset_ms_by_id[sub_stateid].count += sub_ms.count
            subset_ms_by_id[sub_stateid].total_E += sub_ms.total_E
        else:
            subset_ms_by_id[
                sub_stateid
            ] = sub_ms  # create a new key value pair to store this ministate

    subset_microstates = []
    for sub_stateid in subset_ms_by_id:
        sub_ms = subset_ms_by_id[sub_stateid]
        sub_ms.E = sub_ms.average_E = sub_ms.total_E / sub_ms.count
        subset_microstates.append(sub_ms)

    return subset_microstates
<!-- #endraw -->

<!-- #raw -->
h3 = mcce_dir.joinpath("head3.lst")
s2 = mcce_dir.joinpath("step2_out.pdb")
!head -n 20 {h3}
<!-- #endraw -->

<!-- #raw -->
!head -n 200 {s2} | tail -n 20
<!-- #endraw -->

<!-- #raw jupyter={"source_hidden": true} -->
    for ires, res in enumerate(res_of_interest):
        if res in self.free_residue_names:
            conf_select = self.free_residues[ires]
        else:  # this reside is fixed on one or more conformers
            i_fixed = fixed_resnames.index(res)
            conf_select = [self.fixedconfs[i_fixed]]

        iconfs_of_interest.append(conf_select)  # a list of list

    # prepare a list of free residues for grouping microstates
    i_free_res_of_interest = [self.ires_by_iconf[iconfs[0]]
                             for iconfs in iconfs_of_interest
                             if len(iconfs) > 1
                             ]
<!-- #endraw -->

<!-- #raw jupyter={"source_hidden": true} -->
mc_run = ms.selected_MC  # part of pdb name
mc_run

n_sample_size = 4
sample_kind = "deterministic"
ms_sort_by = "energy"

if ms.sampled_cumsum is None:
    ms.get_sampling_params(n_sample_size,
                           kind=sample_kind,
                           sort_by=ms_sort_by
                           )

if ms.sampled_ms is None:
    ms_list = ms.microstates
else:
    ms_list = ms.sampled_ms


all_confs = ms.get_selected_confs(ms.microstates[0])
for i, c in enumerate(all_confs):
    if i >=8:
        break
    print(c[0], c[1])

ms.free_residue_names[:4]
ms.free_residues[:4]

#len(ms.ires_by_iconf) #273
for i, kv in enumerate(ms.ires_by_iconf.items()):
    if i >= 8:
        break
    print(kv)
<!-- #endraw -->

<!-- #raw jupyter={"source_hidden": true} -->
def get_ms_from_smsm(smsm_data: np.ndarray, sel_index: int):
    """
    Retrive a column from the smsm matrix as a Microstate object.
    Args:
        smsm_data (np.ndarray): data from npz file.
    """

    sampling_info = smsm_data["info"][0].split(", ")
    if "deterministic" in sampling_info[1]:
        lst = [info.split("=")[1] for info in sampling_info]
        sort_by = lst[2][1:-1]
        if lst[3] == "False":
            reverse = False
        else:
            reverse = True
        ms_list = ms.sort_microstates(by=sort_by, reverse=reverse)
    else:
        ms_list = ms.microstates

    smsm = smsm_data["smsm"]
    col_index = smsm[0,:].tolist().index(sel_index)

    return ms_list[smsm[0,col_index]]


def get_smsm(ms,
             n_sample_size: int,
             sample_kind: str,
             sort_by: Union[str, None] = None,
             reverse = False,
             seed = None,
             save_to_npz = False) -> np.ndarray:
    """
    Return the Sampled Microstates State Matrix (smsm) for free residues as part of a tuple:
    (info, selection_energies, smsm). If `save_to_npz` is True, the same tuple is saved into
    a numpy `.npz` file containing these three items.
    """

    print(f"Sampling microstates (kind: {sample_kind}) for MC{ms.selected_MC}.")

    ms_list, ms_indices, info = ms.get_sampling_params(n_sample_size,
                                                         kind = sample_kind,
                                                         sort_by = sort_by,
                                                         reverse = reverse,
                                                         seed = seed
                                                         )
    k = sample_kind[0].lower()
    s = sort_by[0].lower()
    if sample_kind.lower() == "random":
        s = ""
        ms_list = ms.microstates

    selection_energies = []
    top_rows = 2
    smsm = np.ones((len(ms.free_residues) + top_rows, n_sample_size), dtype=int) * -1

    for i, x in enumerate(ms_indices):
        sampled_ms = ms_list[x]
        selected_iconfs = ms.get_selected_confs(sampled_ms,
                                                output_val="iconf",
                                                include_fixed = False)
        smsm[0, i] = x                 # row 0 :: ms selection index
        smsm[1, i] = sampled_ms.idx    # row 1 :: selected_ms.idx
        selection_energies.append(sampled_ms.E)

        for r, iconf in enumerate(selected_iconfs, start=2):
            smsm[r, i] = iconf

    info.append(f"{smsm.shape}")
    info.append(datetime.today().strftime("%d-%b-%y %H:%M:%S"))
    if save_to_npz:
        npz_file = ms.msout_file_dir.joinpath(f"smsm{ms.selected_MC}_{k}{s}.npz")
        info.append(npz_file)

        if npz_file.exists():
            npz_file.unlink()

        save_npz(npz_file, info = info, sel_energies = selection_energies, smsm=smsm)

    return  info, selection_energies, smsm
<!-- #endraw -->

```python
ms.free_residue_names[:10]
```

```python
n_sample_size = 4
sample_kind =  "deterministic"
ms_sort_by = "energy"
save = True

info, selection_energies, smsm = sampling.get_smsm(ms, n_sample_size, sample_kind, ms_sort_by, save_to_npz=save)
info
```

```python
fname = sampling.get_output_filename(ms.selected_MC,
                                     sample_kind,
                                     ms_sort_by,
                                     False,
                                     size=n_sample_size)
d_npz_file = ms.msout_file_dir.joinpath(fname)
smsm_data = load_npz(d_npz_file)
info_dict = sampling.smsm_data_info_to_dict(smsm_data["info"])
info_dict
```

```python
len("MET01A0001_002")
xx = 467.14
repr(xx)
sx = f"{xx:,.2f}"
sx
```

```python
sel_indices = smsm_data["selection_indices"]
sel_energies = smsm_data["selection_energies"]
res_names = smsm_data["residue_names"]
smsm = smsm_data["smsm"]

#hdr = f"RESIDUE   | " + "".join([f"{x:<15}" for x in sel_indices])
#hdr = f"RESIDUE   | " + " ".join([f"{x:>14,.2f}" for x in sel_energies])
#print(hdr)
#print("_" * len(hdr))

for r in range(smsm.shape[0]):
    rc = smsm[r,:]
    row_confs = ms.confnames_by_iconfs(rc.tolist())

    print(f"{res_names[r]} |", *row_confs)
    if r > 5:
       break

```

```python
smsm[:5,:5]
res_names[:5]
```

# RESUME HERE

```python
def smsm_to_full_confids(smsm_data):
    """Convert smsm matrix numbers (iconfs) to confIDs and add the
    fixed residues in each selected state.
    """

    sel_indices = smsm_data["selection_indices"]
    sel_energies = smsm_data["selection_energies"]
    res_names = smsm_data["residue_names"]
    smsm = smsm_data["smsm"]

    #hdr = f"RESIDUE   | " + "".join([f"{x:<15}" for x in sel_indices])
    #hdr = f"RESIDUE   | " + " ".join([f"{x:>14,.2f}" for x in sel_energies])
    #print(hdr)
    #print("_" * len(hdr))

    for r in range(smsm.shape[0]):
        rc = smsm[r,:]
        row_confs = ms.confnames_by_iconfs(rc.tolist())
        print(f"{res_names[r]} |", *row_confs)


    for i, c in enumerate(range(smsm.shape[1])):
        rc = smsm[:,i]

        # add fixed res
        all_iconfs = sorted(ms.fixed_iconfs + rc)
        selected_confs = ms.confnames_by_iconfs(all_iconfs)
        print(f"{len(selected_confs)= }")


        if c == 0:
            break
        # write pdb with one_pdb_confs; ms_index from smsm row 1
        #ms_to_pdb(selected_confs: list, ms_index: int, mc_run: int, remark_data: str, step2_path: str, output_folder: str)

```

<!-- #raw -->
# selected MC run:
MC = int(info[1].split("=")[1])

# part of the info related to sampling:
info_0 = info[0].split(", ")
print(info_0)

size = int(info_0[0].split("=")[1])
kind = info_0[1].split("=")[1][1:-1]
by = info_0[2].split('=')[1][1:-1]
reverse = False if info_0[3].split("=")[1].startswith("F") else True

reverse

if v := info_0[4].split("=")[1] == "None":
    seed = None
else:
    seed = int(v)
print(v, seed)
<!-- #endraw -->

```python
!ls -l {ms.msout_file_dir}
```

```python
n_sample_size = 4
sample_kind =  "deterministic"
ms_sort_by = "energy"
save = True
out_dir = Path.cwd()

sampling.pdbs_from_smsm(ms,
                        n_sample_size,
                        sample_kind,
                        sort_by = ms_sort_by,
                        sort_reverse = False,
                        seed = None,
                        output_pdb_format = "standard",
                        output_dir = out_dir,
                        clear_pdbs_folder = False,
                        list_files = True
                       )

```

```python
# Inspect a pdb head:

!head -n 50 ./pdbs_from_ms/mc0_ms104828.pdb
```

<!-- #raw -->
assert len(ms.fixed_residue_names) == len(ms.fixed_iconfs)

for i, (res,conf) in enumerate(zip(ms.fixed_residue_names,ms.fixed_iconfs)):
    print(i, res, conf)
    if i > 10:
        break

for i, (res,confs) in enumerate(zip(ms.free_residue_names, ms.free_residues)):
    print(i, res, confs)
    if i > 10:
        break
<!-- #endraw -->

```python

```

```python

```

```python

```

---
# TODO? PDB conversions

<!-- #region -->
## name.txt
> Rename rule file contains rules to format atom names and residue names.
 The purpose is to unify residue names to 3-char mcce names
 and break some big cofactors into small ones.
 Each line has two fields separated by at least one space. Each field
 is 14 characters long, matching the atom, residue name, chainID and
 sequence number field of a pdb line. The first string will be
 replaced by the second string.
 Symbol "*" in the first string is a wildcard that matchs any character.
 It means "do not replace" in the second string.
 The replace is accumulative in the order of appearing in this file.


The first 3 lines:
```python
0123456789012  01234567890123
***** *******  *****_********
****** *****   ******_*******
******* ***  > *******_****
```
mean: replace a space in position 5, 6 or 7 by "_".

**
<!-- #endregion -->

```python
nametxt = mcce_dir.joinpath("name.txt")
nametxt

# use when running subprocesses from this nbk
nametxt2 = nametxt.relative_to(Path.cwd().parent)
nametxt2
```

```python
import re
from collections import defaultdict


regex = r'^\*{5,8}\s\*{5,}\s\s'
split_tag = "*@@"

mcce_to_pdb_name_dict = defaultdict(list)
with open(nametxt) as f:
    for i, line in enumerate(f):
        if (line.startswith("#")
            or line.startswith(" D")
            or re.match(regex, line)
           ):
            continue
        if len(line) <29:
            continue

        line = line[:30].replace("*  ", split_tag)
        from_str, to_str = line.split(split_tag)[:2]
        mcce_to_pdb_name_dict[to_str].append(from_str)

```

```python
#pp(mcce_to_pdb_name_dict)
```

```python
s = "***** NA******  *****_NA******"
x = "012345678901234567890123456789"
len(s), len(x)
```

<!-- #raw -->
gr_data = subprocess.check_output(
                f"grep -E 'Gehan' ../{nametxt2}",  # | sed -e 's/(//g; s/MONTE_//g; s/)//g'",
                stderr=subprocess.STDOUT,
                shell=True,).decode().splitlines()

gr_data
<!-- #endraw -->

<!-- #raw -->
# in pymccelib.py

# env = Env()
def pdb2mcce(self, pdb):
        """Convert pdb to mcce pdb"""
        atom_exceptions = [" H2 ", " OXT", " HXT"]
        mccelines = []
        lines = [x for x in open(pdb).readlines() if x[:6] == "ATOM  " or x[:6] == "HETATM"]

        icount = 0
        previous_resid = ()
        possible_confs = []
        for line in lines:
            # pdb line
            atomname = line[12:16]
            resname = line[17:20]
            chainid = line[21]
            seqnum = int(line[22:26])
            icode = line[26]
            xyz = line[30:54]

            current_resid = (resname, chainid, seqnum, icode)
            # mcce line, need to add conf_number, radius, charge, conf_type, conf_history
            if current_resid != previous_resid:
                possible_confs = [x.strip() for x in env.tpl[("CONFLIST", resname)].split(",")]
                logging.info("Identified a new residue %s: %s" % (resname, ", ".join(possible_confs)))
                previous_resid = current_resid
            Found = False
            for confname in possible_confs:
                if atomname in env.atomnames[confname]:
                    conf_type = confname[3:5]
                    conf_number = possible_confs.index(confname)
                    cname = confname
                    Found = True
                    break
            if not Found:
                # this atom is not found in all conformers
                if atomname not in atom_exceptions:
                    print("Atom \"%s\" in pdb file %s can not be assigned to any conformer" % (atomname, pdb))
                continue

            key = ("RADIUS", cname, atomname)
            if key in env.tpl:
                radius_str = env.tpl[key]
                rad, _, _ = radius_str.split(",")
                rad = float(rad)
            else:
                rad = 0.0

            key = ("CHARGE", cname, atomname)
            if key in env.tpl:
                charge_str = env.tpl[key]
                crg = float(charge_str)
            else:
                crg = 0.0

            conf_history = "________"
            newline = "ATOM  %5d %4s %s %c%4d%c%03d%s%8.3f    %8.3f      %s%s\n" % \
                      (icount, atomname, resname, chainid, seqnum, icode, conf_number, xyz, rad, crg, conf_type, conf_history)
            mccelines.append(newline)
            icount += 1

        return mccelines

<!-- #endraw -->
