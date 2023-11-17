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

# TUTORIAL: Discovering microstates using `ms_sampling_to_pdbs` executable module
### This tutorial demonstrates the programmatic usage of various functions when the module is imported.
### For usage of the module as a command line interface, see the README.md file. 
## Note:
### If your ms_out folder is zipped, you need to unzip it.
---

__RUN THESE FIRST 3 CELLS without modifications__:

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

```python
# Insert current src dir into sys.path so that modules in ../src can be imported:
add_to_sys_path(Path.cwd(), up=True)
```

```python
from functools import partial
import ms_sampling_to_pdbs as sampling

load_npz = partial(np.load, allow_pickle=True)
save_npz = partial(np.savez, allow_pickle=True)
```

---
---
# CONTENTS

## Paths definitions
## Test executable module
## Info from header and MC files (there is no microstate data in the ms instance, yet)
### _header_ file
### MCx file
## Obtaining the microstates data
## What's in MCx.npz?
## Create the "Sampled Microstates State Matrix" data (sms_data in .npz file)
## What's in a smsmx_*.npz file?
## Matrix from a random sampling
### info contents
### sel_energies (i.e. selection energies) contents
### smsm: the matrix
### Using the top 2 rows of the matrix
## Matrix from a deterministic sampling
### info contents
### sel_energies (i.e. selection energies) contents
### smsm: the matrix
### Getting a microstste object from the matrix
## Returning data when calling `sampling.get_smsm`
## Create a pdb collection from smsm data: use `sampling.pdbs_from_smsm`
## Inspect a pdb head
---
---


---
## Paths definitions

```python
HERE = Path.cwd()  # do not change

test_folder = "granepura_GunnerLab/mcce_data" # you can change this to any other mcce folder
mcce_dir = HERE.parent.parent.joinpath(test_folder)  # do not change
mcce_dir, mcce_dir.exists()
```

```python
!ls -l {mcce_dir}
```

```python
sampling.check_mcce_dir(mcce_dir)  # should run error-free
```

```python
msout_dir = mcce_dir.joinpath("ms_out")

!ls -l {msout_dir}
```

---
## Test executable module

```python
# Need to change pH, Eh according to the file(s) listed in msout_dir
pH = 7
Eh= 0
overwrite = False

ms = sampling.MS(mcce_dir, pH, Eh, overwrite_split_files=overwrite)
print(ms)  # show how MS was instantiated
```

```python
! ls -l {ms.msout_file_dir}
```

---
## Info from header and MC files (there is no microstate data in the ms instance, yet):

```python
print(f"{len(ms.fixed_residue_names) = }; {len(ms.fixed_iconfs) = };\n{len(ms.free_residue_names) = }\n")
print(f"{ms.counts = }, should be 0: {ms.counts==0}\n{ms.selected_MC = }, should be None: {ms.selected_MC is None}")
```

## header file:

```python
hdr = ms.msout_file_dir.joinpath("header")

print("The 3rd line lists the fixed confs\nThe 4th line lists the free residues confs\n")
! head -n4 {hdr}
```

## MC0 file:

```python
mc0 = ms.msout_file_dir.joinpath("MC0")

print("The first line lists the free residues confs selected in the initial state of the MC0 run:\n")
! head -n3 {mc0}
```

---
## Obtaining the microstates data:

```python
# data for MC index = idx
idx = 0
ms.get_mc_data(idx)

print(f"{ms.counts = :,}; {ms.MC_RUNS = :,}")
print(f"{len(ms.microstates) = :,}")
```

```python
# Notice the MCx.npz file created:
!ls -l {ms.msout_file_dir}
```

## What's in MC0.npz?

```python
MC0 = load_npz(ms.msout_file_dir.joinpath("MC0.npz"))
# The files in a MCx.npz file can be accessed like a dictionnary
MC0.files   # ignore 'allow_pickle'

counts = MC0["counts"]           # -> used to populate ms.counts
microstates = MC0["microstates"] # -> np.array used to populate ms.microstates after conversion to list
print(f"{counts = :,}, {len(microstates) = :,}")

# ms.microstates holds a list of Microstates objects:
ms.microstates[:5]
```

```python
# Here are the attibutes of the first microstate:
mc0 = ms.microstates[0]
print(f"{mc0.E = :,.2f}, {mc0.count = :,}\n{mc0.state() = }")
```

---
## Create the "Sampled Microstates State Matrix" data (sms_data in .npz file)

```python
print(sampling.get_smsm.__doc__)
```

## What's in a smsmx_*.npz file?
### info         :: info needed to recover the original microstate
### sel_energies :: the selection energies (floats) of the selected microstates; saved separately to keep the matrix single typed (integers).
### smsm         :: the state matrix



---
## Matrix from a random sampling:

```python
n_sample_size = 3
sample_kind = "random"

# Here, `only_save`=True means create and save the file.
# We're going to reload it shortly.
sampling.get_smsm(ms, n_sample_size, sample_kind, only_save=True)

# file created:
smsm_r_file = sampling.get_output_filename(ms.selected_MC, sample_kind, "", False, size=n_sample_size)
smsm_r_file
```

```python
# Notice the smsmx_*.npz file created:
!ls -l {ms.msout_file_dir}
```

```python
smsm_r_data = load_npz(ms.msout_file_dir.joinpath(smsm_r_file))
smsm_r_data.files  # ignore 'allow_pickle'
```

### info contents:

```python
info_r = smsm_r_data["info"]
len(info_r)
info_r
```

```python
# Here is a helper function to return a dict from the info:
info_dict = sampling.smsm_data_info_to_dict(smsm_r_data["info"])
pp(info_dict)
```

```python
# matrix shape:
r, c = info_dict["smsm_shape"]
r, c

# Check: the number of rows (r) should equal len(ms.free_residues +2) (more on the +2 shortly)
r == len(ms.free_residues) + 2
```

### sel_energies (i.e. selection energies) contents:

```python
energies_r = smsm_r_data["sel_energies"]
energies_r.shape
energies_r  # with ridiculous precision!
```

### smsm: the matrix

```python
smsm0r = smsm_r_data["smsm"]  # a np array
smsm0r.shape, smsm0r.dtype

smsm0r
```

---
## Using the top 2 rows of the matrix:
### First row: the selection index: depends on the kind of sampling
### Second row: the microstate index of creation

The original microstate object can be "recover" with this function:

```python
col_index = 0

original_ms = sampling.get_ms_from_smsm(ms, smsm_r_data, col_index)
print(original_ms)

# Check [1,col_index] == original_ms.idx:
smsm0r[1,col_index] == original_ms.idx
```

#### Note: when sampling kind is random, the original ms can be recovered directly from ms.microstates:

```python
i = smsm0r[0,col_index]
m = ms.microstates[i]
print(m)
```

---
## Matrix from a deterministic sampling:

```python
n_sample_size = 3
sample_kind = "deterministic"
ms_sort_by = "energy"

# Here, `only_save`=True means create and save the file.
# We're going to reload it shortly.
sampling.get_smsm(ms, n_sample_size, sample_kind,ms_sort_by, only_save=True)

# file created:
smsm_d_file = sampling.get_output_filename(ms.selected_MC, sample_kind,  ms_sort_by, False, size=n_sample_size)
smsm_d_file
```

```python
!ls -l {ms.msout_file_dir}
```

```python
smsm_d_data = load_npz(ms.msout_file_dir.joinpath(smsm_d_file))
smsm_d_data.files  # ignore 'allow_pickle'
```

### info contents:

```python
info_d = smsm_d_data["info"]

info_dict_d = sampling.smsm_data_info_to_dict(smsm_d_data["info"])
pp(info_dict_d)
```

### sel_energies (i.e. selection energies) contents:

```python
energies_d = smsm_r_data["sel_energies"]
energies_d.shape
energies_d
```

### smsm: the matrix

```python
smsm0d = smsm_d_data["smsm"]
smsm0d.shape, smsm0d.dtype

smsm0d
```

### Getting a microstate object from the matrix:

```python
col_index = 0

original_ms = sampling.get_ms_from_smsm(ms, smsm_d_data, col_index)
print(original_ms)

# Check [1,col_index] == original_ms.idx:
smsm0d[1,col_index] == original_ms.idx
```

### No direct access to a microstate from ms.microstates if sampling is deterministic:

```python
i = smsm0d[0,col_index]
m = ms.microstates[i]
print(m)

smsm0d[1,col_index] == m.idx
```

---
## Returning data when calling `sampling.get_smsm`
### Simply omit the `only_save` argument (default is False)

```python
n_sample_size = 3
sample_kind =  "deterministic"
ms_sort_by = "energy"
save = True   # up to you

info, selection_energies, smsm = sampling.get_smsm(ms, n_sample_size, sample_kind, ms_sort_by, save_to_npz=save)
```

```python
info_dict = sampling.smsm_data_info_to_dict(info)
pp(info_dict)
```

---
## Create a pdb collection from smsm data: use `sampling.pdbs_from_smsm`

## Note:
#### The `output_dir` can be set to anything: the actual output location wil be `output_dir/pdbs_from_ms`.
#### Yet, it is recommend not to pass it at all (the default is None) then, the pdbs are created "where they belong", i.e. in `/ms_out/pH7eH0ms/pdbs_from_ms`.

The default location can be accessed using the instance:
```
pdb_folder = ms.msout_file_dir.joinpath("pdbs_from_ms")
```

```python
n_sample_size = 3
sample_kind =  "deterministic"
ms_sort_by = "count" #"energy"
save = True
#output_dir = Path.cwd()   # default output location if commented

sampling.pdbs_from_smsm(ms,
                        n_sample_size,
                        sample_kind,
                        sort_by = ms_sort_by,
                        sort_reverse = False,
                        seed = None,
                        output_pdb_format = "standard",
                        #output_dir = output_dir,
                        clear_pdbs_folder = False,
                        list_files = True
                       )

```

```python
pdb_folder = ms.msout_file_dir.joinpath("pdbs_from_ms")

!ls -l {pdb_folder}
```

```python
sel_index = 452162  # pickone that exists!

pdb_name = sampling.get_output_filename(ms.selected_MC, sample_kind,  ms_sort_by, False, sel_index=sel_index, for_pdb=True)
pdb_path = pdb_folder.joinpath(pdb_name)
pdb_path
```

---
## Inspect a pdb head:

```python
!head -n 50 {pdb_path}
```

---
---
