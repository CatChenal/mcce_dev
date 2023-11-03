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

_Cell 1: from jupyterlab template_: run it

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

_Cell 2: from jupyterlab template_: run it

```python
# Insert current src dir into sys.path so that modules in ../src can be imported:
# CHANGE THIS IF NEEDED:

add_to_sys_path(Path.cwd(), up=True)
```

```python
DATA = Path.cwd().parent.joinpath("data")
if DATA.exists():
    DATA
else:
    print(f"Not found: {DATA}")
```

```python
mcce_dir = DATA.joinpath("4lzt")
mcce_dir
msout_dir = mcce_dir.joinpath("ms_out")
!ls -l {msout_dir}

runprm = mcce_dir.joinpath("run.prm.record")
runprm2 = runprm.relative_to(Path.cwd().parent)  # use when running from this nbk
```

```python
import ms_sampling_to_pdbs as sampling
import subprocess
```

```python
#fdir(sampling)
```

---
# Test executable module

```python
which = msout_dir
which = msout_dir.joinpath("pH6eH0ms")
!ls -l {which}
```

```python
mcce_dir = DATA.joinpath("4lzt")
mcce_dir
pH = 6.0
Eh= 0.0
```

```python
ms = sampling.MS(mcce_dir,pH, Eh)
```

```python
#fdir(ms)
```

```python
n_ires = len(ms.ires_by_iconf)
print("n_ires = len(ms.ires_by_iconf) = ", n_ires)
n_confs = len(ms.conformers)
print("n_confs = len(ms.conformers) = ", n_confs)
n_fixed = len(ms.fixed_iconfs)
print("n_fixed = len(ms.fixed_iconfs) = ", n_fixed)
free_confs = n_confs - n_fixed
print("free_confs = n_confs - n_fixed = ", free_confs)
n_free = len(ms.free_residues)
print("n_free = len(ms.free_residues) = ", n_free)

print(f"{ms.MC_RUNS = }, {ms.MC_NITER = :,}")
print(f"{ms.counts = :,}")

print(f"{ms.MC_NITER * ms.MC_RUNS * n_confs  = :,}")
print(f"{ms.MC_NITER * ms.MC_RUNS * free_confs = :,}")
print(f"{ms.MC_NITER * ms.MC_RUNS * n_ires = :,}")
print(f"{ms.MC_NITER * ms.MC_RUNS * n_fixed = :,}")
print(f"{ms.MC_NITER * ms.MC_RUNS * n_free = :,}")
```

```python
n_counts = 0.0
ms_count_values = []

ms_list = ms.sort_microstates(by="energy")

for mc in ms_list:
    # n_counts += ms[1]
    # ms_count_values.append(ms[1])
    n_counts += mc.count
    ms_count_values.append(mc.count)

n_counts = ms.counts
ms_cumsum = np.cumsum([mc.count for mc in ms_list])

ms.counts
n_counts
ms_cumsum[-1]
ms_cumsum2[-1]
len(ms.microstates)
len(ms_count_values)
```

```python
def sample_microstates(
    ms: MS,
    size: int,
    kind: str = "deterministic",
    sort_by: str = "energy",
    reverse: bool = False,
) -> tuple:
    """
    Implement a sampling of MS.microstates depending on `kind`.
    Args:
        ms (MS): An instance of the MS class.
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
                raise ValueError(f"Values for `sort_by` are 'energy' or 'count'; Given: {by}")

        ms_list = ms.sort_microstates(by=sort_by, reverse=reverse)
    else:
        ms_list = ms.microstates

    ms_cumsum = np.cumsum([mc.count for mc in ms_list])

    if kind == "deterministic":
        X = ms.counts - size
        Y = ms.counts / size
        count_selection = np.arange(size, X, Y)
    else:
        rng = np.random.default_rng()
        count_selection = rng.integers(low=1, high=n_counts + 1, size=size)

    return ms_cumsum, count_selection, ms_list
```

```python

```

```python

```

# Getters and Setters

<!-- #raw jupyter={"source_hidden": true} -->
# this works; __eq__ and _hash_ needed for lru_cache
from functools import lru_cache

class Point:

    def __init__(self, x, y, MC = 0):
        self.x = x
        self.y = y
        self.RUNS = 6
        self._selected_MC = MC
        self.do_something_with_MC()

    @property
    def selected_MC(self):
        return self._selected_MC

    @selected_MC.setter
    def selected_MC(self, value):
        if value >= self.RUNS:
             print(f"This value is beyond the range of MONTE_RUNS({RUNS}) used in this simulation: {value}.")
        else:
            self._selected_MC = value
            self.do_something_with_MC()

    def __eq__(self, other):
        return self.selected_MC == other.selected_MC

    def __hash__(self):
        return hash(self.selected_MC)

    def __repr__(self):
        return f"Point({self.x=}, {self.y=}, {self.selected_MC=})"

    @lru_cache(maxsize=None)
    def do_something_with_MC(self):
        print(f"do_something with {self.selected_MC=}")


print("init")
p = Point(1,2)
p

print("selected_MC -> 5")
p.selected_MC = 5
p

print("selected_MC -> 6")
p.selected_MC = 6
p

print("selected_MC -> 1")
p.selected_MC = 1
p
<!-- #endraw -->

<!-- #raw jupyter={"source_hidden": true} -->
# this works

RUNS = 6

class Point:
    x_values = set()
    new_x = False

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __getattr__(self, name: str):
        return self.__dict__[f"_{name}"]

    def __setattr__(self, name, value):
        if name == "x":
            if value in range(RUNS):
                if value not in Point.x_values:
                    Point.new_x = True
                else:
                    Point.new_x = False
                Point.x_values.add(value)
            else:
                print(f"This value is beyond the range of MONTE_RUNS({RUNS}) used in this simulation: {value}.")

        self.__dict__[f"_{name}"] = float(value)

    def __repr__(self):
        return f"Point({self.x=}, {self.y=})"


p = Point(1,2)

p.x, p.y
p.x_values
p.__dict__

p.x = 6
p.y = 10

p.x, p.y
p.x_values
p.__dict__
<!-- #endraw -->

<!-- #raw jupyter={"source_hidden": true} -->
# this works too but unwielding: must call a get_/set_x fn for every assignment

class Person:
    def __init__(self, name, birth_date):
        self._name = name
        self._birth_date = birth_date
        self._changes = []

    def get_name(self):
        return self._name

    def set_name(self, value):
        if self._name != value:
            self._name = value
            self._changes.append(("name",value))
            self.do_this(self._name)

    def get_changes(self):
        return self._changes

    def get_birth_date(self):
        return self._birth_date

    def set_birth_date(self, value, force=False):
        if self._birth_date != value:
            if not force:
                raise AttributeError("can't set birth_date")
            self._birth_date = value
            self._changes.append(("birth_date",value))
            self.do_this(self._birth_date)

    def do_this(self, attrname):
        print(f"I'm doing this for {attrname}")

class Employee(Person):
    def get_name(self):
        return super().get_name().upper()

    def get_birth_date(self):
        return super().get_birth_date()

    def get_changes(self):
        return super()._changes

# examples

jane = Person("Jane Doe", "2000-11-29")
jane.get_name()
jane.get_birth_date()

jane.set_birth_date("2000-10-29", force=True)
jane.get_birth_date()
jane.get_changes()

jane.set_name("Jane Poe")
jane.get_name()
jane.get_changes()

jane.set_birth_date("2000-10-29", force=True)
jane.get_birth_date()
jane.get_changes()
<!-- #endraw -->

---

---

# MCCE - MS Sampling (using test data in ../tests/data/)
---

```python
import base
import mcce_io as io
import ms_sampling as sampling
```

```python
DATA = Path.cwd().parent.joinpath("tests/data")
DATA
```

```python
!ls {DATA}
```

```python
# filepaths of inputs used by MS class:
h3_path = DATA.joinpath("head3.lst")
mcce_output_path = h3_path.parent
mcce_output_path

step2_path = mcce_output_path.joinpath("step2_out.pdb")
msout_path = mcce_output_path.joinpath("ms_out")
msout_path, msout_path.is_dir()
```

```python
opath = None
(not opath) or opath is None
```

```python
# filepaths of outputs:

pH = 5.0
Eh= 0.0
msout_file = io.get_msout_filename(mcce_output_path, pH, Eh)
msout_file

msout_file_dir = msout_file.parent.joinpath(msout_file.stem)
msout_file_dir
```

```python
start_time = time.time()

io.split_msout_file(mcce_output_path, pH, Eh)

end_time = time.time()
print("io.divide_msout_file() took {:.2f} mins".format((end_time - start_time)/60))
```

```python
!ls {msout_file_dir}
```

```python
pdbs_dir = msout_file_dir.joinpath("pdbs_from_ms")
```

```python
!ls {pdbs_dir}
```

<!-- #raw -->
io.clear_folder(pdbs_dir)
!ls {pdbs_dir}
<!-- #endraw -->

# base.MC class

```python
print(base.MS.__doc__)
print(base.MS.__init__.__doc__)
```

```python
# create instance
start_time = time.time()

ms = base.MS(mcce_output_path, pH, Eh)

d = time.time() - start_time
print(f"Loading of base.MS instance took {d/60:.2f} mins or {d:.2f} seconds")
print(ms)
```

```python
# Public vars in MC:
fdir(ms)
```

```python
print(f"ms.counts : {ms.counts:,}")
```

# ms_sampling module

```python
fdir(sampling)
```

```python
n_sample_size = 4
ms_sort_by = "energy"
output_dir = msout_file_dir
mc_run = ms.selected_MC  # part of pdb name
mc_run
```

<!-- #raw jupyter={"source_hidden": true} -->
# check on diff out with diff sort roder: ok
ms_sort_by = "energy"
e_sorted_ms_list = sampling.sort_microstate_list(ms.microstates, by=ms_sort_by)
e_ms_cumsum, e_count_selection = sampling.sample_microstates(n_sample_size, e_sorted_ms_list)
len(e_ms_cumsum), len(e_count_selection)
e_count_selection
assert max(e_ms_cumsum) == ms.counts


ms_sort_by = "count"
c_sorted_ms_list = sampling.sort_microstate_list(ms.microstates, by=ms_sort_by)
c_ms_cumsum, c_count_selection = sampling.sample_microstates(n_sample_size, c_sorted_ms_list)
len(c_ms_cumsum), len(c_count_selection)
c_count_selection


e_ms_cumsum[:5]
e_sorted_ms_list[:5]
c_ms_cumsum[:5]
c_sorted_ms_list[:5]

def check_confs_byorder(sorted_ms_list, count_selection, ms_cumsum):
    for c in count_selection[:3]:
        print("c:",c)
        ms_index = np.where((ms_cumsum - c) > 0)[0][0]
        ms_selection = sorted_ms_list[ms_index]
        print(ms_selection[0])
        print(ms_selection[2]())
        confs_for_pdb = sampling.get_selected_confs(ms, ms_selection)
        print(confs_for_pdb[:5])
        print("call io.ms_to_pdb")


check_confs_byorder(e_sorted_ms_list, e_count_selection, e_ms_cumsum)
check_confs_byorder(c_sorted_ms_list, c_count_selection, c_ms_cumsum)
<!-- #endraw -->

```python
# create pdbs from samples ms
start_time = time.time()

sampling.pdbs_from_ms_samples(ms,
                              mcce_output_path,
                              n_sample_size,
                              ms_sort_by,
                              output_dir,
                              list_files=True)

d = time.time() - start_time
print(f"`sampling.pdbs_from_ms_samples` with sample size={n_sample_size:,} took {d/60:.2f} mins or {d:.2f} seconds")
```

```python
!head -n 20 ../tests/data/ms_out/pH5eH0ms/pdbs_from_ms/mc0_ms1.pdb
```

```python

```
