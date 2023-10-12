import sys
from pathlib import Path


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
        print("Path added to sys.path: {}".format(newp))


# if notebook inside another folder, eg ./notebooks:
nb_folder = "notebooks"
add_to_sys_path(Path.cwd(), Path.cwd().name.startswith(nb_folder))


def get_project_dirs(which=["data", "images"], nb_folder="notebooks"):
    dir_lst = []
    if Path.cwd().name.startswith(nb_folder):
        dir_fn = Path.cwd().parent.joinpath
    else:
        dir_fn = Path.cwd().joinpath

    for d in which:
        DIR = dir_fn(d)
        if not DIR.exists():
            Path.mkdir(DIR)
        dir_lst.append(DIR)
    return dir_lst


# DIR_DATA, DIR_IMG = get_project_dirs()

import numpy as np
import scipy as sp
from scipy import stats as sps
import pandas as pd

# pd.set_option("display.max_colwidth", 200)

import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sb

# plt.ion()
# %matplotlib inline
plt.style.use("seaborn-v0_8-muted")

from pprint import pprint as pp


# Filtered dir() for method discovery:
def fdir(obj, start_with_str="_", exclude=True):
    return [d for d in dir(obj) if not d.startswith(start_with_str) == exclude]


from IPython.core.interactiveshell import InteractiveShell

InteractiveShell.ast_node_interactivity = "all"

from IPython.display import HTML, Markdown  # , IFrame

# for presentations:
# display(HTML("<style>.container { width:100% !important; }</style>"))

# autoreload extension
from IPython import get_ipython

ipython = get_ipython()


def show_versions():
    txt = "<pre><br>"
    txt += f"Python:\t\t{sys.version}<br>"
    txt += f"Python env:\t{Path(sys.prefix).name}<br>"
    txt += f"Numpy:\t\t{np.__version__}<br>"
    txt += f"Scipy:\t\t{sp.__version__}<br>"
    txt += f"Pandas:\t\t{pd.__version__}<br>"
    txt += f"Matplotlib:\t{mpl.__version__}<br>"
    txt += f"Currrent dir: {Path.cwd()}"
    txt += "</pre>"
    div = f"""<div class="alert alert-info"><b>Versions:</b><br>{txt}</div>"""
    return HTML(div)


# All % lines commented for pylint
# if 'autoreload' not in ipython.extension_manager.loaded:
#    %load_ext autoreload

# %autoreload 2
# ..................

no_wmark = False
""" pylint
try:
    %load_ext watermark
    %watermark
except ModuleNotFoundError:
    no_wmark = True

if no_wmark:
    show_versions()
else:
    %watermark -iv
"""
