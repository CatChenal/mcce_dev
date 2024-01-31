"""
Hack to download a (recently updated) /bin/<module> from github GunnerLab/Stable-MCCE
repo as there is no updated conda package on channel 'newbooks' as of 12-10-2023.
Usage:
 1. If <module> already exists in target location, rename it or delete it
 2. Run `python update_mcce_module.py <module> [<target location>]
"""

import argparse
import importlib
from pathlib import Path
import sys


def update_module(module_name: str, output_dir: str = None):
    """
    Update/create a local python module with one downloaded from the github
    repo GunnerLab/Stable-MCCE/bin folder.

    module_name (str): Name of python module file, e.g.: module.py
    output_dir (str): Folder where to save the downloaded module. Current
                      working directory is used if not provided.
    """

    modul = Path(module_name)
    if not modul.suffix:
        modul = Path(modul_name + ".py")

    if modul.suffix != ".py":
        raise ValueError(f"Invalid extension: {modul.suffix}, expected: '.py'")

    help = "Remove or rename it, then run `python update_mcce_module.py <name>` again."

    mdl = modul.stem
    try:
        importlib.import_module(mdl)
        print(f"Import of {mdl} was successful. If you meant to update it:\n{help}")

    except ModuleNotFoundError:
        mdl_name = modul.name

        if output_dir is None:
            output_path = Path(__file__).parent.joinpath(mdl_name)
        else:
            output_path = Path(output_dir).joinpath(mdl_name)

        if output_path.exists():
            print("Module already exists.", help)
            return

        repo_bin = "https://raw.githubusercontent.com/GunnerLab/Stable-MCCE/master/bin/"
        modul_url = repo_bin + mdl_name
        curl_cmd = f"curl {modul_url} -o {str(output_path)}".split()

        import subprocess

        try:
            subprocess.run(curl_cmd, check=True)
        except subprocess.CalledProcessError as cpe:
            print(
                f"Could not download {modul.name} from 'GunnerLab/Stable-MCCE'.\nError:{cpe}"
            )
            raise


if __name__ == "__main__":
    desc = (
        "Download a (recently updated) <module> from github GunnerLab/Stable-MCCE/bin."
    )
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("module_name", help="Name to the python module to download.")
    parser.add_argument(
        "-target_dir", default="", help="Folder path for the downloaded module."
    )
    args = parser.parse_args()

    update_module(args.module_name, args.target_dir)
