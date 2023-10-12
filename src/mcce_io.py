"""
mcce_io: Module contains a collection of functions to find/access mcce output files.
Module to be used in conjuction with ms_analysis.py.
"""

from pathlib import path
import subprocess
import sys
import platform


def get_mcce():
    """Locate mcce executable."""

    sys = platform.system()
    # >>'Windows' | 'Linux' | 'Darwin' (Mac)
    cmd = "which"
    if sys == "Windows":
        cmd = "where"

    # Run the command
    result = subprocess.run([cmd, "mcce"], stdout=subprocess.PIPE)

    # Decode the output from bytes to a string
    output = result.stdout.decode("utf-8")

    # Print the output
    print(output)

    # is code running in an env?
    base_env_prefix = sys.base_prefix
    ex_prefix = sys.exec_prefix
    print(base_env_prefix, ex_prefix)
