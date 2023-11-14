# mcce_dev

Project for ancillary development of [Stable-MCCE repo](https://github.com/GunnerLab/Stable-MCCE).

# Current projects:

## ms_sampling_to_pdb
### -> Provide an executable module with cli to obtain a collection of pdb files created from sampled microstates.
### Platforms: Linux and Linux-like

* Module __doc__ string:
```
    MODULE `ms_sampling_to_pdbs.py` is a command-line executable module that creates the sampled MCCE
    microstates state matrix (smsm) file for each selected Monte Carlo runs, and optionally writes a
    collection of pdb files from the smsm data.

    Note: The smsm data (.npz) file name has this format:
          f"smsm{MC}_{size}_{k}{s}[_rev].npz"

          - MC     index of the selected Monte Carlo run
          - size   the sample size
          - k      initial of the sampling kind, either "d" or "r"
          - s      sort key (for deterministic sampling), either "e" or "c" (energy or count)
          - [_rev] optional indicator if order is reversed (descending)

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

    USAGE:
    -----
    The minimal number of arguments are all the positional arguments: mcce_dir, pH, Eh, sample_size,
    i.e.:
    > ms_sampling_to_pdbs.py a/path/to/mcce/output 7 0 99

    All other arguments have default values, which can be changed by using "-<arg name>" followed by a value,
    i.e.:
    > ms_sampling_to_pdbs.py a/path/to/mcce/output 7 0 99 -MC 2,3          # use Monte Carlo runs 2 and 3; default 0
    > ms_sampling_to_pdbs.py a/path/to/mcce/output 7 0 99 -sampling_kind r # (or random), default is d
    > ms_sampling_to_pdbs.py a/path/to/mcce/output 7 0 99 -sampling_kind r --only_create_smsm  # do not write the pdbs
```
