# mcce_dev

Project for ancillary development of [Stable-MCCE repo](https://github.com/GunnerLab/Stable-MCCE).

# Current projects:

## ms_sampling_to_pdbs
### -> Provide an executable module with cli to obtain a collection of pdb files created from sampled microstates.
### Platforms: Linux and Linux-like

#### The following documentation pertains to using the module as a command line interface:

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

    The remaining arguments are optional and will perform the action they describe only if included:
        --reverse_sort
        --only_create_smsm
        --overwrite_split_files   :: use only once: remove the flag on subsequent runs
        --clear_pdbs_folder
        --list_files
    i.e:
    > ms_sampling_to_pdbs.py a/path/to/mcce/output 7 0 99 --only_create_smsm  # create the smsm file, but not the pdbs
```

# Ad hoc installation to always use the latest update of `ms_sampling_to_pdbs` (which also enables its use programmatically):

  0. Fork [mcce_dev](https://github.com/CatChenal/mcce_dev/tree/main)
  1. Clone the fork locally.
  2. `git pull` any update.
  3. At the command line, run this command (You should get a path ending with /bin/mcce):
     ```
     > which mcce
     ```
     Alternatively, you can use any other bin folder &mdash; e.g. usr/local/bin &mdash; as the mcce excecutable is not required for this project.
  4. cd to that bin directory
  5. Soft-link the ms_sampling_to_pdbs.py module of your clone into it; Run this command:
     ```
      > ln -s <path to your clone>/src/ms_sampling_to_pdbs.py .
     ```
     
From then on, if I push an update, you only have to
repeat step #2.

# Tutorial

Once your are setup as described above, you can go through the tutorial to learn more about the module, but more importantly about
the microstates in any MCCE output folder you have, granted you have run Step 4 with the --ms_out flag.

* [__Tutorial notebook__](https://github.com/CatChenal/mcce_dev/blob/main/notebooks/tutorial.ipynb)
* [__Tutorial markdown__](https://github.com/CatChenal/mcce_dev/blob/main/notebooks/tutorial.md)
   - NOTE: The `tutorial.md` file is a `jupytext` Markdown file paired to the notebook: you can use it instead of the notebook in case you cannot or don't want
to launch jupyter, however you need to disregard its special header section which is used for synchronizing.

