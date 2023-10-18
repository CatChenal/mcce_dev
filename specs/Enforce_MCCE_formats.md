# Enforce MCCE formats


## Ticker tape format (files in ms_out folder)
If you peruse the code in bin/ms_analysis.py, you'll notice that there is a lot of processing done to
exclude comments and blank lines.
Currently, the "header" looks like this:
```
T:298.15,pH:5.00,eH:0.00
METHOD:MONTERUNS
#N_FIXED:FIXED_CONF_ID
65:0 2 4 5 11 17 18 19 20 [...] 263 264 270 275
#N_FREE residues:CONF_IDs for each free residues
63:6 7 8 9 10 ;12 13 14 15 ;24 25 26 27 ;[...] ;301 302 303 304 30
#EVERY MONTERUN START FROM A NEW STATE
#MC:ITER_MONTERUNS
#N_FREE: FREE_CONF_ID,
#ENERGY, COUNT,NEW_CONF

MC:0
```

I propose that it look like this instead:
```
T:298.15,pH:5.00,eH:0.00
METHOD:MONTERUNS
65:0 2 4 5 11 17 18 19 20 [...] 263 264 270 275
63:6 7 8 9 10 ;12 13 14 15 ;24 25 26 27 ;[...] ;301 302 303 304 30
MC:0
```

## Changes in C codebase:
Note: This would require:
  - a bump in version number
  - adding simplified readers of ms_out/files
  - create a switch to use either the simplyfied or old reader function depending on MCCE version number.


### Add `ms_out/format.info`
All those comments/format information should be place an ad-hoc file, e.g. `ms_out/format.info`
created whenever the ms_out folder is created (or retained), so that the ticker tape file would be __pure
data__ without comments or blank lines. (`format.info` is technically called a "metadata" file.)
The content of that file could be:
```
Text files in this folder (ms_out/) store information about microstates created by MCCE Monte Carlo runs
in a compact format. These files are also referred to as "ticker tapes" or "ms_out files".
The files are named using two important experimental variables: pH and Eh, e.g.: "pHXeHYms.txt", where X
and Y are the values of pH, Eh for each titration point defined in run.prm.

The format of these files is the following:

Line 1: Experimental variables for temperature (K), pH and eH; MUST start with "T".
Line 2: The method used to generate the microstates; MUST start with "METHODS".
Line 3: Fixed conformers line: N_FIXED (integer):FIXED_CONF_IDX (space-separated list).
Line 4: Free residues line: N_FREE (integer):FREE_CONF_IDX (grouped & space-separated list). The conformer
        indices are grouped by residue using " ;" as a delimiter.
Line 5: Marker for the first MC run in the file; runs are zero-indexed, currently up to 5. MUST start with "MC:"
Line 6: Initial state (sampled from FREE_CONF_IDX); same line format as Line 3.
Line 7: First comma-separated record out of ITER_MONTERUNS: ENERGY, COUNT, NEW_CONF, i.e: energy of the
        microstate, number of times the state was chosen, which conformer(s) was(were) flipped, if any.
        Every MONTERUN starts with a new state.
Lines [8 to ITER_MONTERUNS]: repeat of Line 7.
Line [ITER_MONTERUNS + 1]: repeat of Line 6
[...]

Note: x_CONF_IDX referred above is the index from head3.lst but starting from 0 (zero-based)..
```

## Changes in python codebase:
  - add simplified reader function for strict format
  - maintaining the current python "read ms_out" code so that ms_out files from older versions can still be processed.
  - create a switch to use either the simplyfied or old reader function depending on MCCE version number.


### Raise error for format violation
Enforcing this format would mean raising errors if the first line does not start with "T", the second
line does not start with 'METHOD', etc., with the error message referencing `format.info`.

When parsing a ms_out file in python, the error message could be the following:
```
# Say `msout_file` store the name of the file being processed:

ERR_MSG_MS_FRMT_pf = """This ms_out file: {} does not comply with the expected format.
No ticker tape file in ms_out should be modified. If you wish to add annotations, comments,
explanations, etc. about the ms_out folder or a specific fle, please create a separate file.
The file ms_out/format.info contains the file format specifications.
"""

[...]

raise ValueError(ERR_MSG_MS_FRMT_pf.format(msout_file)

# Implementation note: the "_pf" ending in message variables indicate it is prepped for the
# print `.format()` function: it contains the empty variable placeholder(s) "{}".

```
