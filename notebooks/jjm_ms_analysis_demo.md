# Microstate Analysis Library Reference

## Microstate Analysis Library Reference

### Import library

Suppose the ms\_analysis.py is in the current working directory or the Python site-packages directory.

```python
from ms_analysis import *
```

### Global constants

Once the library is loaded, two global constants (at temperature 298.15 K) are available:

**ph2Kcal**: Convert ph unit to Kcal/mol  
**Kcal2kT**: Convert Kcal/mol to kT

### Load a microstate file

Go to a working directory. The essential files for microstate analysis are:

- head3.lst file
- ms\_out folder that contains Monte Carlo sampling microstate output

You need to specify which file to load, such as *ms\_out/pH5eH0ms.txt. The name indicates the pH and Eh condition.*

A monte carlo object is reqired to be initialized to hold the microstates with `MC()`

Finally, read the data into the object with `readms()` method.

Example:

```python
cd ~/ms_analysis/4lzt
msfile = "ms_out/pH5eH0ms.txt"
mc = MC()
mc.readms(msfile)
```

Load partial Monte Carlo results. A Monte Carlo sampling is carried out 6 times and is numbered as 0, 1, 2, ..., 5. One can choose to load some of them:

```python
mc.readms(msfile, MC=[1,2])
```


### Data structure:

#### Conformer:

Conformer is a class object.

##### Variables:

- **iconf**: Integer - index of conformer, starting from 0
- **confid**: String - conformer name as in head3.lst
- **resid**: String - unique residue name including name, chain ID and sequence number
- **crg**: Float - net charge


#### Microstate:

Microstate is a class object.

##### Variables:

- **stateid**: String - compressed and encoded string to identify a microstate
- **E**: Float - microstate energy
- **count**: Integer - how many times this microstate is accepted

##### Function:

- **state()**: return a microstate, which is a list of selected conformers

#### Charge\_Microstate:

Charge\_Microstate is a class object. If we only care about residue ionization, we can reduce conformer microstates to charge microstates.

##### Variables:

- **crg\_stateid**: String - compressed and encoded string to identify a charge microstate
- **average\_E**: Float - average charge microstate energy
- **count**: Integer - how many times this charge microstate is accepted

##### Function:

- **state()**: return a charge microstate, which is a list of net charges, in the same order of free residues

#### Subset\_Microstate:

Subset\_Microstate is a class object. If we only care about a selected group of residues, we can group microstates based on the conformer selection of these residues only.

##### Variables:

- **subset\_stateid**: String - compressed and encoded string to identify a subset microstate
- **average\_E**: Float - average subset microstate energy
- **count**: Integer - how many times this subset microstate is accepted

##### Function:

- **state()**: return a subset microstate, which is a list of selected conformers of interested residues

#### Free\_res:

Free\_ress is a class object. It holds information of a free residue.

##### Variables:

- **resid**: String - residue identification name
- **charges**: list of floating point numbers - a list of charge choices


#### MC:

MC is a class object. It holds information of a Monte Carlo microstates output.

##### Variables:

- **T**: Float - Monte Carlo sampling temperature
- **pH**: Float - Monte Carlo sampling pH
- **Eh**: Float - Monte Carlo sampling Eh
- **method**: String - This indicates the microstates output is from either Monte Carlo sampling or Analytical Solution
- **counts**: Integer - Total number of Monte Carlo steps
- **conformers**: A list of Conformer objects that matche the entries in head3.lst
- **iconf\_by\_confname**: A dictionary that returns conformer index number from conformer name
- **fixedconfs**: A list of fixed confomer index numbers
- **free\_residues**: A list of conformer groups (each group is a list of conformer indicies) that make up free residues
- **free\_residue\_names**: A list of free residue names
- **microstates**: A list of Microstate objects. They are accepted microstates.

##### Function:

- **readms(fname, MC=\[\])**: read microstate output file and return a list of microstates. You can optionally choose what parts of Monte Carlo output to load. MC=\[\] means to choose all. MC=\[1,2\] means to choose 1st and 2nd MC runs. The valid numbers are from 0 to 5.
- **get\_occ(microstates)**: Convert a list of microstates to occupancy. It reads in a list of conformers and returns a list of occupancy (0.0 to 1.0) numbers on each conformer. <p class="callout warning">This function does not work on charge microstates or subset microstates.</p>
- **confnames\_by\_iconfs(iconfs)**: Convert a list if conformer indices to a list of conformer names.
- **select\_by\_conformer(microstates, conformer\_in=\[\])**: Select from given microstates if confomer is in the list. Return all if the list is empty. The input conformer\_in is a list of conformer names.
- **select\_by\_energy(microstates, energy\_in=\[\])**: Select from given microstates if the microstates' energy is within the range defined by energy\_in. energy\_in should be given an array with lower bound (inclusive) and a higher bound (exclusive).
- **convert\_to\_charge\_ms()**: Convert all microstates to a list charge microstate objects.
- **convert\_to\_subset\_ms(res\_of\_interest)**: Convert all microstates to a list subset microstate objects. The input res\_of\_interest is a list of residues of interest, in the form of residue names. These residues have to be free residues.

### Functions:

#### get\_erange(microstates)

Get the energy range of given microstates.

##### Input:

- **microstates**: A list of microstates object

##### Output:

A list of two numbers that are lower bound and higher bound of energey

#### bin\_mscounts\_total(microstates, nbins=100, erange=\[\])

Divide microstates into bins based on energy and get the counts of total steps in each bin.

##### Input:

- **microstates**: A list of microstates object
- **nbins**: the number of desired bins. Default value is 100
- **erange**: custom energy range. It is a list of lower bounds of bins

##### Output:

It returns two lists. The first list is the energy range in the form of lower bounds. The second list is number of microstate counts of each bin.

#### bin\_mscounts\_unique(microstates, nbins=100, erange=\[\])

Divide microstates into bins based on energy and get the counts of unique microstates in each bin.

##### Input:

- **microstates**: A list of microstates object
- **nbins**: the number of desired bins. Default value is 100
- **erange**: custom energy range. It is a list of lower bounds of bins

##### Output:

It returns two lists. The first list is the energy range in the form of lower bounds. The second list is number of microstate counts of each bin.

#### get\_count(microstates)

Divide microstates into bins based on energy and get the counts of unique microstates in each bin.

##### Input:

- **microstates**: A list of microstates object
- **nbins**: the number of desired bins. Default value is 100
- **erange**: custom energy range. It is a list of lower bounds of bins

##### Output:

It returns two lists. The first list is the energy range in the form of lower bounds. The second list is number of microstate counts of each bin.

#### average\_e(microstates)

Calculate the average energy of given microstates.

##### Input:

- **microstates**: A list of microstates object

##### Output:

Average energy.

### Code and example:

- Library: [ms\_analysis.py](https://mccewiki.levich.net/attachments/2)
- Demo: [demo.ipynb](https://mccewiki.levich.net/attachments/1)