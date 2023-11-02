I ran the following code against my 4lzt output folder
to verify that NITER * num_free * RUNS = MC.counts, but could not.

 * Note: jj_ms_analysis is the module from your demo without the hard coded paths.
---

```python
from pathlib import Path
import jj_ms_analysis as jjms
import subprocess


# may need to change this:

#DIR_DATA = Path.cwd()
mcce_dir = DIR_DATA.joinpath("4lzt")

fname = mcce_dir.joinpath("ms_out/pH6eH0ms.txt")

jjmc = jjms.MC(mcce_dir)
jjmc.readms(fname)
```
Output:
```
Reading MC:0
Reading MC:1
Reading MC:2
Reading MC:3
Reading MC:4
Reading MC:5
```

```python
runprm = mcce_dir.joinpath("run.prm.record")
runprm2 = runprm.relative_to(Path.cwd().parent)

prmdata = subprocess.check_output(f"grep -E 'MONTE_T)|MONTE_FLIPS|MONTE_RUNS|MONTE_NITER' ../{runprm2}| sed -e 's/(//g; s/MONTE_//g; s/)//g'",
                                  stderr=subprocess.STDOUT,
                                  shell=True
                                  ).decode().splitlines()

T = float(prmdata[0].split()[0])
FLIPS = int(prmdata[1].split()[0])
RUNS =  int(prmdata[2].split()[0])
NITER = int(prmdata[3].split()[0])

counts = 0
for mc in jjmc.microstates:
    counts += mc.count
print(f"Via iterative sum: {counts = :,}, {jjmc.counts = :,}")

print(f"{counts/6=:,.0f}")

jjmc_free = len(jjmc.free_residues)

# MS.counts should be MONTE_NITER x number of free conformers x MC runs.
jjmc_calc = NITER * jjmc_free * RUNS

print(f"NITER({NITER}) * jjmc_free({jjmc_free}) * RUNS({RUNS}) = {jjmc_calc = :,}")
```
Output:
```
Via iterative sum: counts = 3,300,000, jjmc.counts = 3,300,000
counts/6=550,000
NITER(2000) * jjmc_free(63) * RUNS(6) = jjmc_calc = 756,000
```
