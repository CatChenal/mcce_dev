# MCCE TODOs
---
---

# [Apply SOLID principles]()
* Put all (if not, most) io functions into one module
 - Example: mccesteps.py ->io_runprm.py
   - Rename:
    * def export_runprm(runprm) -> def runprm_dict_tofile(runprm) 
  - 

* Put all defaults and constants into one module


---
# GITHUB

# Initial

 * Rename branch from 'master' to 'main'
 * Create a 'dev' branch instead of having a different repo
 * Setup branch protection rules if not present
  - https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-protected-branches/managing-a-branch-protection-rule
  - https://johnnymetz.com/posts/disable-direct-push-to-main-branch/
 * Use semantic versioning system: X.Y.Z
  - Add __version__ in __init__.py?

 * Setup Changelog w/ PRs
 * Create a project page to facilitate:
  - Roadmap (& discussion)
  - List all TODOs related to coding, docs, packaging
  - Task assignment to members
  - Progress
 * Document installation process: 
  - Add specific link to Installation in README (https://gunnerlab.github.io/Stable-MCCE/quick/)
  - In "Quick Install": For python development, install packages inside an activated conda envir.
  - Under "Python and modules", module `requests` is missing from the list
    => all python requirements should be in a requiremnts.txt file.
 * Python3 version: which is in use?
    - suggested minimum: 3.8: https://devguide.python.org/versions/
 * Add linters for C/C++ and python:
  - https://github.com/caramelomartins/awesome-linters#cc
  - https://github.com/caramelomartins/awesome-linters#python
 * Add unit tests for python and C
  - E.g. running lysozyme through MCCE should be packaged as a standard test so it can be run after installation with make: `make lysozyme_test`.
    - https://codeahoy.com/learn/cprogramming/ch46/
 * conda channel name: change 'newbooks' (?) to 'gunnerlab'
  - https://docs.conda.io/projects/conda/en/22.9.x/user-guide/tasks/create-custom-channels.html

---

---
# Misc Resources:
 * https://edcarp.github.io/introduction-to-conda-for-data-scientists/04-sharing-environments/index.html

