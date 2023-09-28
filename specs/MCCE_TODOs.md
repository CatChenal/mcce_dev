# MCCE TODOs
---
---
# General coding guidelines:

## [Apply SOLID principles](https://gist.github.com/dmmeteo/f630fa04c7a79d3c132b9e9e5d037bfd)
  * S: Single Responsibility Principle
  * O: Open-Closed Principle
  * L: Liskov Substitution Principle
  * I: Interface Segregation Principle
  * D: Dependency Inversion Principle

  + DRY: Don't Repeat Yourself

## Use a good IDE, e.g.:
 * VS Code
 * pycharm
 * JupyterLab

## TEST!
  * Implement tests at the function or module level
  * Always include `pytest` in your environment for testing python code
  * In C, [setup tests to run with `make`](https://codeahoy.com/learn/cprogramming/ch46/)

## Use code QA tools:
 ### In python:
   * [Black, code formater](https://github.com/psf/black)
    > Blackened code looks the same regardless of the project you're reading. Formatting becomes transparent after a while and you can focus on the content instead. Black makes code review faster by producing the smallest diffs possible.

   * [Pylint](https://github.com/pylint-dev/pylint)
     > Pylint analyses your code without actually running it. It checks for errors, enforces a coding standard, looks for code smells, and can make suggestions about how the code could be refactored.

   * [Flake8: Your Tool For Style Guide Enforcement](Flake8: Your Tool For Style Guide Enforcement)
 
 ### In C/C++:
   * [CPPCheck](https://sourceforge.net/projects/cppcheck/reviews/)

 ### In notebooks:
   * [NbQA](https://nbqa.readthedocs.io/en/latest/index.html)
   * Also use:
     - [JupyterLab](https://github.com/mwouts/jupytext)
     - [Jupytext](https://github.com/mwouts/jupytext)


## Documentation:
 * Use "docstrings" for modules and functions (see style to follow in link below).
 * Use an automated code documenter such as [Doxygen](https://www.doxygen.nl/manual/docblocks.html)
  >> The wiki could link to it under "API".
 * [Python comment block example](https://www.doxygen.nl/manual/docblocks.html#pythonblocks)
 * [C++ code block example](https://www.doxygen.nl/manual/docblocks.html)

## Use Git repos and Git projects
 * Good way to manage your (revertible) changes
 * Create a [project](https://docs.github.com/en/issues/planning-and-tracking-with-projects/learning-about-projects/about-projects) for each repo you have: split tasks into small chuncks, and assigned them to yourself or a collaborator. Smile when you see your accomplishments!


---
---
# GITHUB - MCCE code base

# Initial

 * Rename branch from 'master' to 'main'
 * Create a 'dev' branch instead of having a different repo
 * Setup branch protection rules if not present
  - https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-protected-branches/managing-a-branch-protection-rule
  - https://johnnymetz.com/posts/disable-direct-push-to-main-branch/
 * Use semantic versioning system: X.Y.Z
  - Add __version__ in __init__.py?

## Setup Changelog w/ PRs
## Create a project page to facilitate:
  - Roadmap (& discussion)
  - List all TODOs related to coding, docs, packaging
  - Task assignment to members
  - Progress
 * => This is where each item in this document should be

 * Document installation instructions: 
  - Add specific link to Installation in README (https://gunnerlab.github.io/Stable-MCCE/quick/)
  - In "Quick Install": For python development, install packages inside an activated conda envir.
  - Under "Python and modules", module `requests` is missing from the list
    => all python requirements should be in a requirements.txt file.

 * Python3 version: which is in use?
    - suggested minimum: 3.8: https://devguide.python.org/versions/

 * Add linters for C/C++ and python:
  - https://github.com/caramelomartins/awesome-linters#cc
  - https://github.com/caramelomartins/awesome-linters#python


## Put all (if not, most) io functions into one module
  * Example: mccesteps.py & all other read/write functions: -> mcce_io.py

## Put all defaults and constants into one module, e.g.: mcce_config.py

## Use type hints in python code

## Rename:
  * def export_runprm(runprm) -> def dict_to_file(d: dict, which='runprm')


 * Add unit tests for python and C
  - E.g. running lysozyme through MCCE should be packaged as a standard test so it can be run after installation with make: `make lysozyme_test`.
    - https://codeahoy.com/learn/cprogramming/ch46/
 * conda channel name: change 'newbooks' (?) to 'gunnerlab'
  - https://docs.conda.io/projects/conda/en/22.9.x/user-guide/tasks/create-custom-channels.html

---

---
# Misc Resources:
 * [Conda for scientists](https://edcarp.github.io/introduction-to-conda-for-data-scientists/04-sharing-environments/index.html)
 * [Notebooks tips & tricks](https://www.dataquest.io/blog/jupyter-notebook-tips-tricks-shortcuts/)

