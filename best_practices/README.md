# Best Practices for Software Development
Note: Most of these __strong__ recommendations will be implemented for linux-like paltforms (Linux, Apple, WSL2), and assume that conda is the package manager in use.


[WIP}: The links in the the following contents will be active when the underlaying documents are created: this is a road map.

# 1. Customize your `.bashrc` file
It make navigation easier and safer (among many other things).
 * <a href="https://github.com/CatChenal/mcce_dev/blob/main/best_practices/Customize_bashrc.md" target="_blank">Customize_bashrc (md)</a>
 * <a href="https://github.com/CatChenal/mcce_dev/blob/main/best_practices/Customize_bashrc.ipynb" target="_blank">Customize_bashrc (nb)</a>

# 2. Use the `python-dotenv` library
Keep your local environemnt and credential private... And use `.gitignore` to make sure they really are.
 * [Use_python-dotenv]()

# 3. Organize your project folder like a pro
There are many tools relying on the [cookiecutter](https://cookiecutter.readthedocs.io/en/latest/installation.html#install-cookiecutter) approach, notably [cookicutter-cms](https://github.com/MolSSI/cookiecutter-cms)
> Python-centric Cookiecutter for Molecular Computational Chemistry Packages
 * [Project_organization]()

# 4. Tips and tricks when working with notebooks
Recommended: [jupyterlab](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html) and [nb_conda_kernels](https://github.com/Anaconda-Platform/nb_conda_kernels) installed in the __base__ conda environment, and [ipkernel](https://ipython.readthedocs.io/en/5.x/install/kernel_install.html) installed in each dedicated `conda/envs`.
 * [Notebooks_tips]()

# 5. SOLID Principles
Best practices for OOP development with a solid acronym. Courtesy of [dmmeteo](https://github.com/dmmeteo)'s gist
 * <a href="https://github.com/CatChenal/mcce_dev/blob/main/best_practices/SOLID.md" target="_blank">SOLID Principles (md)</a>
 * <a href="https://github.com/CatChenal/mcce_dev/blob/main/best_practices/SOLID.ipynb" target="_blank">SOLID Principles (nb)</a>
