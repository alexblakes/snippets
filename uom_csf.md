# UoM CSF

## Installing Python packages to conda environments with pip.

1. Activate a new (or existing) conda environment.
2. Ensure that pip is installed into that environment (i.e. `conda install pip`).
3. Double check that the conda environment's pip executable is used (`which -a pip`).

To be pedantic about this you could run:
```bash
"${CONDA_PREFIX}/bin/pip" install <package>
```

To be super pedantic you could run:
```bash
"${CONDA_PREFIX}/bin/python" -m pip install <package>
```
4. Run `pip install <package>` (or one of the above commands)
5. Check the install with `conda list <package>`

Sources:
- [Anaconda best practices for using pip with conda](https://www.anaconda.com/blog/using-pip-in-a-conda-environment)
- [Bits of Analytics blog post on exactly this issue](https://bitsofanalytics.org/posts/pip-conda-local-dev/pip_conda_local_dev.html)

## Install local packages to conda environments with pip
1. Read the instructions above
2. Navigate to the package directory
3. Perform an editable install with `pip install -e .`

See these resources for more information on packaging python projects:
- [The python packaging user guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
- [A simple example `pyproject.toml`](https://github.com/space-physics/lowtran/blob/main/pyproject.toml)
- [This useful blog post and references on pyproject.toml files](https://www.scivision.dev/python-minimal-package/)
- [Local project installs (pip docs)](https://pip.pypa.io/en/stable/topics/local-project-installs/)
- [Editable installs (setuptools docs)](https://setuptools.pypa.io/en/latest/userguide/development_mode.html)
