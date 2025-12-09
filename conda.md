# Conda snippets

## Data science conda .yml
```bash
name: "<some_name>"
channels:
  - conda-forge # conda-conda-forge (GEL)
  - bioconda # conda-bioconda (GEL)
  - defaults # conda-main (GEL)
dependencies:
  - black
  - gtfparse
  - ipykernel
  - ipython
  - jupyterlab
  - matplotlib
  - nb_conda_kernels
  - numpy
  - pandas
  - pip
  - python=3.13
  - ruff
  - scikit-learn
  - scipy
  - seaborn
  - setuptools
  - statannotations
  - pip:
    - adjusttext
    - pandas-checks
```
