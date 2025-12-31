# Using pixi in GERE

## Install pixi
Ironically, we have to install pixi through conda.

```bash
conda create -n pixi pixi
```

Activate the conda environment to access pixi

```bash
conda activate pixi
```

## Set pixi home directory
We want to access the same pixi configuration whether we run pixi from the GERE desktop or from the cluster.

To acheive this, the pixi configuration files need to be in a directory accessible to both filesystems. The `/re_gecip/` directories acheive this purpose. 

The `PIXI_HOME` environment variable can't be updated in the `.bashrc`:

```bash
# .bashrc
export PIXI_HOME="/re_gecip/enhanced_interpretation/AlexBlakes/.pixi" # Fails

```

Instead we can manually create this directory...

```bash
mkdir "/re_gecip/enhanced_interpretation/AlexBlakes/.pixi"
```

... then create softlinks from here to our `HOME` directories:

```bash
# Repeat this from both the GERE desktop and an interactive node
ln -s "/re_gecip/enhanced_interpretation/AlexBlakes/.pixi" ~/.pixi
```

## Global config
The GERE is air-gapped. We need to specify mirrors for conda and an index-url for pypi.

Add the following to the `config.toml` in the pixi "home" directory, created above:

```toml
# Pixi global config

default-channels = ["conda-conda-forge", "conda-bioconda"] # Order matters

[mirrors]
"https://conda.anaconda.org" = ["https://artifactory.aws.gel.ac:443/artifactory/"]

[pypi-config]
index-url = "https://artifactory.aws.gel.ac/artifactory/api/pypi/pypi/simple"
```

## Mapping conda and pypi packages
In order that pypi installs should not cause duplication / overwriting of conda packages, a mapping between the two is needed.