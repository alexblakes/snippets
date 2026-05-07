# Using pixi in GERE
This is [Pixi](https://github.com/prefix-dev/pixi/).
Here's their [docs](https://pixi.prefix.dev/latest/).
Using Pixi has greatly improved my GEL experience. This guide describes the installation and setup of Pixi in the GERE.

## Table of contents
TBA

## Create a "home" for Pixi
We want to access the same pixi configuration whether we run pixi from the GERE desktop or from the HPC.

To achieve this, the pixi configuration files need to be in a directory accessible to both filesystems. The `/re_gecip/` directories achieve this purpose.

Make a directory called `.pixi` in your personal `/re_gecip/` directory. For example:
```bash
mkdir -p /re_gecip/enhanced_interpretation/<my_gel_username>/.pixi # PIXI_HOME
```

Now, make a soft link from this directory to your `$HOME` directory:
```bash
# Repeat this from both the GERE desktop and an interactive node
ln -s "/re_gecip/enhanced_interpretation/<my_gel_username>/.pixi" ~/.pixi
```


## Install Pixi
Pixi is not installed by default in the GERE. Here's how you get it.

Copy the compressed archive to your `$PIXI_HOME` directory, then unpack it:
```bash
cp /re_gecip/shared_allGeCIPs/AlexBlakes/pixi-x86_64-linux.tar.gz $PIXI_HOME # v0.63
cd $PIXI_HOME
tar -xzf pixi-x86_64-linux.tar.gz
```

Pixi should now be installed! Try it from the current directory (`$PIXI_HOME`):
```bash
pixi -h
```

## Access Pixi from anywhere
The Pixi binary is now available in your `/re_gecip/` directory. Now we need to tell the shell how to find it, even when you're not in that directory.

Add the following lines to your `~/.bashrc` in both the Desktop and HPC:
```bash
# Add this to the end of ~/.bashrc
# Repeat this in both the GERE Desktop and an interactive session on the HPC

# Set PIXI_HOME environment variable
# >>> EDIT THIS PATH TO YOUR OWN RE_GECIP DIRECTORY <<<
export PIXI_HOME="/re_gecip/enhanced_interpretation/<my_gel_username>/.pixi"

# Add PIXI_HOME and pixi-installed binaries to PATH
PATH="${PIXI_HOME}:${PIXI_HOME}/bin:${PATH}"
export PATH
```

In both the Desktop and HPC, restart your shell:
```bash
bash # Restart the shell
```

Check that `$PIXI_HOME` has been set:
```bash
echo $PIXI_HOME # Should print the directory we made above.
```

Now try to run Pixi from your home directory:
```bash
cd 
pixi -h
```

Not working? Double check that pixi is in your `PATH`:
```bash
echo $PATH 

# Or 
which -a pixi
```

Missing from your `PATH`? Restart your shell:
```bash
bash # Restart the shell
```

Still missing? Double check your `.bashrc` for typos in the earlier configuration step. Don't forget to apply any changes to the `.bashrc` files in both the GERE Desktop and HPC, and to restart your shell with `bash`.


> [!TIP]
> Note that your GERE desktop and the HPC have different $HOME directories, each with its own `.bashrc`. Changes you make in one aren't reflected in the other.
> 
> One work-around is to keep a "central" `.bashrc` in your /re_gecip/ folder (e.g. `/re_gecip/enhanced_interpretation/<your_gel_username>`), and to have softlinks to this file (each named `.bashrc`) in the Desktop and HPC `$HOME` directories. That way, you only have to change the central file for the effects to be consistent in both contexts. We did the same thing for our `$PIXI_HOME` directory, above.

## Installing with conda
Alternatively, one can install pixi through conda... But it slightly defeats the point to have nested package managers...
```bash
conda create -n pixi pixi

conda activate pixi
```

## Global configuration (config.toml)
The GERE is air-gapped. We need to specify mirrors for conda and an index-url for pypi.

pixi searches for global config settings in `~/.pixi/config.toml` by default. Create this file:
```bash
touch "${PIXI_HOME}/config.toml"
```

Add the following lines to that file:
```toml
# Global pixi configuration in ${PIXI_HOME}/config.toml

default-channels = ["conda-conda-forge", "conda-bioconda"] # Order matters

[mirrors]
"https://conda.anaconda.org" = ["https://artifactory.aws.gel.ac:443/artifactory/"]

[pypi-config]
index-url = "https://artifactory.aws.gel.ac/artifactory/api/pypi/pypi/simple"
```

## Mapping conda and pypi packages
In order that pypi installs should not cause duplication / overwriting of conda packages, a mapping between the two is needed.

Because the GERE is airgapped, this mapping must be provided locally.

A recommended mapping is available [here](https://github.com/prefix-dev/parselmouth/blob/main/files/compressed_mapping.json).

Copy and paste this `compressed_mapping.json` into the pixi "home" directory, described above.

```json
{"21cmfast": "21cmfast", "2dfatmic": null, "4ti2": null, ...}
```

> [!TIP]
> You can instead provide an empty JSON file (`{}`). If this empty JSON is provided, Pixi will install packages from both PyPi and Conda, rather than from Conda alone.

## Global installs with pixi
If you are just interested in having newer versions of tools available globally, try Pixi global installs. For example:
```bash
pixi global install bcftools
bcftools --version
```

Here are a few tools that may improve your GERE experience:
```bash
pixi global add ripgrep bat ensembl-vep bcftools bedtools starship parallel
```

Never `module load bcftools` again!

## Per-project usage of Pixi
More idiomatic and reproducible use of Pixi involves setting up a virtual environment for each project / package.

A little more setup is required.

Configuration for each pixi project is set in `pixi.toml` in the root of that project.

Start a new project with
```bash 
pixi init <my_project>
``` 

A new `pixi.toml` file is automatically added to the `<my_project>` directory. Add the following lines to it:

```toml
[workspace] 
# Append these lines to the `workspace` table
# The conda-pypi mapping is provided for each channel.
conda-pypi-map = { conda-conda-forge = "/re_gecip/enhanced_interpretation/<my_gel_username>/.pixi/compressed_mapping.json", conda-bioconda = "/re_gecip/enhanced_interpretation/<my_gel_username>/.pixi/compressed_mapping.json" }

[system-requirements] # See https://pixi.prefix.dev/dev/workspace/system_requirements/
# Obtained with `ldd --version`
libc = { family = "glibc", version = "2.26" }
```

> [!WARNING]
> Unfortunately, these settings cannot be configured in the global config at present. (See GitHub issues [#1795](https://github.com/prefix-dev/pixi/issues/1795) and [#3199](https://github.com/prefix-dev/pixi/issues/3199)). So this needs to be repeated for each project.

Now try adding Python to your list of package requirements:
```bash
pixi add python=3.13
pixi run python --version
```

PyPi packages are installed with the `--pypi` flag:
```bash
pixi add --pypi statannotations
```


## Running pixi in VS Code
The GERE uses an old version of VS Code, with an out-of-date Python extension which doesn't recognise pixi virtual environments.

To work around this, open a terminal in the GERE desktop, navigate to your pixi workspace, and run:

```bash
pixi run code .
```

Ta-da!
