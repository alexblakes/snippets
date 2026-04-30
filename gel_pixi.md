# Using pixi in GERE
This is [Pixi](https://github.com/prefix-dev/pixi/).
Here's their [GitHub](https://pixi.prefix.dev/latest/).
Using Pixi has greatly improved my GEL experience.

## Set pixi home directory
We want to access the same pixi configuration whether we run pixi from the GERE desktop or from the HPC.

To acheive this, the pixi configuration files need to be in a directory accessible to both filesystems. The `/re_gecip/` directories acheive this purpose.

For example:
```bash
mkdir -p /re_gecip/enhanced_interpretation/AlexBlakes/.pixi
```

The `PIXI_HOME` environment variable can be set in your `~/.bash_profile`, or `~/.bashrc`. For example:

```bash
# ~/.bashrc
export PIXI_HOME="/re_gecip/enhanced_interpretation/AlexBlakes/.pixi"
```

Run `bash` to restart your shell. Check that `$PIXI_HOME` has been set with `echo $PIXI_HOME`

**Optionally**, we can then create softlinks from here to our `$HOME` directories:

```bash
# Repeat this from both the GERE desktop and an interactive node
ln -s "/re_gecip/enhanced_interpretation/AlexBlakes/.pixi" ~/.pixi
```

pixi searches for global config settings in `~/.pixi/config.toml` by default.

## Install pixi
pixi is not currently available in the GERE. 

Copy the compressed archive to your home directory, then unpack it.

```bash
cp /re_gecip/machine_learning/CRDG_Shared/tools/pixi-x86_64-linux.tar.gz $PIXI_HOME
cd $PIXI_HOME
tar -xzf pixi-x86_64-linux.tar.gz
```

Alternatively, we can install pixi through conda.

```bash
conda create -n pixi pixi

conda activate pixi
```

Now add pixi to your path

```bash
#.bashrc
PATH="${PIXI_HOME}/pixi:${PATH}"
export PATH
```

## Add binaries in $PIXI_HOME to your $PATH
```
# ~/.bashrc
PATH="${PATH}:${PIXI_HOME}/bin"
export PATH
```
Note that your GERE desktop and the HPC have different $HOME directories, each with it's own `.bashrc`. Changes you make in one aren't reflected in the other.

You should make these changes separately in each `.bashrc` file.

Alternatively, I keep a "central" `.bashrc` in my /re_gecip/ folder (`/re_gecip/enhanced_interpretation/AlexBlakes`). I have softlinks to this file in my Desktop and HPC `$HOME` directories. That way, I only have to change the central file, for the effects to be consistent in both contexts. We did the same thing for our `$PIXI_HOME` directory, above.

## Global configuration (config.toml)
The GERE is air-gapped. We need to specify mirrors for conda and an index-url for pypi.

Add the following to the `config.toml` in the `PIXI_HOME` directory, described above:

```toml
# ${PIXI_HOME}/config.toml

# Global pixi configuration

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
# compressed_mapping.json
{"21cmfast": "21cmfast", "2dfatmic": null, "4ti2": null, ...}
```

This mapping is referenced in the pixi manifest for each pixi directory, described below.

## Global installs with pixi
If you are just interested in having newer versions of tools available globally, try Pixi global installs. For example:
```bash
pixi global install bcftools
bcftools --version
```

Never `module load bcftools` again!

## Customising the pixi manifest (pixi.toml)
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
conda-pypi-map = { conda-conda-forge = "/re_gecip/enhanced_interpretation/AlexBlakes/.pixi/compressed_mapping.json", conda-bioconda = "/re_gecip/enhanced_interpretation/AlexBlakes/.pixi/compressed_mapping.json" }

[system-requirements] # See https://pixi.prefix.dev/dev/workspace/system_requirements/
# Obtained with `ldd --version`
libc = { family = "glibc", version = "2.26" }
```

Unfortuntaly, these settings cannot be configured in the global config at present. (See GitHub issues [#1795](https://github.com/prefix-dev/pixi/issues/1795) and [#3199](https://github.com/prefix-dev/pixi/issues/3199)). So this needs to be repeated for each project.

## Running VS Code
The GERE uses an old version of VS Code, with an out-of-date Python extension which doesn't recognise pixi virtual environments.

To work around this, open a terminal in the GERE desktop, navigate to your pixi workspace, activate the pixi conda environment, and run:

```bash
pixi run code .
```

Ta-da!
