# deltares_snakemake: Snakemake training Deltares
This repository contains material to discover Snakemake and its functionalities as a smart and reproducible workflow manager.

It uses Jupyter notebooks and requires a Python installation.

Two exercises are prepared around two workflows:

- a simple workflow for file manipulation in Windows: [notebooks/1_first_workflow.ipynb](https://github.com/hcwinsemius/deltares_snakemake/blob/main/notebooks/1_first_workflow.ipynb)
- a more advanced workflow with a simple weather generator and hydrological model: [notebooks/2_rainfall_runoff_workflow.ipynb](https://github.com/hcwinsemius/deltares_snakemake/blob/main/notebooks/2_rainfall_runoff_workflow.ipynb)

## Installation
### Python and conda/mamba
You'll need Python 3.8 or greater and a package manager such as conda or mamba. These package managers help you to install (Python) packages and 
[manage environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) such that different installations do not conflict.

We recommend using the [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) Python distribution. This installs Python and the 
[mamba package manager](https://github.com/mamba-org/mamba). [Miniforge](https://github.com/conda-forge/miniforge) and 
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) will install Python and the [conda package manager](https://docs.conda.io/en/latest/).

### Download the content of deltares_snakemake repository
To run the exercices, you will need to download the content of the deltares_snakemake repository locally. You can either do a [manual download](https://github.com/hcwinsemius/deltares_snakemake/archive/refs/heads/main.zip) and extract the content of the dowloaded ZIP folder or clone the repository locally:

``` bash
git clone https://github.com/hcwinsemius/deltares_snakemake.git
```

### Install the python dependencies in a separate Python environment
The last step is to install all the python dependencies required to run the notebooks, including of course snakemake. All required dependencies can be found
in the [environment.yml](https://github.com/hcwinsemius/deltares_snakemake/blob/main/environment.yml) file. 

First navigate into the deltares_snakemake folder (where the environment.yml file is located). Create a new snakeymakey environment using the environment.yml file 
(you can exchange mamba/conda in the exmaple below):

``` bash
cd deltares_snakemake
mamba env create -f environment.yml
```

## Running the exercises
We have prepared two iPython notebook examples to discover snakemake in the notebooks folder:

- 1_first_workflow.ipynb: a simple workflow for file manipulation in Windows to discover how to write and run a Snakefile, determine rule dependency and use a configuration file.
- 2_rainfall_runoff_workflow.ipynb: a more advanced workflow with a simple weather generator and hydrological model to discover generalization and parallelization in Snakemake with wildcards and Paramspace.

You can open and run the notebooks in Visual Studio Code or as a Jupyter notebook:

``` bash
conda activate snakeymakey
cd notebooks
jupyter notebook
```
