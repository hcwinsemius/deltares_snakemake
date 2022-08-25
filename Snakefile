import shutil
import os
from snakemake.utils import Paramspace
import pandas as pd

# Path to the Snakemake config file (options for basins to build, data catalog, model output directory)
configfile: "config/snake_config.yml"
# when configfile is set, the file will be read and stored in variable "config"

# you can also define a parameter space for things that should be varied and that are easier stored as a csv file
paramspace = Paramspace(pd.read_csv("config/parameters.csv", sep=","), filename_params="*")
print(paramspace.wildcard_pattern)

# The master rule defines what eventually will be expected as output, basically it is a rule with only inputs
rule all:
    input:
    	# input_csv = config["input_csv"]
        output_csv = config["output_dir"] + "river_flow.csv"
    	

# rule to make
rule weather_generator:
    # there is not input section because this rule does not require any inputs
    output:
        input_csv = config["input_csv"]
    # parameters can be anything like strings, numbers, etc.
    params:
    	start_date = config["start_date"],
    	end_date = config["end_date"],
    	multiplier = config["multiplier"],
    	evapo = config["evapo"],
    	dry_day_fraction = config["dry_day_fraction"]
    script:
        "scripts/weather_generator.py"


# Rule to run the calibrated wflow model
rule run_simple_model:
    input:
        input_csv = config["input_csv"]

    output:
        output_csv = config["output_dir"] + "river_flow.csv"
    params:
        param_names = paramspace.wildcard_pattern,
        param_set = paramspace.instance
    script:
        "scripts/simple_hydrology.py"

rule plot_output:
	# TODO
