# This first line is used to save the python code below into a file called "src/create_file.py"

######################################################################################
####################### Content of the python script #################################
import os
import pandas as pd
import matplotlib.pyplot as plt

# Main function
def plot_timeseries(
    csv_files: list, variables: list, outputfile: str, **kwargs
) -> None:
    """
    Reads different csv files in csv_files list and prepares a plot per variable in variables list.
    The plot will be saved to outputfile.

    kwargs can be used to pass specific arguments into the pandas.read_csv function.
    """
    # Initialiase plots
    print(
        f"Initialising plots for variables {variables} based on csv files: {csv_files}"
    )
    n = len(variables)
    fig, axes = plt.subplots(n, 1, sharex=True, figsize=(15, n * 4))
    axes = [axes] if n == 1 else axes

    # Loop over csv files
    for fn in csv_files:
        # Read file
        print(f"Reading csv file {fn}")
        df = pd.read_csv(fn, **kwargs)

        # Loop over varibales to add to plot if they are present in the csv file
        for i in range(len(variables)):
            var = variables[i]
            if var in df.columns:
                df[var].plot.line(
                    ax=axes[i],
                    x="time",
                    label=os.path.basename(fn),
                )
    # Set axis label/legend
    for i in range(len(variables)):
        axes[i].set_ylabel(variables[i])
        axes[i].legend()

    # Save plot
    print(f"Saving plot in {outputfile}")
    fig.savefig(outputfile, bbox_inches="tight", dpi=200)


# Get the snakemake arguments
csv_files = [
    "data/precip_evapo_weagen1.csv",
    "output/river_flow_weagen1.csv",
    "data/precip_evapo_weagen2.csv",
    "output/river_flow_weagen2.csv",
]
variables = ["p", "ep", "Q"]
outputfile = "output/plot_results_example.png"

# Run our function
plot_timeseries(
    csv_files=csv_files,
    variables=variables,
    outputfile=outputfile,
)
