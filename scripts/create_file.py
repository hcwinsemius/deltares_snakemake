# This first line is used to save the python code below into a file called "src/create_file.py"

######################################################################################
####################### Content of the python script #################################
from pathlib import Path

# Main function
def create_file(ouputfile: Path, lines: list) -> None:
    """
    Create a file with each line contant stored as items in ``lines`` list. Saves it to ``outputfile``.
    """
    with open(outputfile, "w") as file:
        for line in lines:
            file.write(f"{line}\n")


# Get the snakemake arguments
outputfile = Path(str(snakemake.output))
line1 = snakemake.params.line1
line2 = snakemake.params.line2

# Run our function
create_file(outputfile, lines=[line1, line2])
