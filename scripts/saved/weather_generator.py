import os
import hydromt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# SnakeMake arguments are automatically passed WITHOUT importing a lib. So this .py file is basicallly used as a template.
start_date = snakemake.params.start_date
end_date = snakemake.params.end_date
fn = os.path.abspath(snakemake.output.weather_csv)
evapo = snakemake.params.evapo
multiplier = snakemake.params.multiplier
dry_day_frac = snakemake.params.dry_day_fraction
jpg_out = os.path.join(os.path.split(fn)[0], "weather_plot.jpg")

# now we make some random data, by making some random number between
dates = pd.date_range(start_date, end_date)

print(f"Creating weather data from {start_date} to {end_date}")
# now we create a nice short random precipitation series. Put this into a pandas dataframe and write to a csv file
precip = np.maximum((np.random.rand(len(dates)) - dry_day_frac) * multiplier, 0)
pot_evapo = np.ones(len(dates)) * evapo
df = pd.DataFrame({"p": precip, "ep": pot_evapo}, index=dates)
df.index.name = "date"
print(f"Writing weather data to {fn}")
df.to_csv(fn)

print(f"Plotting weather data in {jpg_out}")
f = plt.figure(figsize=(16, 9))
ax = plt.subplot(111)
df.plot(ax=ax)
f.savefig(jpg_out, bbox_inches="tight", dpi=200)
