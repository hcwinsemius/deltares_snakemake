import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_markov_rainfall(prev_rainfall, p_01, p_10, thres=1., multiplier=20.):
    """
    Function to generate daily rainfall from previous day rainfall and markov parameters

    Parameters
    ----------
    prev_rainfall : float
        Amount of rainfall from the previous day in [mm day-1]
    p_01 : float(0-1)
        Probability of a wet day following on a dry day
    p_10 : float(0-1)
        Probability of a dry day following on a wet day
    thres : float, optional
        Threshold defining amount of rainfall [mm day-1] that is considered a dry (lt) or wet (ge) day rainfall amount
    multiplier : float, optional
        Multiplier used to convert a randomly drawn number (0-1) into a mm day-1 rainfall amount.
    Returns
    -------
    rainfall : float
        Amount of rainfall at the present day

    """
    if prev_rainfall >= thres:
        wet = np.random.choice(np.array([True, False]), p=np.array([1-p_10, p_10]))
    else:
        wet = np.random.choice(np.array([True, False]), p=np.array([p_01, 1-p_01]))
    if wet:
        return np.random.rand() * multiplier  # can be under thres, which makes this into a rainfall day that can transition more easily to a dry day.
    else:
        return 0.


# SnakeMake arguments are automatically passed WITHOUT importing a lib. So this .py file is basicallly used as a template.
start_date = snakemake.params.start_date
end_date = snakemake.params.end_date
fn = os.path.abspath(snakemake.output.weather_csv)
evapo = snakemake.params.evapo
multiplier = snakemake.params.multiplier

# probabilities from dry to wet and wet to dry
p_01 = snakemake.params.p_01
p_10 = snakemake.params.p_10
thres = snakemake.params.rainfall_thres

# dry_day_frac = snakemake.params.dry_day_fraction
jpg_out = os.path.join(os.path.split(fn)[0], "weather_plot.jpg")

# now we make some random data, by making some random number between
dates = pd.date_range(start_date, end_date)

print(f"Creating weather data from {start_date} to {end_date}")
# now we create a nice precipitation series. Put this into a pandas dataframe and write to a csv file
precip = [0.]
for n in range(len(dates)-1):
    precip.append(
        get_markov_rainfall(
            precip[-1],
            p_01,
            p_10,
            thres=thres,
            multiplier=multiplier
        )
    )
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
