import numpy as np
import pandas as pd
import copy
import os

def hyd_model(storage_1, storage_2, precip, pot_evapo, interception, storage_1_max, alpha, k_fast, k_slow):
    """
    Simple bucket model with two time transient stores (storage_1, storage_2), two forcings (precip and pot_evapo)
    and 5 parameters (interception, storage_1_max, alpha, k_fast, k_slow)
    """
    # explicit solution defined below

    # store the previous time step storage for water balance computation
    _storage_2 = copy.deepcopy(storage_2)
    _storage_1 = copy.deepcopy(storage_1)

    # reduce rainfall with interception threshold
    p_eff = np.maximum(precip - interception, 0)
    e_I = precip - p_eff
    # potential transpiraiton
    tp = pot_evapo - e_I
    # separate rainfall into upper and lower store
    p_to_storage_1 = alpha * p_eff
    p_to_storage_2 = (1 - alpha) * p_eff
    
    # first treat lower store
    storage_2 += p_to_storage_2
    # outflow from storage_2
    Q2 = storage_2 / k_slow
    storage_2 -= Q2
    
    # now treat upper store
    storage_1 += p_to_storage_1
    Q1_overflow = np.maximum(storage_1 - storage_1_max, 0)
    storage_1 -= Q1_overflow

    # calculate transpiration as linear function of storage relative to max storage
    t = storage_1 / storage_1_max * tp
    storage_1 -= t
    # total evapo
    E = t + e_I

    # calculate outflow from drainage
    Q1_outflow = storage_1 / k_fast
    storage_1 -= Q1_outflow
    # total flow
    Q = Q1_outflow + Q1_overflow + Q2
    wb = (_storage_1 - storage_1) + (_storage_2 - storage_2) + p - E - Q
    return storage_1, storage_2, Q, E, wb

# SnakeMake arguments are automatically passed WITHOUT importing a lib. So this .py file is basicallly used as a template.
input_csv = snakemake.input.input_csv
output_csv = os.path.abspath(snakemake.output.output_csv)

# get the parameters needed to run hydrological model, keep them as dict
param_set = snakemake.params.param_set
param_names = snakemake.params.param_names

# make a nice dict of pars
params = {k: v for k, v in zip(param_names, param_set)}

print(params)

# get the storage_1 and storage_2 from params
storage_1 = params["storage_1"]
storage_2 = params["storage_2"]
del params["storage_1"]
del params["storage_2"]

# read dataframe

print(f"Reading forcing from {input_csv}")
df = pd.read_csv(input_csv)

# go through all rows
Qs  =[]
Es = []
wbs = []
print(f"Running model with the following parameters: {params}")
for n, r in df.iterrows():
    precip, pot_evapo = r["p"], r["ep"]
#     print(type(p), type(ep))
    storage_1, storage_2, Q, E, wb = hyd_model(
        storage_1,
        storage_2,
        precip,
        pot_evapo,
        **params
    )
    Qs.append(Q)
    Es.append(E)
    wbs.append(wb)

# convert outputs in numpy arrays
Qs = np.array(Qs)
wbs = np.array(wbs)
Es = np.array(Es)

# make a nice pandas dataframe
df = pd.DataFrame({"Q": Qs, "E": Es, "wb": wbs}, index=df.index.dates)
df.index.name = "date"

# store in csv output file
print(f"Writing outputs to {output_csv}")
df.to_csv(output_csv)

