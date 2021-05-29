#!/usr/bin/env python3

# PLIP ANALYSIS
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import os
import pandas as pd
from shutil import copyfile

# adjust args accordingly - here as an example for SEH in corresponding dir
if __name__ == '__main__':

    # copy most recent version of PLIPAnalyzer to directory
    copyfile("../../scripts/python/PLIPAnalyzer.py", "PLIPAnalyzer.py")

    # import PLIPAnalyzer
    from PLIPAnalyzer import PLIPAnalyzer as pa

    # analysis
    df = pd.read_csv("data_train.csv")
    PDB_ENTRIES = list(df["PDB_ENTRY"])
    print("First 5 PDB entries:")
    print(PDB_ENTRIES[:5])
    result = pa(PDB_ENTRIES, path = "structures", normalize = False)
    r = result.save("results/soluble_epoxide_hydrolase")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/soluble_epoxide_hydrolase_freq.csv")
    result.plot("Soluble Epoxide Hydrolase", filename = "results/soluble_epoxide_hydrolase.jpg")
    result.plot("Soluble Epoxide Hydrolase", filename = "results/soluble_epoxide_hydrolase.png")
    # normalized results
    result = pa(PDB_ENTRIES, path = "structures", normalize = True)
    r = result.save("results/soluble_epoxide_hydrolase_normalized")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/soluble_epoxide_hydrolase_normalized_freq.csv")
    result.plot("Soluble Epoxide Hydrolase", filename = "results/soluble_epoxide_hydrolase_normalized.jpg")
    result.plot("Soluble Epoxide Hydrolase", filename = "results/soluble_epoxide_hydrolase_normalized.png")

    # clean up
    os.remove("PLIPAnalyzer.py")
