#!/usr/bin/env python3

# PLIP ANALYSIS
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import os
import pandas as pd
from shutil import copyfile

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
    result = pa(PDB_ENTRIES, path = "structures")
    r = result.save("results/coagulation_factor_xa")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/coagulation_factor_xa_freq.csv")
    result.plot("Coagulation Factor XA", filename = "results/coagulation_factor_xa.jpg")
    result.plot("Coagulation Factor XA", filename = "results/coagulation_factor_xa.png")

    # clean up
    os.remove("PLIPAnalyzer.py")
