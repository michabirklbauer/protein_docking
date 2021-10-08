#!/usr/bin/env python3

# PIA ANALYSIS
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import os
import pandas as pd
from shutil import copyfile

if __name__ == '__main__':

    # copy most recent version of PIA to directory
    copyfile("../../scripts/python/PIA.py", "PIA.py")

    # import PIA
    from PIA import PIA as pa

    # analysis
    df = pd.read_csv("data_train.csv")
    PDB_ENTRIES = list(df["PDB_ENTRY"])
    print("First 5 PDB entries:")
    print(PDB_ENTRIES[:5])
    result = pa(PDB_ENTRIES, path = "structures", normalize = False)
    r = result.save("results/monoamine_oxidase_B")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/monoamine_oxidase_B_freq.csv")
    result.plot("Monoamine Oxidase B", filename = "results/monoamine_oxidase_B.jpg")
    result.plot("Monoamine Oxidase B", filename = "results/monoamine_oxidase_B.png")
    # normalized results
    result = pa(PDB_ENTRIES, path = "structures", normalize = True)
    r = result.save("results/monoamine_oxidase_B_normalized")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/monoamine_oxidase_B_normalized_freq.csv")
    result.plot("Monoamine Oxidase B", filename = "results/monoamine_oxidase_B_normalized.jpg")
    result.plot("Monoamine Oxidase B", filename = "results/monoamine_oxidase_B_normalized.png")

    # clean up
    os.remove("PIA.py")
