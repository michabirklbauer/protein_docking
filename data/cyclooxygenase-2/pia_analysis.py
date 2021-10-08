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
    r = result.save("results/cyclooxygenase-2")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/cyclooxygenase-2_freq.csv")
    result.plot("Cyclooxygenase-2", filename = "results/cyclooxygenase-2.jpg")
    result.plot("Cyclooxygenase-2", filename = "results/cyclooxygenase-2.png")
    # normalized analysis
    result = pa(PDB_ENTRIES, path = "structures", normalize = True)
    r = result.save("results/cyclooxygenase-2_normalized")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/cyclooxygenase-2_normalized_freq.csv")
    result.plot("Cyclooxygenase-2", filename = "results/cyclooxygenase-2_normalized.jpg")
    result.plot("Cyclooxygenase-2", filename = "results/cyclooxygenase-2_normalized.png")

    # analysis human
    df = pd.read_csv("data_human.csv")
    PDB_ENTRIES = list(df["PDB_ENTRY"])
    print("First 5 PDB entries:")
    print(PDB_ENTRIES[:5])
    result = pa(PDB_ENTRIES, path = "structures", normalize = False)
    r = result.save("results_human/cyclooxygenase-2_human")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results_human/cyclooxygenase-2_human_freq.csv")
    result.plot("Cyclooxygenase-2 (Homo sapiens)", filename = "results_human/cyclooxygenase-2_human.jpg")
    result.plot("Cyclooxygenase-2 (Homo sapiens)", filename = "results_human/cyclooxygenase-2_human.png")
    result = pa(PDB_ENTRIES, path = "structures", normalize = True)
    r = result.save("results_human/cyclooxygenase-2_human_normalized")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results_human/cyclooxygenase-2_human_normalized_freq.csv")
    result.plot("Cyclooxygenase-2 (Homo sapiens)", filename = "results_human/cyclooxygenase-2_human_normalized.jpg")
    result.plot("Cyclooxygenase-2 (Homo sapiens)", filename = "results_human/cyclooxygenase-2_human_normalized.png")

    # clean up
    os.remove("PIA.py")
