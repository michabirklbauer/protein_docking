#!/usr/bin/env python3

# SDF WORKFLOW
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

"""
DESCRIPTION
Examplary workflow for cleaning a PDB file from existing ligands and calculate
frequencies of interactions with ligands from a docking result [which was saved
in SDF file format].
"""

import os
from shutil import copyfile

if __name__ == '__main__':

    # copy most recent version of PLIPAnalyzer to directory
    copyfile("../../scripts/python/PLIPAnalyzer.py", "PLIPAnalyzer.py")

    # import PLIPAnalyzer
    from PLIPAnalyzer import PLIPAnalyzer as pa
    from PLIPAnalyzer import Preparation as prep

    # analysis
    p = prep()
    # clean PDB file
    pdb = p.remove_ligands("6hgv.pdb", "6hgv_cleaned.pdb")
    # read ligands from docking SDF file
    ligands = p.get_ligands("results_vs_6hgv_6A_Gold.sdf")
    # write ligands into PDB files
    structures = p.add_ligands_multi("6hgv_cleaned.pdb", "structures", ligands)
    # PLIPAnalyzer
    result = pa(structures, path = "current")
    r = result.save("results/results_vs_6hgv_6A_Gold")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/results_vs_6hgv_6A_Gold_freq.csv")
    result.plot("results_vs_6hgv_6A_Gold", filename = "results/results_vs_6hgv_6A_Gold.jpg")
    result.plot("results_vs_6hgv_6A_Gold", filename = "results/results_vs_6hgv_6A_Gold.png")

    # clean up
    os.remove("PLIPAnalyzer.py")
