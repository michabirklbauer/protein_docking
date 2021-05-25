#!/usr/bin/env python3

# VS WORKFLOW
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

"""
DESCRIPTION
Examplary workflow for cleaning a PDB file from existing ligands and calculate
frequencies of interactions with ligands from a docking result [which was saved
in SDF file format] and comparison between active and inactive ligands.
"""

import os
import json
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
    ligands = p.get_ligands("sEH_6hgv_results.sdf")
    # get ligand names
    sdf_metainfo = p.get_sdf_metainfo("sEH_6hgv_results.sdf")
    ligand_names = sdf_metainfo["names"]
    # write ligands into PDB files
    structures = p.add_ligands_multi("6hgv_cleaned.pdb", "structures", ligands)
    # get subsets
    actives_idx, inactives_idx = p.actives_inactives_split("sEH_6hgv_results.sdf")
    actives_structures = [structures[i] for i in actives_idx]
    actives_names = [ligand_names[i] for i in actives_idx]
    inactives_structures = [structures[i] for i in inactives_idx]
    inactives_names = [ligand_names[i] for i in inactives_idx]

    # save subsets
    with open("subsets.json", "w") as f:
        json.dump({"actives_idx": actives_idx,
                   "actives_structures": actives_structures,
                   "actives_names": actives_names,
                   "inactives_idx": inactives_idx,
                   "inactives_structures": inactives_structures,
                   "inactives_names": inactives_names}, f)
        f.close()

    # PLIPAnalyzer - actives
    result = pa(actives_structures, ligand_names = actives_names, poses = "best", path = "current")
    r = result.save("results/sEH_6hgv_results_actives")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/sEH_6hgv_results_actives_freq.csv")
    result.plot("sEH_6hgv_results_actives", filename = "results/sEH_6hgv_results_actives.jpg")
    result.plot("sEH_6hgv_results_actives", filename = "results/sEH_6hgv_results_actives.png")


    # PLIPAnalyzer - inactives
    result = pa(inactives_structures, ligand_names = inactives_names, poses = "best", path = "current")
    r = result.save("results/sEH_6hgv_results_inactives")
    print("Result saved in:")
    print(r)
    r = result.to_csv("results/sEH_6hgv_results_inactives_freq.csv")
    result.plot("sEH_6hgv_results_inactives", filename = "results/sEH_6hgv_results_inactives.jpg")
    result.plot("sEH_6hgv_results_inactives", filename = "results/sEH_6hgv_results_inactives.png")

    # clean up
    os.remove("PLIPAnalyzer.py")
