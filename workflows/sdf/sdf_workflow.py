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
import json
from shutil import copyfile

if __name__ == '__main__':

    # copy most recent version of PIA to directory
    copyfile("../../scripts/python/PIA.py", "PIA.py")

    # import PIA
    from PIA import PIA as pa
    from PIA import Preparation as prep

    enabled_1 = True
    enabled_2 = True
    enabled_3 = True
    enabled_4 = True

    # analysis
    p = prep()
    # clean PDB file
    pdb = p.remove_ligands("6hgv.pdb", "6hgv_cleaned.pdb")
    # read ligands from docking SDF file
    ligands = p.get_ligands("results_vs_6hgv_6A_Gold.sdf")
    # get ligand names
    sdf_metainfo = p.get_sdf_metainfo("results_vs_6hgv_6A_Gold.sdf")
    ligand_names = sdf_metainfo["names"]
    # write ligands into PDB files
    structures = p.add_ligands_multi("6hgv_cleaned.pdb", "structures", ligands)

    if enabled_1:
        # PIA
        result = pa(structures, path = "current", normalize = False)
        r = result.save("results/results_vs_6hgv_6A_Gold")
        print("Result saved in:")
        print(r)
        r = result.to_csv("results/results_vs_6hgv_6A_Gold_freq.csv")
        result.plot("results_vs_6hgv_6A_Gold", filename = "results/results_vs_6hgv_6A_Gold.jpg")
        result.plot("results_vs_6hgv_6A_Gold", filename = "results/results_vs_6hgv_6A_Gold.png")

    if enabled_2:
        # PIA - normalized
        result = pa(structures, path = "current")
        r = result.save("results/results_normalized_vs_6hgv_6A_Gold")
        print("Result saved in:")
        print(r)
        r = result.to_csv("results/results_normalized_vs_6hgv_6A_Gold_freq.csv")
        result.plot("results_vs_6hgv_6A_Gold", filename = "results/results_normalized_vs_6hgv_6A_Gold.jpg")
        result.plot("results_vs_6hgv_6A_Gold", filename = "results/results_normalized_vs_6hgv_6A_Gold.png")

    if enabled_3:
        # PIA - best poses only
        result = pa(structures, ligand_names = ligand_names, poses = "best", path = "current", normalize = False)
        r = result.save("results/results_best_vs_6hgv_6A_Gold")
        print("Result saved in:")
        print(r)
        r = result.to_csv("results/results_best_vs_6hgv_6A_Gold_freq.csv")
        result.plot("results_best_vs_6hgv_6A_Gold", filename = "results/results_best_vs_6hgv_6A_Gold.jpg")
        result.plot("results_best_vs_6hgv_6A_Gold", filename = "results/results_best_vs_6hgv_6A_Gold.png")
        # save attributes to json
        with open("pdb_entry_results.json", "w") as f:
            json.dump(result.pdb_entry_results, f)
            f.close()
        with open("best_pdb_entries.json", "w") as f:
            json.dump({"best_pdb_entries": result.best_pdb_entries}, f)
            f.close()

    if enabled_4:
        # PIA - best poses only, normalized
        result = pa(structures, ligand_names = ligand_names, poses = "best", path = "current")
        r = result.save("results/results_best_normalized_vs_6hgv_6A_Gold")
        print("Result saved in:")
        print(r)
        r = result.to_csv("results/results_best_normalized_vs_6hgv_6A_Gold_freq.csv")
        result.plot("results_best_vs_6hgv_6A_Gold", filename = "results/results_best_normalized_vs_6hgv_6A_Gold.jpg")
        result.plot("results_best_vs_6hgv_6A_Gold", filename = "results/results_best_normalized_vs_6hgv_6A_Gold.png")

    # clean up
    os.remove("PIA.py")
