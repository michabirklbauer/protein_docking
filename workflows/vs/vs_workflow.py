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

    # copy most recent version of PIA to directory
    copyfile("../../scripts/python/PIA.py", "PIA.py")

    # import PIA
    from PIA import PIA as pa
    from PIA import Preparation as prep
    from PIA import Comparison as paComp

    enabled_1 = True
    enabled_2 = True
    enabled_3 = True
    enabled_4 = True

    if enabled_1:
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

    if enabled_2:
        # PIA - actives
        result = pa(actives_structures, ligand_names = actives_names, poses = "best", path = "current")
        r = result.save("results/sEH_6hgv_results_actives")
        print("Result saved in:")
        print(r)
        r = result.to_csv("results/sEH_6hgv_results_actives_freq.csv")
        result.plot("sEH_6hgv_results_actives", filename = "results/sEH_6hgv_results_actives.jpg")
        result.plot("sEH_6hgv_results_actives", filename = "results/sEH_6hgv_results_actives.png")


    if enabled_3:
        # PIA - inactives
        result = pa(inactives_structures, ligand_names = inactives_names, poses = "best", path = "current")
        r = result.save("results/sEH_6hgv_results_inactives")
        print("Result saved in:")
        print(r)
        r = result.to_csv("results/sEH_6hgv_results_inactives_freq.csv")
        result.plot("sEH_6hgv_results_inactives", filename = "results/sEH_6hgv_results_inactives.jpg")
        result.plot("sEH_6hgv_results_inactives", filename = "results/sEH_6hgv_results_inactives.png")

    if enabled_4:
        # Compare actives vs inactives
        comp = paComp("Actives", "Inactives",
                      "results/sEH_6hgv_results_actives_frequencies.json",
                      "results/sEH_6hgv_results_inactives_frequencies.json",
                      True, True)
        comp.plot("Comparison sEH_6hgv_results.sdf: Actives vs. Inactives", filename = "results/sEH_6hgv_results_comparison.jpg")
        comp.plot("Comparison sEH_6hgv_results.sdf: Actives vs. Inactives", filename = "results/sEH_6hgv_results_comparison.png")

    # clean up
    os.remove("PIA.py")
