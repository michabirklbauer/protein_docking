#!/usr/bin/env python3

# SCORING WORKFLOW - DPP4
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

"""
DESCRIPTION
"""

import os
import json
from shutil import copyfile

if __name__ == '__main__':

    # copy most recent version of PIA to directory
    copyfile("../../../scripts/python/PIA.py", "PIA.py")

    # import PIA
    from PIA import PIA as pa
    from PIA import Preparation as prep
    from PIA import Comparison as paComp
    from PIA import Scoring as paScore

    enabled_1 = False
    enabled_2 = True

    # combine sdf files
    with open("Dpp4_resultsactives.sdf", "r") as f:
        actives = f.read()
        f.close()

    with open("Dpp4_results_decoys.sdf", "r") as f:
        inactives = f.read()
        f.close()

    content = actives + inactives

    with open("Dpp4_results.sdf", "w") as f:
        f.write(content)
        f.close()

    if enabled_1:
        # prepare files
        p = prep()
        # clean PDB file
        pdb = p.remove_ligands("2g5t.pdb", "2g5t_cleaned.pdb")
        # read ligands from docking SDF file
        ligands = p.get_ligands("Dpp4_results.sdf")
        # get ligand names
        sdf_metainfo = p.get_sdf_metainfo("Dpp4_results.sdf")
        ligand_names = sdf_metainfo["names"]
        # write ligands into PDB files
        structures = p.add_ligands_multi("2g5t_cleaned.pdb", "structures", ligands)
        result = pa(structures, ligand_names = ligand_names, poses = "best", path = "current")
        # save results
        with open("structures.json", "w") as f:
                json.dump(result.pdb_entry_results, f)
                f.close()
        with open("structures_best.json", "w") as f:
                json.dump(result.best_pdb_entries, f)
                f.close()

    if enabled_2:
        # Scoring - generate Datasets and feature info
        s = paScore("structures.json", "structures_best.json", is_file = [True, True])
        df_train, df_val, df_test = s.generate_datasets()
        ia_info = s.get_actives_inactives()
        comp_train = s.compare(partition = "train")
        comp_val = s.compare(partition = "val")
        comp_test = s.compare(partition = "test")
        features = s.get_feature_information(filename = "features.csv")
        comp_train.plot("DPP4: Actives vs. Inactives - Training Set", filename = "results/Dpp4_results_comparison_train.jpg")
        comp_train.plot("DPP4: Actives vs. Inactives - Training Set", filename = "results/Dpp4_results_comparison_train.png")
        comp_val.plot("DPP4: Actives vs. Inactives - Validation Set", filename = "results/Dpp4_results_comparison_val.jpg")
        comp_val.plot("DPP4: Actives vs. Inactives - Validation Set", filename = "results/Dpp4_results_comparison_val.png")
        comp_test.plot("DPP4: Actives vs. Inactives - Test Set", filename = "results/Dpp4_results_comparison_test.jpg")
        comp_test.plot("DPP4: Actives vs. Inactives - Test Set", filename = "results/Dpp4_results_comparison_test.png")

    # clean up
    os.remove("Dpp4_results.sdf")
    os.remove("PIA.py")
