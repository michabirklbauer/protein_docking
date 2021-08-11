#!/usr/bin/env python3

# SCORING WORKFLOW - VITAMIN E WORKFLOW6.1
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

    # copy most recent version of PLIPAnalyzer to directory
    copyfile("../../../scripts/python/PLIPAnalyzer.py", "PLIPAnalyzer.py")

    # import PLIPAnalyzer
    from PLIPAnalyzer import PLIPAnalyzer as pa
    from PLIPAnalyzer import Preparation as prep
    from PLIPAnalyzer import Comparison as paComp
    from PLIPAnalyzer import Scoring as paScore

    enabled_1 = True
    enabled_2 = True

    if enabled_1:
        # prepare files
        p = prep()
        # clean PDB file
        pdb = p.remove_ligands("6ncf.pdb", "6ncf_cleaned.pdb")
        # read ligands from docking SDF file
        ligands_info = p.get_labeled_names("Vitamin_E_6_bearbeitet.sdf", ">=", 1000)
        ligands = ligands_info["ligands"]
        # get labeled ligand names
        ligand_names = ligands_info["names"]
        # write ligands into PDB files
        structures = p.add_ligands_multi("6ncf_cleaned.pdb", "structures", ligands)
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
        comp_train.plot("Vitamin E - W6.1: Actives vs. Inactives - Training Set", filename = "results/vite_6ncf_results_comparison_train.jpg")
        comp_train.plot("Vitamin E - W6.1: Actives vs. Inactives - Training Set", filename = "results/vite_6ncf_results_comparison_train.png")
        comp_val.plot("Vitamin E - W6.1: Actives vs. Inactives - Validation Set", filename = "results/vite_6ncf_results_comparison_val.jpg")
        comp_val.plot("Vitamin E - W6.1: Actives vs. Inactives - Validation Set", filename = "results/vite_6ncf_results_comparison_val.png")
        comp_test.plot("Vitamin E - W6.1: Actives vs. Inactives - Test Set", filename = "results/vite_6ncf_results_comparison_test.jpg")
        comp_test.plot("Vitamin E - W6.1: Actives vs. Inactives - Test Set", filename = "results/vite_6ncf_results_comparison_test.png")

    # clean up
    os.remove("PLIPAnalyzer.py")
