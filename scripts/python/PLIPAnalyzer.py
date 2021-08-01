#!/usr/bin/env python3

# PLIP ANALYZER
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

version = "0.5.1"
date = "20210709"

import json
import warnings
import numpy as np
import pandas as pd
import traceback as tb
from rdkit import Chem
from biopandas.pdb import PandasPdb
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from plip.structure.preparation import PDBComplex
from plip.basic.config import biolip_list

#### ----------------------------- COFACTORS ------------------------------ ####

# Common cofactors found in PDB.
# List taken from:
# http://www.ebi.ac.uk/thornton-srv/databases/CoFactor/
# excluded = ["Coenzyme M", "Factor F430"]
# 3 letter codes from:
# https://www.ebi.ac.uk/pdbe-srv/pdbechem/
# Notes: Heme? (--> see hemes), Menaquinone? (MQ7 to MQ9 included), MIO?,
#        Molybdopterin? (--> see molybdopterins), Orthoquinone residues?,
#        Ubiquinone? (UQ1 to U10 included)
cofactors = ["B1Z", "ASC", "BIO", "BTN", "COA", "COZ", "TP7", "DPM", "FAD", "FMN",
             "GSH", "LPA", "MQ7", "MQ8", "MQ9", "NAD", "PNS", "PLP", "PQQ", "SAM",
             "THF", "TPP", "TPQ", "U10", "UQ1", "UQ2", "UQ7", "UQ8", "UQ9"]

# Heme cofactors
# list of shortcodes of heme-like structures from:
# https://www.ebi.ac.uk/pdbe-srv/pdbechem/
# query: "Molecule Name like Heme"
hemes = ["1FH", "2FH", "522", "89R", "DDH", "DHE", "HAS", "HDD", "HDE", "HDM",
         "HEA", "HEB", "HEC", "HEM", "HEO", "HES", "HEV", "HP5", "MH0", "MHM",
         "N7H", "NTE", "OBV", "SRM", "VER", "VOV"]

# Molybdopterin cofactors
# list of shortcodes of molybdopterin-like structures from:
# https://www.ebi.ac.uk/pdbe-srv/pdbechem/
# query: "Molecule Name like Molybdopterin"
molybdopterins = ["2MD", "MGD", "MSS", "MTQ", "MTV", "PCD"]

# a list of userdefined cofactors (should be extended as needed)
# entries from data/cofactor_info.yaml
userdef_cofactors = ["NAP", "NDP", "COH"]

exclusion_list = biolip_list + cofactors + hemes + molybdopterins + userdef_cofactors

#### ------------------------------ CLASSES ------------------------------- ####

class NoInteractionsWarning(UserWarning):
    """
    -- DESCRIPTION --
    Displayed if no interactions are found during analysis.
    """
    pass

class CalculationErrorWarning(UserWarning):
    """
    -- DESCRIPTION --
    Displayed if one or more but not all of the results are empty.
    """
    pass

class NotInPartitionWarning(UserWarning):
    """
    -- DESCRIPTION --
    Displayed if an entry or interaction is found that does not belong to any
    data partition.
    """
    pass

class ParseError(RuntimeError):
    """
    -- DESCRIPTION --
    Raised if file could not be parsed properly.
    """
    pass

class Preparation:
    """to be implemented"""

    # empty constructor
    def __init__(self):
        pass

    # get molecule names from sdf file -- safer to use get_sdf_metainfo() instead!
    def get_sdf_names(self,
                      sdf_file):

        """
        -- DESCRIPTION --
        """

        with open(sdf_file, "r") as f:
            content = f.read()
            f.close()

        names = []
        mols = content.split("$$$$")
        for mol in mols:
            mol_cleaned = mol.strip()
            if mol_cleaned != "":
                names.append(mol_cleaned.split("\n")[0].strip())

        return names

    # get GOLD fitness from sdf file -- safer to use get_sdf_metainfo() instead!
    def get_sdf_fitness(self,
                        sdf_file):

        """
        -- DESCRIPTION --
        """

        with open(sdf_file, "r") as f:
            content = f.read()
            f.close()

        scores = []
        mols = content.split("$$$$")
        for mol in mols:
            mol_cleaned = mol.strip()
            if mol_cleaned != "":
                tmp = mol_cleaned.split("> <Gold.Goldscore.Fitness>")[1]
                score = float(tmp.split(">")[0])
                scores.append(score)

        return scores

    # get GOLD IC50 5-LO values
    def get_sdf_IC50s(self,
                      sdf_file,
                      NA = 1000.0):

        """
        -- DESCRIPTION --
        """

        with open(sdf_file, "r") as f:
            content = f.read()
            f.close()

        ic50s = []
        NA_counter = 0
        mols = content.split("$$$$")
        for mol in mols:
            mol_cleaned = mol.strip()
            if mol_cleaned != "":
                tmp = mol_cleaned.split("> <IC50 5-LO>")[1]
                try:
                    ic50 = float(tmp.split(">")[0])
                except ValueError as e:
                    NA_counter = NA_counter + 1
                    ic50 = NA
                ic50s.append(ic50)

        if NA_counter != 0:
            print("NAs encountered in " + NA_counter + " molecules. IC50 value substituted as " + str(NA) + " (defined in get_sdf_IC50s).")
        return ic50s

    # get names and GOLD fitness from sdf file
    def get_sdf_metainfo(self,
                         sdf_file):

        """
        -- DESCRIPTION --
        """

        names = self.get_sdf_names(sdf_file)
        scores = self.get_sdf_fitness(sdf_file)

        if len(names) == len(scores):
            return {"names": names, "fitness": scores}
        else:
            raise ParseError("SDF file could not be parsed [nr of names != nr of scores].")

    # label names as active or decoy depending on IC50
    def get_labeled_names(self,
                          sdf_file,
                          condition_operator = "==",
                          condition_value = 1000):

        """
        -- DESCRIPTION --
        """

        names = self.get_sdf_names(sdf_file)
        ic50s = self.get_sdf_IC50s(sdf_file)

        new_names = []

        if len(names) != len(ic50s):
            raise ParseError("SDF file could not be parsed [nr of names != nr of IC50s].")
        else:
            for i, name in enumerate(names):
                name_split = name.split("|")
                if condition_operator == "==":
                    if ic50s[i] == condition_value:
                        new_name = name_split[0] + "|" + name_split[1] + "_decoy|" + "|".join(name_split[2:])
                    else:
                        new_name = name_split[0] + "|" + name_split[1] + "_active|" + "|".join(name_split[2:])
                elif condition_operator == "<=":
                    if ic50s[i] <= condition_value:
                        new_name = name_split[0] + "|" + name_split[1] + "_decoy|" + "|".join(name_split[2:])
                    else:
                        new_name = name_split[0] + "|" + name_split[1] + "_active|" + "|".join(name_split[2:])
                elif condition_operator == "<":
                    if ic50s[i] < condition_value:
                        new_name = name_split[0] + "|" + name_split[1] + "_decoy|" + "|".join(name_split[2:])
                    else:
                        new_name = name_split[0] + "|" + name_split[1] + "_active|" + "|".join(name_split[2:])
                elif condition_operator == ">=":
                    if ic50s[i] >= condition_value:
                        new_name = name_split[0] + "|" + name_split[1] + "_decoy|" + "|".join(name_split[2:])
                    else:
                        new_name = name_split[0] + "|" + name_split[1] + "_active|" + "|".join(name_split[2:])
                elif condition_operator == ">":
                    if ic50s[i] > condition_value:
                        new_name = name_split[0] + "|" + name_split[1] + "_decoy|" + "|".join(name_split[2:])
                    else:
                        new_name = name_split[0] + "|" + name_split[1] + "_active|" + "|".join(name_split[2:])
                elif condition_operator == "!=":
                    if ic50s[i] != condition_value:
                        new_name = name_split[0] + "|" + name_split[1] + "_decoy|" + "|".join(name_split[2:])
                    else:
                        new_name = name_split[0] + "|" + name_split[1] + "_active|" + "|".join(name_split[2:])
                else:
                    new_name = name
                new_names.append(new_name)

            return new_names

    # get indices of active and inactive molecules in sdf file
    def actives_inactives_split(self,
                                sdf_file):

        """
        -- DESCRIPTION --
        """

        actives = []
        inactives = []
        names = self.get_sdf_metainfo(sdf_file)["names"]

        for i, name in enumerate(names):
            if "inactive" in name or "decoy" in name:
                inactives.append(i)
            else:
                actives.append(i)

        return [actives, inactives]

    # remove ligands from a pdb file
    def remove_ligands(self,
                       input_file,
                       output_file):

        """
        -- DESCRIPTION --
        """

        mol = PDBComplex()
        mol.load_pdb(input_file)
        for ligand in mol.ligands:
            mol.characterize_complex(ligand)

        ligands = [str(lig.split(":")[0]).strip() for lig in mol.interaction_sets.keys()]

        pdb = PandasPdb().read_pdb(input_file)

        def remove_conect(pdb, residue_name):

            atoms = pdb.df["HETATM"][pdb.df["HETATM"]["residue_name"] == residue_name]["atom_number"]
            atoms_str = [str(atom) for atom in atoms]

            df = pdb.df["OTHERS"]
            idx = []
            for index, row in df.iterrows():
                if row["record_name"] == "CONECT":
                    if any(a in row["entry"] for a in atoms_str):
                        idx.append(index)
            return idx

        def remove_hetatm(pdb, residue_name):
            df = pdb.df["HETATM"]
            idx = []
            for index, row in df.iterrows():
                if residue_name in str(row["residue_name"]):
                    idx.append(index)
            return idx

        for ligand in ligands:
            pdb.df["OTHERS"].drop(remove_conect(pdb, ligand), inplace = True)
            pdb.df["HETATM"].drop(remove_hetatm(pdb, ligand), inplace = True)

        pdb.to_pdb(output_file)

        return pdb

    # read ligands from an sdf file
    def get_ligands(self, sdf_file):

        """
        -- DESCRIPTION --
        """

        return [ x for x in Chem.SDMolSupplier(sdf_file)]

    # add ligands to a pdb structure
    def add_ligands(self,
                    input_file,
                    output_file,
                    ligands,
                    ligand_IDs = None,
                    ligand_chain_IDs = None,
                    verbose = 1):

        """
        -- DESCRIPTION --
        """

        if ligand_IDs is None:
            ligand_IDs = ["LG" + str(x) for x in range(1, len(ligands) + 1)]

        if ligand_chain_IDs is None:
                ligand_chain_ids = ["A" for x in ligands]

        pdb = Chem.MolFromPDBFile(input_file)

        skipped_ligands = []
        for i, ligand in enumerate(ligands):

            # try merging pdb structure and ligand
            try:
                merged = Chem.CombineMols(pdb, ligand)
            except Exception as e:
                print("Error merging structure with ligand " + ligand_IDs[i] + "! Skipping...")
                skipped_ligands.append(str(ligand_IDs[i]))
                tb.print_exc()
                continue

            Chem.MolToPDBFile(merged, output_file)

            if verbose:
                print("Adding ligand " + ligand_IDs[i] + " to file " + output_file + "!")

            pdb = PandasPdb().read_pdb(output_file)
            pdb.df["HETATM"].loc[pdb.df["HETATM"]["residue_name"] == "UNL", "chain_id"] = ligand_chain_ids[i]
            pdb.df["HETATM"].replace("UNL", ligand_IDs[i], inplace = True)
            pdb.to_pdb(output_file)

            pdb = Chem.MolFromPDBFile(output_file)

        if len(skipped_ligands) > 0:
            print("Number of skipped ligands: " + str(len(skipped_ligands)))
            print("Skipped ligand IDs: " + ", ".join(skipped_ligands))
            if len(skipped_ligands) == len(ligands):
                print("No ligand could be added to the structure! No file was created!")
                return 2
            else:
                return 1
        else:
            return 0

    # write ligands to multiple files
    def add_ligands_multi(self,
                          input_file,
                          output_path,
                          ligands,
                          ligands_per_file = 1,
                          **kwargs):

        """
        -- DESCRIPTION --
        """

        output_files = []
        for i in range(0, len(ligands), ligands_per_file):
            l_idx = i
            r_idx = i + ligands_per_file
            current_ligands = ligands[l_idx:r_idx]
            if output_path not in ["", ".", "current"]:
                output_file = output_path + "/" + input_file.split(".")[0]
            else:
                output_file = input_file.split(".")[0]
            output_file = output_file + str(int(i / ligands_per_file + 1)) + ".pdb"
            output_files.append(output_file)
            self.add_ligands(input_file, output_file, current_ligands, **kwargs)

        return output_files

class Comparison:
    """
    -- DESCRIPTION --
    """

    name_a = None
    name_b = None
    frequencies_a = None
    frequencies_b = None
    comparison = None

    def __init__(self,
                 name_a,
                 name_b,
                 frequencies_a,
                 frequencies_b,
                 is_file_a = False,
                 is_file_b = False):

        """
        -- DESCRIPTION --
        """

        if is_file_a:
            with open(frequencies_a, "r", encoding = "utf-8") as f:
                frequencies_a = json.load(f)
                f.close()

        if is_file_b:
            with open(frequencies_b, "r", encoding = "utf-8") as f:
                frequencies_b = json.load(f)
                f.close()

        self.name_a = name_a
        self.name_b = name_b
        self.frequencies_a = frequencies_a
        self.frequencies_b = frequencies_b

        keys_a = set(frequencies_a.keys())
        keys_b = set(frequencies_b.keys())
        keys = keys_a.union(keys_b)

        comparison = {}
        for key in keys:
            if key in frequencies_a and key in frequencies_b:
                comparison[key] = {"difference": abs(frequencies_a[key] - frequencies_b[key]),
                                   name_a: frequencies_a[key],
                                   name_b: frequencies_b[key]}
            elif key in frequencies_a and key not in frequencies_b:
                comparison[key] = {"difference": abs(frequencies_a[key] - 0),
                                   name_a: frequencies_a[key],
                                   name_b: 0}
            elif key in frequencies_b and key not in frequencies_a:
                comparison[key] = {"difference": abs(0 - frequencies_b[key]),
                                   name_a: 0,
                                   name_b: frequencies_b[key]}
            else:
                warnings.warn("Oops! Key not in either dict! This should not happen!", UserWarning)

        self.comparison = dict(sorted(comparison.items(), key = lambda x: x[1]["difference"], reverse = True))

    def plot(self,
             title = None,
             filename = None,
             width = 20,
             height = 5,
             sig_digits = 2,
             label_offset = 0.05):

        """
        -- DESCRIPTION --
        """

        if title is None:
            title = "Comparison: " + self.name_a + " vs. " + self.name_b

        values_a = []
        values_b = []

        for key in self.comparison:
            values_a.append(self.comparison[key][self.name_a])
            values_b.append(self.comparison[key][self.name_b])

        fig = plt.figure(figsize = (width, height))
        ax = fig.add_axes([0,0,1,1])
        x = np.arange(len(values_a))
        ax.bar(x - 0.2, values_a, 0.4)
        ax.bar(x + 0.2, values_b, 0.4)
        plt.xticks(x, list(self.comparison.keys()), rotation = "vertical")
        xlocs, xlabs = plt.xticks()
        for i, v in enumerate(values_a):
            plt.text(xlocs[i] - 0.2, v + label_offset, str(round(v, sig_digits)), horizontalalignment = "center", rotation = 90)
        for i, v in enumerate(values_b):
            plt.text(xlocs[i] + 0.2, v + label_offset, str(round(v, sig_digits)), horizontalalignment = "center", rotation = 90)
        plt.title(title)
        plt.xlabel("Interaction")
        plt.ylabel("Relative Frequency")
        plt.legend([self.name_a, self.name_b])
        if filename is not None:
            fig.savefig(filename, bbox_inches = "tight", dpi = 150)
        plt.show()

        return fig

class Scoring:
    """
    -- DESCRIPTION --
    """

    pdb_entry_results = None
    best_pdb_entries = None

    # internal train/val/test keys - set by generate_datasets()
    entries_train_keys = None
    entries_val_keys = None
    entries_test_keys = None

    # constructor to set input files
    def __init__(self,
                 pdb_entry_results,
                 best_pdb_entries,
                 is_file = [False, False]):

        """
        -- DESCRIPTION --
        """

        if is_file[0]:
            with open(pdb_entry_results, "r", encoding="utf-8") as f:
                self.pdb_entry_results = json.load(f)
                f.close()
        else:
            self.pdb_entry_results = pdb_entry_results

        if is_file[1]:
            with open(best_pdb_entries, "r", encoding="utf-8") as f:
                self.best_pdb_entries = json.load(f)
                f.close()
        else:
            self.best_pdb_entries = best_pdb_entries

    # generate training, validation and test data
    def generate_datasets(self,
                          train_output = "data_train.csv",
                          val_output = "data_val.csv",
                          test_output = "data_test.csv"):

        """
        -- DESCRIPTION --
        """

        # dummy label encoder
        def get_label(input):
            if "inactive" in input or "decoy" in input:
                return "inactive"
            else:
                return "active"

        # extract ligand names
        ligand_names = []
        for entry in self.pdb_entry_results:
            ligand_names.append(self.pdb_entry_results[entry]["ligand_name"])

        # filter for unique ligand names
        u_names = []
        for l_name in ligand_names:
            n = "|".join(l_name.split("|")[:-1])
            u_names.append(n)
        unique_names = list(set(u_names))

        # split ligand names by train_size in train, val, test
        names_train_val, names_test = train_test_split(unique_names, train_size = 0.8, random_state = 42, shuffle = True)
        names_train, names_val = train_test_split(names_train_val, train_size = 0.8, random_state = 1337, shuffle = True)

        # get train, val and test entries
        entries_train = []
        entries_train_keys = []
        entries_val = []
        entries_val_keys = []
        entries_test = []
        entries_test_keys = []
        for entry in self.best_pdb_entries:
            l_name = self.pdb_entry_results[entry]["ligand_name"]
            l_name_short = "|".join(l_name.split("|")[:-1])
            if l_name_short in names_train:
                entries_train.append(self.pdb_entry_results[entry])
                entries_train_keys.append(entry)
            elif l_name_short in names_val:
                entries_val.append(self.pdb_entry_results[entry])
                entries_val_keys.append(entry)
            elif l_name_short in names_test:
                entries_test.append(self.pdb_entry_results[entry])
                entries_test_keys.append(entry)
            else:
                warnings.warn("Oops! " + l_name_short + " not in any data partition! This should not happen!", NotInPartitionWarning)

        # set internal train/val/test keys
        self.entries_train_keys = entries_train_keys
        self.entries_val_keys = entries_val_keys
        self.entries_test_keys = entries_test_keys

        # get interactions in the training partition
        training_interactions = []
        for entry in entries_train:
            interactions = entry["interactions"]
            for hc in interactions["Hydrophobic_Contacts"]:
                training_interactions.append("Hydrophobic_Interaction:" + hc)
            for sb in interactions["Salt_Bridges"]:
                training_interactions.append("Salt_Bridge:" + sb)
            for hb in interactions["Hydrogen_Bonds"]:
                training_interactions.append("Hydrogen_Bond:" + hb)
            for ps in interactions["Pi_Stacking"]:
                training_interactions.append("Pi-Stacking:" + ps)
            for pc in interactions["Pi_Cation_Interactions"]:
                training_interactions.append("Pi-Cation_Interaction:" + pc)
            for hab in interactions["Halogen_Bonds"]:
                training_interactions.append("Halogen_Bond:" + hab)
            for wb in interactions["Water_Bridges"]:
                training_interactions.append("Water_Bridge:" + wb)
            for mc in interactions["Metal_Complexes"]:
                training_interactions.append("Metal_Complex:" + mc)
        interactions_train = list(set(training_interactions))

        # -- generate training dataset --
        data_train = {
            "INDEX": list(range(1, len(entries_train) + 1)),
            "NAME": ["dummy_name" for entry in entries_train]
            }

        for i in interactions_train:
            data_train[i] = [0 for entry in entries_train]

        for i, entry in enumerate(entries_train):
            data_train["NAME"][i] = "|".join(entry["ligand_name"].split("|")[:-1])
            interactions = entry["interactions"]
            for hc in interactions["Hydrophobic_Contacts"]:
                data_train["Hydrophobic_Interaction:" + hc][i] = data_train["Hydrophobic_Interaction:" + hc][i] + 1
            for sb in interactions["Salt_Bridges"]:
                data_train["Salt_Bridge:" + sb][i] = data_train["Salt_Bridge:" + sb][i] + 1
            for hb in interactions["Hydrogen_Bonds"]:
                data_train["Hydrogen_Bond:" + hb][i] = data_train["Hydrogen_Bond:" + hb][i] + 1
            for ps in interactions["Pi_Stacking"]:
                data_train["Pi-Stacking:" + ps][i] = data_train["Pi-Stacking:" + ps][i] + 1
            for pc in interactions["Pi_Cation_Interactions"]:
                data_train["Pi-Cation_Interaction:" + pc][i] = data_train["Pi-Cation_Interaction:" + pc][i] + 1
            for hab in interactions["Halogen_Bonds"]:
                data_train["Halogen_Bond:" + hab][i] = data_train["Halogen_Bond:" + hab][i] + 1
            for wb in interactions["Water_Bridges"]:
                data_train["Water_Bridge:" + wb][i] = data_train["Water_Bridge:" + wb][i] + 1
            for mc in interactions["Metal_Complexes"]:
                data_train["Metal_Complex:" + mc][i] = data_train["Metal_Complex:" + mc][i] + 1

        df_train = pd.DataFrame(data_train)

        df_train["LABEL"] = df_train["NAME"].apply(lambda x: get_label(x))

        if df_train.loc[df_train["NAME"] == "dummy_name"].shape[0] == 0:
            df_train.to_csv(train_output, index = False)
        else:
            warnings.warn("Found dummy variables in dataset! Not saving to " + train_output, UserWarning)

        # -- generate validation dataset --
        data_val = {
            "INDEX": list(range(1, len(entries_val) + 1)),
            "NAME": ["dummy_name" for entry in entries_val]
            }

        for i in interactions_train:
            data_val[i] = [0 for entry in entries_val]

        for i, entry in enumerate(entries_val):
            data_val["NAME"][i] = "|".join(entry["ligand_name"].split("|")[:-1])
            interactions = entry["interactions"]
            for hc in interactions["Hydrophobic_Contacts"]:
                key = "Hydrophobic_Interaction:" + hc
                if key in data_val:
                    data_val[key][i] = data_val[key][i] + 1
            for sb in interactions["Salt_Bridges"]:
                key = "Salt_Bridge:" + sb
                if key in data_val:
                    data_val[key][i] = data_val[key][i] + 1
            for hb in interactions["Hydrogen_Bonds"]:
                key = "Hydrogen_Bond:" + hb
                if key in data_val:
                    data_val[key][i] = data_val[key][i] + 1
            for ps in interactions["Pi_Stacking"]:
                key = "Pi-Stacking:" + ps
                if key in data_val:
                    data_val[key][i] = data_val[key][i] + 1
            for pc in interactions["Pi_Cation_Interactions"]:
                key = "Pi-Cation_Interaction:" + pc
                if key in data_val:
                    data_val[key][i] = data_val[key][i] + 1
            for hab in interactions["Halogen_Bonds"]:
                key = "Halogen_Bond:" + hab
                if key in data_val:
                    data_val[key][i] = data_val[key][i] + 1
            for wb in interactions["Water_Bridges"]:
                key = "Water_Bridge:" + wb
                if key in data_val:
                    data_val[key][i] = data_val[key][i] + 1
            for mc in interactions["Metal_Complexes"]:
                key = "Metal_Complex:" + mc
                if key in data_val:
                    data_val[key][i] = data_val[key][i] + 1

        df_val = pd.DataFrame(data_val)

        df_val["LABEL"] = df_val["NAME"].apply(lambda x: get_label(x))

        if df_val.loc[df_val["NAME"] == "dummy_name"].shape[0] == 0:
            df_val.to_csv(val_output, index = False)
        else:
            warnings.warn("Found dummy variables in dataset! Not saving to " + val_output, UserWarning)

        # -- generate test dataset --
        data_test = {
            "INDEX": list(range(1, len(entries_test) + 1)),
            "NAME": ["dummy_name" for entry in entries_test]
            }

        for i in interactions_train:
            data_test[i] = [0 for entry in entries_test]

        for i, entry in enumerate(entries_test):
            data_test["NAME"][i] = "|".join(entry["ligand_name"].split("|")[:-1])
            interactions = entry["interactions"]
            for hc in interactions["Hydrophobic_Contacts"]:
                key = "Hydrophobic_Interaction:" + hc
                if key in data_test:
                    data_test[key][i] = data_test[key][i] + 1
            for sb in interactions["Salt_Bridges"]:
                key = "Salt_Bridge:" + sb
                if key in data_test:
                    data_test[key][i] = data_test[key][i] + 1
            for hb in interactions["Hydrogen_Bonds"]:
                key = "Hydrogen_Bond:" + hb
                if key in data_test:
                    data_test[key][i] = data_test[key][i] + 1
            for ps in interactions["Pi_Stacking"]:
                key = "Pi-Stacking:" + ps
                if key in data_test:
                    data_test[key][i] = data_test[key][i] + 1
            for pc in interactions["Pi_Cation_Interactions"]:
                key = "Pi-Cation_Interaction:" + pc
                if key in data_test:
                    data_test[key][i] = data_test[key][i] + 1
            for hab in interactions["Halogen_Bonds"]:
                key = "Halogen_Bond:" + hab
                if key in data_test:
                    data_test[key][i] = data_test[key][i] + 1
            for wb in interactions["Water_Bridges"]:
                key = "Water_Bridge:" + wb
                if key in data_test:
                    data_test[key][i] = data_test[key][i] + 1
            for mc in interactions["Metal_Complexes"]:
                key = "Metal_Complex:" + mc
                if key in data_test:
                    data_test[key][i] = data_test[key][i] + 1

        df_test = pd.DataFrame(data_test)

        df_test["LABEL"] = df_test["NAME"].apply(lambda x: get_label(x))

        if df_test.loc[df_test["NAME"] == "dummy_name"].shape[0] == 0:
            df_test.to_csv(test_output, index = False)
        else:
            warnings.warn("Found dummy variables in dataset! Not saving to " + test_output, UserWarning)

        # return dataframes
        return [df_train, df_val, df_test]

    # get active and inactive structures from the training dataset
    def get_actives_inactives(self,
                              partition = "train"):

        """
        -- DESCRIPTION --
        """

        active_structures = []
        active_names = []
        inactive_structures = []
        inactive_names = []

        if partition == "val":
            entries_keys = self.entries_val_keys
        elif partition == "test":
            entries_keys = self.entries_test_keys
        else:
            entries_keys = self.entries_train_keys

        if entries_keys is not None:
            for key in entries_keys:
                if "inactive" in self.pdb_entry_results[key]["ligand_name"] or "decoy" in self.pdb_entry_results[key]["ligand_name"]:
                    inactive_structures.append(key)
                    inactive_names.append(self.pdb_entry_results[key]["ligand_name"])
                else:
                    active_structures.append(key)
                    active_names.append(self.pdb_entry_results[key]["ligand_name"])
            return {"active_structures": active_structures,
                    "active_names": active_names,
                    "inactive_structures": inactive_structures,
                    "inactive_names": inactive_names}
        else:
            warnings.warn("You need to generate the datasets first!", UserWarning)
            return 1

    # compare actives and inactives
    def compare(self,
                partition = "train",
                path = "current",
                **kwargs):

        """
        -- DESCRIPTION --
        """

        partition_info = self.get_actives_inactives(partition = partition)

        # get result for active ligands
        result_actives = PLIPAnalyzer(partition_info["active_structures"],
                                      ligand_names = partition_info["active_names"],
                                      path = path, **kwargs)

        # get result for inactive ligands
        result_inactives = PLIPAnalyzer(partition_info["inactive_structures"],
                                      ligand_names = partition_info["inactive_names"],
                                      path = path, **kwargs)

        # get comparison
        comparison = Comparison("Actives", "Inactives",
                                result_actives.i_frequencies,
                                result_inactives.i_frequencies)

        # return comparison
        return comparison

    # get information about interactions (active freq, inactive freq, diff)
    def get_feature_information(self,
                                filename = None,
                                **kwargs):

        """
        -- DESCRIPTION --
        """

        comparison = self.compare(partition = "train", **kwargs).comparison

        INTERACTIONS = []
        DIFFERENCES = []
        ACTIVE_FREQ = []
        INACTIVE_FREQ = []

        for key in comparison.keys():
            INTERACTIONS.append(key)
            DIFFERENCES.append(comparison[key]["difference"])
            ACTIVE_FREQ.append(comparison[key]["Actives"])
            INACTIVE_FREQ.append(comparison[key]["Inactives"])

        features = pd.DataFrame({"INDEX": list(range(1, len(INTERACTIONS) + 1)),
                                 "INTERACTION": INTERACTIONS,
                                 "DIFFERENCE": DIFFERENCES,
                                 "ACTIVE_FREQUENCY": ACTIVE_FREQ,
                                 "INACTIVE_FREQUENCY": INACTIVE_FREQ})

        if filename is not None:
            features.to_csv(filename, index = False)

        return features


class PLIPAnalyzer:
    """
    -- DESCRIPTION --
    """

    nr_structures = None
    normalized = None
    i_frequencies = None
    i_structures = None
    result = None

    # additional attributes
    pdb_entry_results = None
    best_pdb_entries = None

    # constructor with plip analysis for docking sdf files with multiple poses
    def __init__(self,
                 list_of_pdb_entries,
                 ligand_names = None,
                 poses = "all",
                 path = "current",
                 chain = "A",
                 exclude = ["LIG", "HOH"],
                 excluded_ligands = exclusion_list,
                 discard_exceeding_hc = True,
                 normalize = True,
                 verbose = 1):

        """
        -- DESCRIPTION --
        """

        INTERACTION_FREQ = {}
        INTERACTION_STRUC = {}
        RESULT = {}

        # if ligand names not given, generate unique ones for every structure
        if ligand_names is None:
            ligand_names = ["LIG" + str(i) + "|ligand|sdf|1|dock1" for i in range(1, len(list_of_pdb_entries) + 1)]

        # get nr of unique ligands
        unique_ligand_names = []
        for ligand_name in ligand_names:
            unique_ligand_names.append("|".join(ligand_name.split("|")[:-1]))

        # set attributes
        self.nr_structures = len(set(unique_ligand_names))
        self.normalized = normalize

        # constant for absolute / normalized frequencies
        if normalize:
            c = 1 / self.nr_structures
        else:
            c = 1

        # dictionary for all analyzed structures
        pdb_entry_results = {}

        for i, pdb_entry in enumerate(list_of_pdb_entries):

            # interactions for this entry
            e_Salt_Bridges = []
            e_Hydrogen_Bonds = []
            e_Pi_Stacking = []
            e_Pi_Cation_Interactions = []
            e_Hydrophobic_Contacts = []
            e_Halogen_Bonds = []
            e_Water_Bridges = []
            e_Metal_Complexes = []

            # load pdb file
            if path not in ["", ".", "current"]:
                filename = path + "/" + pdb_entry
            else:
                filename = pdb_entry
            if filename.split(".")[-1] != "pdb":
                filename = filename + ".pdb"
            mol = PDBComplex()
            if verbose:
                print("Analyzing file: ", filename)
            mol.load_pdb(filename)

            # get interactions
            # --> according to plipcmd.process_pdb()
            for ligand in mol.ligands:
                mol.characterize_complex(ligand)

            # iterate over interaction sets
            for key in mol.interaction_sets:
                iHet_ID, iChain, iPosition = key.split(":")
                # discard suspicious ligands
                # e.g. see -> https://github.com/pharmai/plip/blob/master/plip/basic/config.py
                # for complete exclusion list --> see top of file <COFACTORS>
                # this is expandable via the 'excluded_ligands' parameter
                if iHet_ID.strip().upper() in excluded_ligands:
                    continue
                # discard uninteressted chains
                if iChain != chain:
                    continue
                interaction = mol.interaction_sets[key]

                # get interaction residues
                # SALT BRIDGES
                tmp_salt_bridges = interaction.saltbridge_lneg + interaction.saltbridge_pneg
                Salt_Bridges = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in tmp_salt_bridges if i.restype not in exclude]
                e_Salt_Bridges = e_Salt_Bridges + Salt_Bridges
                # HYDROGEN BONDS
                tmp_h_bonds = interaction.hbonds_pdon + interaction.hbonds_ldon
                Hydrogen_Bonds = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in tmp_h_bonds if i.restype not in exclude]
                e_Hydrogen_Bonds = e_Hydrogen_Bonds + Hydrogen_Bonds
                # PI STACKING
                Pi_Stacking = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.pistacking if i.restype not in exclude]
                e_Pi_Stacking = e_Pi_Stacking + Pi_Stacking
                # PI CATION INTERACTION
                tmp_pication = interaction.pication_laro + interaction.pication_paro
                Pi_Cation_Interactions = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in tmp_pication if i.restype not in exclude]
                e_Pi_Cation_Interactions = e_Pi_Cation_Interactions + Pi_Cation_Interactions
                # HYDROPHOBIC CONTACTS
                Hydrophobic_Contacts = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.hydrophobic_contacts if i.restype not in exclude]
                e_Hydrophobic_Contacts = e_Hydrophobic_Contacts + Hydrophobic_Contacts
                # HALOGEN BONDS
                Halogen_Bonds = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.halogen_bonds if i.restype not in exclude]
                e_Halogen_Bonds = e_Halogen_Bonds + Halogen_Bonds
                # WATER BRIDGES
                Water_Bridges = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.water_bridges if i.restype not in exclude]
                e_Water_Bridges = e_Water_Bridges + Water_Bridges
                # METAL COMPLEXES
                Metal_Complexes = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.metal_complexes if i.restype not in exclude]
                e_Metal_Complexes = e_Metal_Complexes + Metal_Complexes

            # add analyzed structure to results via pdb_entry name
            # discard "copies" of hydrophobic contacts or not
            if discard_exceeding_hc:
                nr_of_interactions = len(e_Salt_Bridges) + len(e_Hydrogen_Bonds) + len(e_Pi_Stacking) + len(e_Pi_Cation_Interactions) + len(set(e_Hydrophobic_Contacts)) + len(e_Halogen_Bonds) + len(e_Water_Bridges) + len(e_Metal_Complexes)
                pdb_entry_results[pdb_entry] = {"ligand_name": ligand_names[i],
                                                "nr_of_interactions": nr_of_interactions,
                                                "interactions": {"Salt_Bridges": e_Salt_Bridges,
                                                                 "Hydrogen_Bonds": e_Hydrogen_Bonds,
                                                                 "Pi_Stacking": e_Pi_Stacking,
                                                                 "Pi_Cation_Interactions": e_Pi_Cation_Interactions,
                                                                 "Hydrophobic_Contacts": list(set(e_Hydrophobic_Contacts)),
                                                                 "Halogen_Bonds": e_Halogen_Bonds,
                                                                 "Water_Bridges": e_Water_Bridges,
                                                                 "Metal_Complexes": e_Metal_Complexes}
                                               }
            else:
                nr_of_interactions = len(e_Salt_Bridges) + len(e_Hydrogen_Bonds) + len(e_Pi_Stacking) + len(e_Pi_Cation_Interactions) + len(e_Hydrophobic_Contacts) + len(e_Halogen_Bonds) + len(e_Water_Bridges) + len(e_Metal_Complexes)
                pdb_entry_results[pdb_entry] = {"ligand_name": ligand_names[i],
                                                "nr_of_interactions": nr_of_interactions,
                                                "interactions": {"Salt_Bridges": e_Salt_Bridges,
                                                                 "Hydrogen_Bonds": e_Hydrogen_Bonds,
                                                                 "Pi_Stacking": e_Pi_Stacking,
                                                                 "Pi_Cation_Interactions": e_Pi_Cation_Interactions,
                                                                 "Hydrophobic_Contacts": e_Hydrophobic_Contacts,
                                                                 "Halogen_Bonds": e_Halogen_Bonds,
                                                                 "Water_Bridges": e_Water_Bridges,
                                                                 "Metal_Complexes": e_Metal_Complexes}
                                               }

        # get best pose for each unique ligand -- based on GOLD sdf naming schema!
        pdb_entry_results_sorted = dict(sorted(pdb_entry_results.items(), key = lambda x: x[1]["nr_of_interactions"], reverse = True))
        best_pdb_entries_dict = {}
        for key in pdb_entry_results_sorted.keys():
            best_pdb_entries_dict_key = "|".join(pdb_entry_results_sorted[key]["ligand_name"].split("|")[:-1])
            if best_pdb_entries_dict_key not in best_pdb_entries_dict:
                best_pdb_entries_dict[best_pdb_entries_dict_key] = key
        best_pdb_entries = list(best_pdb_entries_dict.values())

        # select structures for best / all poses
        if poses == "best":
            entries_to_be_used = best_pdb_entries
        else:
            entries_to_be_used = list_of_pdb_entries

        for pdb_entry in entries_to_be_used:

            # get interactions from dictonary
            Salt_Bridges = pdb_entry_results[pdb_entry]["interactions"]["Salt_Bridges"]
            Hydrogen_Bonds = pdb_entry_results[pdb_entry]["interactions"]["Hydrogen_Bonds"]
            Pi_Stacking = pdb_entry_results[pdb_entry]["interactions"]["Pi_Stacking"]
            Pi_Cation_Interactions = pdb_entry_results[pdb_entry]["interactions"]["Pi_Cation_Interactions"]
            Hydrophobic_Contacts = pdb_entry_results[pdb_entry]["interactions"]["Hydrophobic_Contacts"]
            Halogen_Bonds = pdb_entry_results[pdb_entry]["interactions"]["Halogen_Bonds"]
            Water_Bridges = pdb_entry_results[pdb_entry]["interactions"]["Water_Bridges"]
            Metal_Complexes = pdb_entry_results[pdb_entry]["interactions"]["Metal_Complexes"]

            # build dictonary
            for sb_residue in Salt_Bridges:
                k = "Salt_Bridge:" + sb_residue
                if k in INTERACTION_FREQ:
                    INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + c
                else:
                    INTERACTION_FREQ[k] = c
                if k in INTERACTION_STRUC:
                    INTERACTION_STRUC[k].append(pdb_entry)
                else:
                    INTERACTION_STRUC[k] = [pdb_entry]
                if k in RESULT:
                    RESULT[k]["frequency"] = RESULT[k]["frequency"] + c
                    RESULT[k]["structures"].append(pdb_entry)
                else:
                    RESULT[k] = {"frequency": c, "structures": [pdb_entry]}

            for hb_residue in Hydrogen_Bonds:
                k = "Hydrogen_Bond:" + hb_residue
                if k in INTERACTION_FREQ:
                    INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + c
                else:
                    INTERACTION_FREQ[k] = c
                if k in INTERACTION_STRUC:
                    INTERACTION_STRUC[k].append(pdb_entry)
                else:
                    INTERACTION_STRUC[k] = [pdb_entry]
                if k in RESULT:
                    RESULT[k]["frequency"] = RESULT[k]["frequency"] + c
                    RESULT[k]["structures"].append(pdb_entry)
                else:
                    RESULT[k] = {"frequency": c, "structures": [pdb_entry]}

            for ps_residue in Pi_Stacking:
                k = "Pi-Stacking:" + ps_residue
                if k in INTERACTION_FREQ:
                    INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + c
                else:
                    INTERACTION_FREQ[k] = c
                if k in INTERACTION_STRUC:
                    INTERACTION_STRUC[k].append(pdb_entry)
                else:
                    INTERACTION_STRUC[k] = [pdb_entry]
                if k in RESULT:
                    RESULT[k]["frequency"] = RESULT[k]["frequency"] + c
                    RESULT[k]["structures"].append(pdb_entry)
                else:
                    RESULT[k] = {"frequency": c, "structures": [pdb_entry]}

            for pc_residue in Pi_Cation_Interactions:
                k = "Pi-Cation_Interaction:" + pc_residue
                if k in INTERACTION_FREQ:
                    INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + c
                else:
                    INTERACTION_FREQ[k] = c
                if k in INTERACTION_STRUC:
                    INTERACTION_STRUC[k].append(pdb_entry)
                else:
                    INTERACTION_STRUC[k] = [pdb_entry]
                if k in RESULT:
                    RESULT[k]["frequency"] = RESULT[k]["frequency"] + c
                    RESULT[k]["structures"].append(pdb_entry)
                else:
                    RESULT[k] = {"frequency": c, "structures": [pdb_entry]}

            for hc_residue in Hydrophobic_Contacts:
                k = "Hydrophobic_Interaction:" + hc_residue
                if k in INTERACTION_FREQ:
                    INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + c
                else:
                    INTERACTION_FREQ[k] = c
                if k in INTERACTION_STRUC:
                    INTERACTION_STRUC[k].append(pdb_entry)
                else:
                    INTERACTION_STRUC[k] = [pdb_entry]
                if k in RESULT:
                    RESULT[k]["frequency"] = RESULT[k]["frequency"] + c
                    RESULT[k]["structures"].append(pdb_entry)
                else:
                    RESULT[k] = {"frequency": c, "structures": [pdb_entry]}

            for halogenb_residue in Halogen_Bonds:
                k = "Halogen_Bond:" + halogenb_residue
                if k in INTERACTION_FREQ:
                    INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + c
                else:
                    INTERACTION_FREQ[k] = c
                if k in INTERACTION_STRUC:
                    INTERACTION_STRUC[k].append(pdb_entry)
                else:
                    INTERACTION_STRUC[k] = [pdb_entry]
                if k in RESULT:
                    RESULT[k]["frequency"] = RESULT[k]["frequency"] + c
                    RESULT[k]["structures"].append(pdb_entry)
                else:
                    RESULT[k] = {"frequency": c, "structures": [pdb_entry]}

            for wb_residue in Water_Bridges:
                k = "Water_Bridge:" + wb_residue
                if k in INTERACTION_FREQ:
                    INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + c
                else:
                    INTERACTION_FREQ[k] = c
                if k in INTERACTION_STRUC:
                    INTERACTION_STRUC[k].append(pdb_entry)
                else:
                    INTERACTION_STRUC[k] = [pdb_entry]
                if k in RESULT:
                    RESULT[k]["frequency"] = RESULT[k]["frequency"] + c
                    RESULT[k]["structures"].append(pdb_entry)
                else:
                    RESULT[k] = {"frequency": c, "structures": [pdb_entry]}

            for mc_residue in Metal_Complexes:
                k = "Metal_Complexation:" + mc_residue
                if k in INTERACTION_FREQ:
                    INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + c
                else:
                    INTERACTION_FREQ[k] = c
                if k in INTERACTION_STRUC:
                    INTERACTION_STRUC[k].append(pdb_entry)
                else:
                    INTERACTION_STRUC[k] = [pdb_entry]
                if k in RESULT:
                    RESULT[k]["frequency"] = RESULT[k]["frequency"] + c
                    RESULT[k]["structures"].append(pdb_entry)
                else:
                    RESULT[k] = {"frequency": c, "structures": [pdb_entry]}

        # print warnings if no interactions are found or results do not agree
        if not INTERACTION_FREQ or not INTERACTION_STRUC or not RESULT:
            if not INTERACTION_FREQ and not INTERACTION_STRUC and not RESULT:
                warnings.warn("There are no found interactions!", NoInteractionsWarning)
            else:
                warnings.warn("It seems like there was an error during calculation! Results may not be correct or complete.", CalculationErrorWarning)

        self.i_frequencies = dict(sorted(INTERACTION_FREQ.items(), key = lambda x: x[1], reverse = True))
        self.i_structures = dict(sorted(INTERACTION_STRUC.items(), key = lambda x: len(x[1]), reverse = True))
        self.result = dict(sorted(RESULT.items(), key = lambda x: x[1]["frequency"], reverse = True))
        self.pdb_entry_results = pdb_entry_results
        self.best_pdb_entries = best_pdb_entries

    # save results as json with a given prefix
    def save(self, prefix, save_pdb_entry_results = False):

        """
        -- DESCRIPTION --
        """

        filename_1 = prefix + "_result.json"
        with open(filename_1, "w") as f:
            json.dump(self.result, f)
            f.close()

        filename_2 = prefix + "_frequencies.json"
        with open(filename_2, "w") as f:
            json.dump(self.i_frequencies, f)
            f.close()

        filename_3 = prefix + "_structures.json"
        with open(filename_3, "w") as f:
            json.dump(self.i_structures, f)
            f.close()

        if save_pdb_entry_results:
            filename_4 = prefix + "_pdb_entry_results.json"
            with open(filename_4, "w") as f:
                json.dump(self.pdb_entry_results, f)
                f.close()
            return [filename_1, filename_2, filename_3, filename_4]
        else:
            return [filename_1, filename_2, filename_3]

    # save frequencies as csv
    def to_csv(self, filename):

        """
        -- DESCRIPTION --
        """

        frequencies_csv = "Interaction;Frequency\n"
        for key in self.i_frequencies:
            frequencies_csv = frequencies_csv+ str(key) + ";" + str(self.i_frequencies[key]) + "\n"

            with open(filename, "w") as f:
                f.write(frequencies_csv)
                f.close()

        return frequencies_csv

    # plot frequencies with matplotlib
    def plot(self,
             title,
             filename = None,
             width = 20,
             height = 5,
             sig_digits = 4,
             label_offset = None):

        """
        -- DESCRIPTION --
        """

        # set label offset if not specified
        if label_offset is None:
            if self.normalized:
                label_offset = 1 / self.nr_structures
            else:
                label_offset = 1

        fig = plt.figure(figsize = (width, height))
        ax = fig.add_axes([0,0,1,1])
        plt.xticks(rotation = "vertical")
        ax.bar(self.i_frequencies.keys(), self.i_frequencies.values())
        xlocs, xlabs = plt.xticks()
        for i, v in enumerate(self.i_frequencies.values()):
            plt.text(xlocs[i], v + label_offset, str(round(v, sig_digits)), horizontalalignment = "center", rotation = 90)
        plt.title(title)
        plt.xlabel("Interaction")
        if self.normalized:
            plt.ylabel("Relative Frequency")
        else:
            plt.ylabel("Absolute Frequency")
        if filename is not None:
            fig.savefig(filename, bbox_inches = "tight", dpi = 150)
        plt.show()

        return fig
