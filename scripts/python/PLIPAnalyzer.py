#!/usr/bin/env python3

# PLIP ANALYZER
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

version = "0.2.0"
date = "20210503"

import json
import warnings
import numpy as np
import pandas as pd
import traceback as tb
from rdkit import Chem
from biopandas.pdb import PandasPdb
from matplotlib import pyplot as plt
from plip.structure.preparation import PDBComplex
from plip.basic.config import biolip_list

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

class Preparation:
    """to be implemented"""

    # empty constructor
    def __init__(self):
        pass

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
    """to be implemented"""
    pass

class PLIPAnalyzer:
    """
    -- DESCRIPTION --
    """

    i_frequencies = None
    i_structures = None
    result = None

    # constructor with plip analysis
    def __init__(self,
                 list_of_pdb_entries,
                 path = "current",
                 chain = "A",
                 exclude = ["LIG", "HOH"],
                 verbose = 1):

        """
        -- DESCRIPTION --
        """

        INTERACTION_FREQ = {}
        INTERACTION_STRUC = {}
        RESULT = {}

        for pdb_entry in list_of_pdb_entries:

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
                if iHet_ID.strip().upper() in biolip_list:
                    continue
                # discard uninteressted chains
                if iChain != chain:
                    continue
                interaction = mol.interaction_sets[key]

                # get interaction residues
                # SALT BRIDGES
                tmp_salt_bridges = interaction.saltbridge_lneg + interaction.saltbridge_pneg
                Salt_Bridges = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in tmp_salt_bridges if i.restype not in exclude]
                # HYDROGEN BONDS
                tmp_h_bonds = interaction.hbonds_pdon + interaction.hbonds_ldon
                Hydrogen_Bonds = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in tmp_h_bonds if i.restype not in exclude]
                # PI STACKING
                Pi_Stacking = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.pistacking if i.restype not in exclude]
                # PI CATION INTERACTION
                tmp_pication = interaction.pication_laro + interaction.pication_paro
                Pi_Cation_Interactions = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in tmp_pication if i.restype not in exclude]
                # HYDROPHOBIC CONTACTS
                Hydrophobic_Contacts = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.hydrophobic_contacts if i.restype not in exclude]
                # HALOGEN BONDS
                Halogen_Bonds = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.halogen_bonds if i.restype not in exclude]
                # WATER BRIDGES
                Water_Bridges = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.water_bridges if i.restype not in exclude]
                # METAL COMPLEXES
                Metal_Complexes = [''.join([str(i.restype), str(i.resnr), str(i.reschain)]) for i in interaction.metal_complexes if i.restype not in exclude]

                # build dictonary
                for sb_residue in Salt_Bridges:
                    k = "Salt_Bridge:" + sb_residue
                    if k in INTERACTION_FREQ:
                        INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                    else:
                        INTERACTION_FREQ[k] = 1
                    if k in INTERACTION_STRUC:
                        INTERACTION_STRUC[k].append(pdb_entry)
                    else:
                        INTERACTION_STRUC[k] = [pdb_entry]
                    if k in RESULT:
                        RESULT[k]["frequency"] = RESULT[k]["frequency"] + 1
                        RESULT[k]["structures"].append(pdb_entry)
                    else:
                        RESULT[k] = {"frequency": 1, "structures": [pdb_entry]}

                for hb_residue in Hydrogen_Bonds:
                    k = "Hydrogen_Bond:" + hb_residue
                    if k in INTERACTION_FREQ:
                        INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                    else:
                        INTERACTION_FREQ[k] = 1
                    if k in INTERACTION_STRUC:
                        INTERACTION_STRUC[k].append(pdb_entry)
                    else:
                        INTERACTION_STRUC[k] = [pdb_entry]
                    if k in RESULT:
                        RESULT[k]["frequency"] = RESULT[k]["frequency"] + 1
                        RESULT[k]["structures"].append(pdb_entry)
                    else:
                        RESULT[k] = {"frequency": 1, "structures": [pdb_entry]}

                for ps_residue in Pi_Stacking:
                    k = "Pi-Stacking:" + ps_residue
                    if k in INTERACTION_FREQ:
                        INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                    else:
                        INTERACTION_FREQ[k] = 1
                    if k in INTERACTION_STRUC:
                        INTERACTION_STRUC[k].append(pdb_entry)
                    else:
                        INTERACTION_STRUC[k] = [pdb_entry]
                    if k in RESULT:
                        RESULT[k]["frequency"] = RESULT[k]["frequency"] + 1
                        RESULT[k]["structures"].append(pdb_entry)
                    else:
                        RESULT[k] = {"frequency": 1, "structures": [pdb_entry]}

                for pc_residue in Pi_Cation_Interactions:
                    k = "Pi-Cation_Interaction:" + pc_residue
                    if k in INTERACTION_FREQ:
                        INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                    else:
                        INTERACTION_FREQ[k] = 1
                    if k in INTERACTION_STRUC:
                        INTERACTION_STRUC[k].append(pdb_entry)
                    else:
                        INTERACTION_STRUC[k] = [pdb_entry]
                    if k in RESULT:
                        RESULT[k]["frequency"] = RESULT[k]["frequency"] + 1
                        RESULT[k]["structures"].append(pdb_entry)
                    else:
                        RESULT[k] = {"frequency": 1, "structures": [pdb_entry]}

                for hc_residue in Hydrophobic_Contacts:
                    k = "Hydrophobic_Interaction:" + hc_residue
                    if k in INTERACTION_FREQ:
                        INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                    else:
                        INTERACTION_FREQ[k] = 1
                    if k in INTERACTION_STRUC:
                        INTERACTION_STRUC[k].append(pdb_entry)
                    else:
                        INTERACTION_STRUC[k] = [pdb_entry]
                    if k in RESULT:
                        RESULT[k]["frequency"] = RESULT[k]["frequency"] + 1
                        RESULT[k]["structures"].append(pdb_entry)
                    else:
                        RESULT[k] = {"frequency": 1, "structures": [pdb_entry]}

                for halogenb_residue in Halogen_Bonds:
                    k = "Halogen_Bond:" + halogenb_residue
                    if k in INTERACTION_FREQ:
                        INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                    else:
                        INTERACTION_FREQ[k] = 1
                    if k in INTERACTION_STRUC:
                        INTERACTION_STRUC[k].append(pdb_entry)
                    else:
                        INTERACTION_STRUC[k] = [pdb_entry]
                    if k in RESULT:
                        RESULT[k]["frequency"] = RESULT[k]["frequency"] + 1
                        RESULT[k]["structures"].append(pdb_entry)
                    else:
                        RESULT[k] = {"frequency": 1, "structures": [pdb_entry]}

                for wb_residue in Water_Bridges:
                    k = "Water_Bridge:" + wb_residue
                    if k in INTERACTION_FREQ:
                        INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                    else:
                        INTERACTION_FREQ[k] = 1
                    if k in INTERACTION_STRUC:
                        INTERACTION_STRUC[k].append(pdb_entry)
                    else:
                        INTERACTION_STRUC[k] = [pdb_entry]
                    if k in RESULT:
                        RESULT[k]["frequency"] = RESULT[k]["frequency"] + 1
                        RESULT[k]["structures"].append(pdb_entry)
                    else:
                        RESULT[k] = {"frequency": 1, "structures": [pdb_entry]}

                for mc_residue in Metal_Complexes:
                    k = "Metal_Complexation:" + mc_residue
                    if k in INTERACTION_FREQ:
                        INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                    else:
                        INTERACTION_FREQ[k] = 1
                    if k in INTERACTION_STRUC:
                        INTERACTION_STRUC[k].append(pdb_entry)
                    else:
                        INTERACTION_STRUC[k] = [pdb_entry]
                    if k in RESULT:
                        RESULT[k]["frequency"] = RESULT[k]["frequency"] + 1
                        RESULT[k]["structures"].append(pdb_entry)
                    else:
                        RESULT[k] = {"frequency": 1, "structures": [pdb_entry]}

        # print warnings if no interactions are found or results do not agree
        if not INTERACTION_FREQ or not INTERACTION_STRUC or not RESULT:
            if not INTERACTION_FREQ and not INTERACTION_STRUC and not RESULT:
                warnings.warn("There are no found interactions!", NoInteractionsWarning)
            else:
                warnings.warn("It seems like there was an error during calculation! Results may not be correct or complete.", CalculationErrorWarning)

        self.i_frequencies = dict(sorted(INTERACTION_FREQ.items(), key = lambda x: x[1], reverse = True))
        self.i_structures = dict(sorted(INTERACTION_STRUC.items(), key = lambda x: len(x[1]), reverse = True))
        self.result = dict(sorted(RESULT.items(), key = lambda x: x[1]["frequency"], reverse = True))

    # save results as json with a given prefix
    def save(self, prefix):

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
             label_offset = 1):

        """
        -- DESCRIPTION --
        """

        fig = plt.figure(figsize = (width, height))
        ax = fig.add_axes([0,0,1,1])
        plt.xticks(rotation = "vertical")
        ax.bar(self.i_frequencies.keys(), self.i_frequencies.values())
        xlocs, xlabs = plt.xticks()
        for i, v in enumerate(self.i_frequencies.values()):
            plt.text(xlocs[i], v + label_offset, str(v), horizontalalignment = "center")
        plt.title(title)
        plt.xlabel("Interaction")
        plt.ylabel("Absolute Frequency")
        if filename is not None:
            fig.savefig(filename, bbox_inches = "tight", dpi = 150)
        plt.show()

        return fig
