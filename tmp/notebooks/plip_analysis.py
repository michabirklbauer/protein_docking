#!/usr/bin/env python3

# PLIP ANALYZER
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from plip.structure.preparation import PDBComplex
from plip.basic.config import biolip_list

def analyze_structures(list_of_pdb_entries,
                       path = "structures",
                       chain = "A",
                       exclude = ["LIG", "HOH"]):

    INTERACTION_FREQ = {}
    INTERACTION_STRUC = {}
    RESULT = {}

    for pdb_entry in list_of_pdb_entries:

        # load pdb file
        filename = path + "/" + pdb_entry + ".pdb"
        mol = PDBComplex()
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

    return [RESULT, INTERACTION_FREQ, INTERACTION_STRUC]
