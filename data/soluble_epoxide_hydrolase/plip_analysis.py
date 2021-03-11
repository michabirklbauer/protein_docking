#!/usr/bin/env python3

# PLIP ANALYZER
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/t0xic-m/
# micha.birklbauer@gmail.com

from plip.structure.preparation import PDBComplex
from plip.basic.config import biolip_list

def analyze_structures(list_of_pdb_entries, path = "structures", chain = "A"):

    INTERACTION_LIST = {}

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
            Salt_Bridges = [''.join([str(i.resnr), str(i.reschain)]) for i in tmp_salt_bridges if i.restype not in ['LIG', 'HOH']]
            # HYDROGEN BONDS
            tmp_h_bonds = interaction.hbonds_pdon + interaction.hbonds_ldon
            Hydrogen_Bonds = [''.join([str(i.resnr), str(i.reschain)]) for i in tmp_h_bonds if i.restype not in ['LIG', 'HOH']]
            # PI STACKING
            Pi_Stacking = [''.join([str(i.resnr), str(i.reschain)]) for i in interaction.pistacking if i.restype not in ['LIG', 'HOH']]
            # PI CATION INTERACTION
            tmp_pication = interaction.pication_laro + interaction.pication_paro
            Pi_Cation_Interactions = [''.join([str(i.resnr), str(i.reschain)]) for i in tmp_pication if i.restype not in ['LIG', 'HOH']]
            # HYDROPHOBIC CONTACTS
            Hydrophobic_Contacts = [''.join([str(i.resnr), str(i.reschain)]) for i in interaction.hydrophobic_contacts if i.restype not in ['LIG', 'HOH']]
            # HALOGEN BONDS
            Halogen_Bonds = [''.join([str(i.resnr), str(i.reschain)]) for i in interaction.halogen_bonds if i.restype not in ['LIG', 'HOH']]
            # WATER BRIDGES
            Water_Bridges = [''.join([str(i.resnr), str(i.reschain)]) for i in interaction.water_bridges if i.restype not in ['LIG', 'HOH']]
            # METAL COMPLEXES
            Metal_Complexes = [''.join([str(i.resnr), str(i.reschain)]) for i in interaction.metal_complexes if i.restype not in ['LIG', 'HOH']]

            # build dictonary
            for sb_residue in Salt_Bridges:
                k = "Salt_Bridge:" + sb_residue
                if k in INTERACTION_LIST:
                    INTERACTION_LIST[k] = INTERACTION_LIST[k] + 1
                else:
                    INTERACTION_LIST[k] = 1

            for hb_residue in Hydrogen_Bonds:
                k = "Hydrogen_Bond:" + hb_residue
                if k in INTERACTION_LIST:
                    INTERACTION_LIST[k] = INTERACTION_LIST[k] + 1
                else:
                    INTERACTION_LIST[k] = 1

            for ps_residue in Pi_Stacking:
                k = "Pi-Stacking:" + ps_residue
                if k in INTERACTION_LIST:
                    INTERACTION_LIST[k] = INTERACTION_LIST[k] + 1
                else:
                    INTERACTION_LIST[k] = 1

            for pc_residue in Pi_Cation_Interactions:
                k = "Pi-Cation_Interaction:" + pc_residue
                if k in INTERACTION_LIST:
                    INTERACTION_LIST[k] = INTERACTION_LIST[k] + 1
                else:
                    INTERACTION_LIST[k] = 1

            for hc_residue in Hydrophobic_Contacts:
                k = "Hydrophobic_Interaction:" + hc_residue
                if k in INTERACTION_LIST:
                    INTERACTION_LIST[k] = INTERACTION_LIST[k] + 1
                else:
                    INTERACTION_LIST[k] = 1

            for halogenb_residue in Halogen_Bonds:
                k = "Halogen_Bond:" + halogenb_residue
                if k in INTERACTION_LIST:
                    INTERACTION_LIST[k] = INTERACTION_LIST[k] + 1
                else:
                    INTERACTION_LIST[k] = 1

            for wb_residue in Water_Bridges:
                k = "Water_Bridge:" + wb_residue
                if k in INTERACTION_LIST:
                    INTERACTION_LIST[k] = INTERACTION_LIST[k] + 1
                else:
                    INTERACTION_LIST[k] = 1

            for mc_residue in Metal_Complexes:
                k = "Metal_Complexation:" + mc_residue
                if k in INTERACTION_LIST:
                    INTERACTION_LIST[k] = INTERACTION_LIST[k] + 1
                else:
                    INTERACTION_LIST[k] = 1

    return INTERACTION_LIST
