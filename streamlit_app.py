#!/usr/bin/env python3

# PLIP ANALYZER - STREAMLIT WEBUI
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

"""
#####################################################
##                                                 ##
##               -- STREAMLIT APP --               ##
##                                                 ##
#####################################################
"""

import os
import shutil
import streamlit as st
import traceback as tb
from rdkit import Chem
from biopandas.pdb import PandasPdb
from matplotlib import pyplot as plt
from plip.structure.preparation import PDBComplex
from plip.basic.config import biolip_list

# see PLIPAnalyzer for reference
def remove_ligands(input_file, output_file):

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

    return

# see PLIPAnalyzer for reference
def get_ligands(sdf_file):
    return [ x for x in Chem.SDMolSupplier(sdf_file)]

# see PLIPAnalyzer for reference
def add_ligands(input_file, output_file, ligands, ligand_IDs = None, ligand_chain_IDs = None, verbose = 1):

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

# main streamlit application
def main():

    st.title("Interaction Frequency Analyzer")

    st.write("""
    Calculate interaction frequencies accross different protein-ligand complexes based on PLIP.
    """)

    st.write("""
    ### Upload a PDB base structure!
    """)

    # PDB file uploader
    pdb_file = st.file_uploader("Upload a PDB File:", type=["pdb"])

    st.write("""
    ### Upload ligands in the form of a SDF file!
    """)

    # SDF file uploader
    sdf_file = st.file_uploader("Upload a SDF File:", type=["sdf"])

    analyze = st.button("Analyze!")

    progress_bar = st.progress(0.0)

    if analyze:
        if pdb_file is not None and sdf_file is not None:

            # make temporary directory for structures
            os.mkdir("streamlit_structures")

            # write uploaded files to tmp directory
            with open("streamlit_structures/pdb_file.pdb", "wb") as f1:
                f1.write(pdb_file.getbuffer())
            with open("streamlit_structures/sdf_file.sdf", "wb") as f2:
                f2.write(sdf_file.getbuffer())

            # remove ligands from base PDB structure
            remove_ligands("streamlit_structures/pdb_file.pdb", "streamlit_structures/pdb_file_cleaned.pdb")
            # get ligands from SDF file
            ligands = get_ligands("streamlit_structures/sdf_file.sdf")
            # add ligands to PDB structure - one ligand per file currently
            with st.spinner("Step 1/2: Adding ligands to structure..."):
                structures = []
                ligands_per_file = 1
                for i in range(0, len(ligands), ligands_per_file):
                    l_idx = i
                    r_idx = i + ligands_per_file
                    current_ligands = ligands[l_idx:r_idx]
                    output_file = "streamlit_structures/pdb_file_cleaned" + str(int(i / ligands_per_file + 1)) + ".pdb"
                    structures.append(output_file)
                    add_ligands("streamlit_structures/pdb_file_cleaned.pdb", output_file, current_ligands)
                    progress = (i + 1) / (len(ligands) / ligands_per_file)
                    progress_bar.progress(progress)

            # PLIPAnalyzer-Step - see PLIPAnalyzer for reference
            INTERACTION_FREQ = {}
            chain = "A"
            exclude = ["LIG", "HOH"]

            with st.spinner("Step 2/2: Analyzing structures..."):
                progress_bar.progress(0.0)
                for i, pdb_entry in enumerate(structures):

                    # load pdb file
                    mol = PDBComplex()
                    mol.load_pdb(pdb_entry)

                    # get interactions
                    # --> according to plipcmd.process_pdb()
                    for ligand in mol.ligands:
                        mol.characterize_complex(ligand)

                    # iterate over interaction sets
                    for key in mol.interaction_sets:
                        iHet_ID, iChain, iPosition = key.split(":")
                        print(key)
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

                        for hb_residue in Hydrogen_Bonds:
                            k = "Hydrogen_Bond:" + hb_residue
                            if k in INTERACTION_FREQ:
                                INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                            else:
                                INTERACTION_FREQ[k] = 1

                        for ps_residue in Pi_Stacking:
                            k = "Pi-Stacking:" + ps_residue
                            if k in INTERACTION_FREQ:
                                INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                            else:
                                INTERACTION_FREQ[k] = 1

                        for pc_residue in Pi_Cation_Interactions:
                            k = "Pi-Cation_Interaction:" + pc_residue
                            if k in INTERACTION_FREQ:
                                INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                            else:
                                INTERACTION_FREQ[k] = 1

                        for hc_residue in Hydrophobic_Contacts:
                            k = "Hydrophobic_Interaction:" + hc_residue
                            if k in INTERACTION_FREQ:
                                INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                            else:
                                INTERACTION_FREQ[k] = 1

                        for halogenb_residue in Halogen_Bonds:
                            k = "Halogen_Bond:" + halogenb_residue
                            if k in INTERACTION_FREQ:
                                INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                            else:
                                INTERACTION_FREQ[k] = 1

                        for wb_residue in Water_Bridges:
                            k = "Water_Bridge:" + wb_residue
                            if k in INTERACTION_FREQ:
                                INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                            else:
                                INTERACTION_FREQ[k] = 1

                        for mc_residue in Metal_Complexes:
                            k = "Metal_Complexation:" + mc_residue
                            if k in INTERACTION_FREQ:
                                INTERACTION_FREQ[k] = INTERACTION_FREQ[k] + 1
                            else:
                                INTERACTION_FREQ[k] = 1

                    progress = (i + 1) / len(structures)
                    progress_bar.progress(progress)

            frequencies = dict(sorted(INTERACTION_FREQ.items(), key = lambda x: x[1], reverse = True))

            # remove tmp directory
            shutil.rmtree("streamlit_structures")

            # plot result
            fig = plt.figure(figsize = (20, 5))
            ax = fig.add_axes([0,0,1,1])
            plt.xticks(rotation = "vertical")
            ax.bar(frequencies.keys(), frequencies.values())
            xlocs, xlabs = plt.xticks()
            for i, v in enumerate(frequencies.values()):
                plt.text(xlocs[i], v + 1, str(v), horizontalalignment = "center")
            plt.title(str(sdf_file.name))
            plt.xlabel("Interaction")
            plt.ylabel("Absolute Frequency")

            st.pyplot(fig)

if __name__ == "__main__":

    main()
