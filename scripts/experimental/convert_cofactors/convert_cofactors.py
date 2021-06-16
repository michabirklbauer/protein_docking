#!/usr/bin/env python3

# CONVERT COFACTORS
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

"""
DESCRIPTION
A function to re-label HETATM entries of a cofactor in a PDB file to ATOM
entries to mark cofactor as part of the protein. Entries are re-labelled
successfully but cofactor is still detected as ligand (unwanted behaviour).
"""

from biopandas.pdb import PandasPdb

def convert_cofactors(input_file, cofactor_shortcode, output_file):
    pdb = PandasPdb().read_pdb(input_file)
    pdb.df["HETATM"]["record_name"] = pdb.df["HETATM"]["residue_name"].apply(lambda x: "ATOM" if x == cofactor_shortcode else "HETATM")
    pdb.to_pdb(output_file)
    return pdb
