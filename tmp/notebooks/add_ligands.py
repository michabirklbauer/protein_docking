from rdkit import Chem
from biopandas.pdb import PandasPdb

def get_ligands(sdf_file):
    return [ x for x in Chem.SDMolSupplier(sdf_file)]

def add_ligands(pdb_file, outpath, ligands, ligand_chain_ids = None, verbose = 1):

    if ligand_chain_ids is None:
            ligand_chain_ids = ["A" for x in ligands]

    pdb = Chem.MolFromPDBFile(pdb_file)

    i = 1
    for ligand in ligands:
        merged = Chem.CombineMols(pdb, ligand)
        Chem.MolToPDBFile(merged, outpath)

        ligand_name = "LG" + str(i)
        if verbose:
            print("Adding ligand " + ligand_name + " to file " + outpath + "!")

        pdb = PandasPdb().read_pdb(outpath)
        pdb.df["HETATM"].loc[pdb.df["HETATM"]["residue_name"] == "UNL", "chain_id"] = ligand_chain_ids[i-1]
        pdb.df["HETATM"].replace("UNL", ligand_name, inplace = True)
        pdb.to_pdb(outpath)

        pdb = Chem.MolFromPDBFile(outpath)
        i = i + 1
