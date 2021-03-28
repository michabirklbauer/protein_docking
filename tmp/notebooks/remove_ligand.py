from biopandas.pdb import PandasPdb
from plip.structure.preparation import PDBComplex

def clean_pdb(filepath, outpath):

    mol = PDBComplex()
    mol.load_pdb(filepath)
    for ligand in mol.ligands:
        mol.characterize_complex(ligand)

    ligands = [str(lig.split(":")[0]).strip() for lig in mol.interaction_sets.keys()]

    pdb = PandasPdb().read_pdb(filepath)

    def remove_conect(residue_name, pdb):

        atoms = pdb.df["HETATM"][pdb.df["HETATM"]["residue_name"] == residue_name]["atom_number"]
        atoms_str = [str(atom) for atom in atoms]

        df = pdb.df["OTHERS"]
        idx = []
        for index, row in df.iterrows():
            if row["record_name"] == "CONECT":
                if any(a in row["entry"] for a in atoms_str):
                    idx.append(index)
        return idx

    def remove_hetatm(residue_name, pdb):
        df = pdb.df["HETATM"]
        idx = []
        for index, row in df.iterrows():
            if residue_name in str(row["residue_name"]):
                idx.append(index)
        return idx

    for ligand in ligands:
        pdb.df["OTHERS"].drop(remove_conect(ligand, pdb), inplace = True)
        pdb.df["HETATM"].drop(remove_hetatm(ligand, pdb), inplace = True)

    pdb.to_pdb(outpath)
