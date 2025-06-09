import os
import json

from rdkit import Chem


DIR = os.path.dirname(__file__)
SUBGROUPS = [
    "[CX4H3]",
    "[CX3H2v4]",
    "[CX2H1v4]",
    "[!R;CX4H2]",
    "[!R;CX4H]",
    "[!R;CX4H0]",
    "[!R;CX3H1v4]",
    "[!R;CX3H0v4]",
    "[!R;CX2H0;$([CX2H0](=*)(=*))]",
    "[!R;CX2H0;$([CX2H0](#*))]",
    "[R;CX4H2]",
    "[R;CX4H]",
    "[R;CX4H0]",
    "[R;CX3H1v4]",
    "[R;CX3H0v4]",
    "[R;CX2H0;$([CX2H0](=*)(=*))]",
    "[R;CX2H0;$([CX2H0](#*))]",
    "[R;cX4h2]",
    "[R;cX4h]",
    "[R;cX3h1v4]",
    "[R;cX3h0v4]",
    "[FX1H0]",
    "[ClX1H0]",
    "[BrX1H0]",
    "[IX1H0]",
    "[OX2H1]",
    "[!R;OX2H0]",
    "[R;OX2H0]",
    "[R;oX2h0]",
    "[OX1H0v2]",
    "[NX3H2v3]",
    "[NX3H1v3;!R]",
    "[NX3H1v3;R]",
    "[nX3h1v3;R]",
    "[NX3H0v3;!R]",
    "[NX3H0v3;R]",
    "[nX3h0v3;R]",
    "[NX2H0v3;!R]",
    "[NX2H0v3;R]",
    "[nX2h0v3;R]",
    "[NX1H0v3]",
    "[SX2H1v2]",
    "[SX2H0v2;!R]",
    "[SX2H0v2;R]",
    "[sX2h0v2;R]",
    "[SX1H0v2]",
    "[OX1H0v1]",
    "[NX3H0v4]",
]
with open(os.path.join(DIR, "InChIKey_to_index.json")) as json_file:
    InChIKey_to_index = json.load(json_file)


def have_COSMO(SMILES: str) -> bool:
    mol = Chem.MolFromSmiles(SMILES)
    InChIKey = Chem.MolToInchiKey(mol)
    if InChIKey in InChIKey_to_index:
        return True
    else:
        return False


def is_only_consist_with(SMILES: str):
    mol = Chem.MolFromSmiles(SMILES)
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if num_heavy_atoms < 3:

        return False
    if num_heavy_atoms > 25:

        return False

    for subgroup in SUBGROUPS:
        patt = Chem.MolFromSmarts(subgroup)
        matches = mol.GetSubstructMatches(patt)
        num_heavy_atoms -= len(matches)
        if num_heavy_atoms == 0:

            return True

    return False
