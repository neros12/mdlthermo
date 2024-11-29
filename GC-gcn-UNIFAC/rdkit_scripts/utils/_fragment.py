from typing import Tuple, List

import numpy as np
from scipy.linalg import fractional_matrix_power
from rdkit.Chem import Descriptors, rdchem, MolFromSmarts, MolFromSmiles, MolToSmiles


SMARTS = [
    "[CX4H3]",  # -CH3
    "[CX3H2v4]",  # =CH2
    "[CX2H1v4]",  # ≡CH
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
]


def fragment_mol(mol: rdchem.Mol) -> List[Tuple[Tuple[int], int]]:

    # H2O
    if MolToSmiles(mol) == "O" or MolToSmiles(mol) == "[2H]O[2H]":

        return [((0,), 49)]

    fragments = []
    heavy_atoms = np.arange(0, Descriptors.HeavyAtomCount(mol), 1)

    # -NO2
    if mol.HasSubstructMatch(MolFromSmarts("[NX3v4;!R](=O)[OX1]")) == True:
        candidates = mol.GetSubstructMatches(MolFromSmarts("[NX3v4;!R](=O)[OX1]"))
        for candidate in candidates:
            if set(candidate) & set(heavy_atoms) == set(candidate):
                fragments.append((candidate, 46))
                heavy_atoms = np.setdiff1d(heavy_atoms, candidate)
    # -COOH
    if mol.HasSubstructMatch(MolFromSmarts("[CH0v4;!R](=O)[OX2H1]")) == True:
        candidates = mol.GetSubstructMatches(MolFromSmarts("[CH0v4;!R](=O)[OX2H1]"))
        for candidate in candidates:
            if set(candidate) & set(heavy_atoms) == set(candidate):
                fragments.append((candidate, 47))
                heavy_atoms = np.setdiff1d(heavy_atoms, candidate)
    # -COO-
    if mol.HasSubstructMatch(MolFromSmarts("[CH0v4;!R](=O)[OX2H0]")) == True:
        candidates = mol.GetSubstructMatches(MolFromSmarts("[CH0v4;!R](=O)[OX2H0]"))
        for candidate in candidates:
            if set(candidate) & set(heavy_atoms) == set(candidate):
                fragments.append((candidate, 48))
                heavy_atoms = np.setdiff1d(heavy_atoms, candidate)

    # Others
    for index, subgroups in enumerate(SMARTS):
        if mol.HasSubstructMatch(MolFromSmarts(subgroups)) == True:
            candidates = mol.GetSubstructMatches(MolFromSmarts(subgroups))
            for candidate in candidates:
                if set(candidate) & set(heavy_atoms) == set(candidate):
                    fragments.append((candidate, index))
                    heavy_atoms = np.setdiff1d(heavy_atoms, candidate)

    if len(heavy_atoms) != 0:
        raise ValueError("Unfragmented atom has been found!")

    return fragments


def get_feature_matrces(SMILES: str, max_num_group=25) -> Tuple[np.ndarray, np.ndarray]:
    mol = MolFromSmiles(SMILES)
    if mol == None:
        raise ValueError("Unsupported SMILES format!")

    fragments = fragment_mol(mol)

    atom_to_group_idx = {}
    for i in range(len(fragments)):
        for j in range(len(fragments[i][0])):
            atom_to_group_idx[fragments[i][0][j]] = i

    # Node Feature Matrix
    nfm = np.zeros((max_num_group, len(SMARTS) + 4))
    for i in range(len(fragments)):
        nfm[i, fragments[i][1]] = 1

    # Edge Feature Matrix
    efm = np.zeros((len(fragments), len(fragments)))
    for i in range(len(fragments)):
        atom_idxs = fragments[i][0]
        for atom_idx in atom_idxs:
            target_atom = mol.GetAtomWithIdx(atom_idx)
            for neigbor_atom in target_atom.GetNeighbors():
                neigbor_atom_idx = neigbor_atom.GetIdx()
                if neigbor_atom_idx not in atom_idxs:
                    efm[i, atom_to_group_idx[neigbor_atom_idx]] += 1
    np.fill_diagonal(efm, 1)
    degmat = np.diag(np.sum(efm, axis=1))
    degmat_m0p5 = fractional_matrix_power(degmat, -0.5)
    efm = np.matmul(np.matmul(degmat_m0p5, efm), degmat_m0p5)
    padding = max_num_group - len(fragments)
    efm = np.pad(efm, ((0, padding), (0, padding)), "constant", constant_values=0.0)

    return nfm, efm
