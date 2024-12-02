import json
from os.path import join as opj
from pathlib import Path
from typing import Tuple

import numpy as np
from rdkit.Chem import Descriptors, MolFromSmarts, MolFromSmiles
from scipy.linalg import fractional_matrix_power


FILE_DIR = Path(__file__).parent

SMARTS_list = [
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
]


def fragment_molecule(SMILES: str) -> dict:
    """
    Fragment the molecule into groups.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.

    Returns
    -------
    dict
        Dictionary of the groups.
    """
    mol = MolFromSmiles(SMILES)
    group_index = {}
    heavy_atoms = np.arange(0, Descriptors.HeavyAtomCount(mol), 1)

    if len(heavy_atoms) < 3:
        raise ValueError("Unsupported Molecule")

    # Find -NO2 group
    if mol.HasSubstructMatch(MolFromSmarts("[NX3v4;!R](=O)[OX1]")):
        Candi = mol.GetSubstructMatches(MolFromSmarts("[NX3v4;!R](=O)[OX1]"))

        for a in range(len(Candi)):
            if set(Candi[a]) & set(heavy_atoms) == set(Candi[a]):
                group_index[len(group_index)] = [Candi[a], 46]
                heavy_atoms = np.setdiff1d(heavy_atoms, Candi[a])

    # Find -COOH group
    if mol.HasSubstructMatch(MolFromSmarts("[CH0v4;!R](=O)[OX2H1]")):
        Candi = mol.GetSubstructMatches(MolFromSmarts("[CH0v4;!R](=O)[OX2H1]"))

        for a in range(len(Candi)):
            if set(Candi[a]) & set(heavy_atoms) == set(Candi[a]):
                group_index[len(group_index)] = [Candi[a], 47]
                heavy_atoms = np.setdiff1d(heavy_atoms, Candi[a])

    # Find -COO- group
    if mol.HasSubstructMatch(MolFromSmarts("[CH0v4;!R](=O)[OX2H0]")):
        Candi = mol.GetSubstructMatches(MolFromSmarts("[CH0v4;!R](=O)[OX2H0]"))

        for a in range(len(Candi)):
            if set(Candi[a]) & set(heavy_atoms) == set(Candi[a]):
                group_index[len(group_index)] = [Candi[a], 48]
                heavy_atoms = np.setdiff1d(heavy_atoms, Candi[a])

    # Find Other groups
    for index, Subgroups in enumerate(SMARTS_list):
        if mol.HasSubstructMatch(MolFromSmarts(Subgroups)):
            Candi = mol.GetSubstructMatches(MolFromSmarts(Subgroups))

            for a in range(len(Candi)):
                if set(Candi[a]) & set(heavy_atoms) == set(Candi[a]):
                    group_index[len(group_index)] = [Candi[a], index]
                    heavy_atoms = np.setdiff1d(heavy_atoms, Candi[a])

    if len(heavy_atoms) != 0:
        raise ValueError("Unsupported Molecule!")

    return group_index


def _get_input_matrices(SMILES: str, max_atom=25) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert the molecule to input matrices.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.
    max_atom : int, optional
        Maximum number of atoms in the molecule. The default is 25.

    Raises
    ------
    ValueError
        If the molecule is not valid.

    Returns
    -------
    nfm : numpy.ndarray
        Node feature matrix.
    efm : numpy.ndarray
        Edge feature matrix.
    """
    # Check if the molecule is valid
    try:
        mol = MolFromSmiles(SMILES)
    except ValueError:
        raise ValueError("Wrong SMILES format")

    if mol is None:
        raise ValueError("Molecule is not valid.")

    group_index_dict = fragment_molecule(SMILES)

    atom_to_group = {}
    for i in group_index_dict:
        for j in group_index_dict[i][0]:
            atom_to_group[j] = i

    # Get node feature matrix
    nfm = np.zeros((max_atom, len(SMARTS_list) + 3))
    for group_index in group_index_dict:
        nfm[group_index, group_index_dict[group_index][1]] = 1

    # Get edge feature matrix
    efm = np.zeros((len(group_index_dict), len(group_index_dict)))
    for group_index in group_index_dict:
        atom_idxs = group_index_dict[group_index][0]

        for atom_idx in atom_idxs:
            target_atom = mol.GetAtomWithIdx(atom_idx)

            for neigbor_atom in target_atom.GetNeighbors():
                neigbor_atom_idx = neigbor_atom.GetIdx()

                if neigbor_atom_idx not in atom_idxs:
                    efm[group_index, atom_to_group[neigbor_atom_idx]] += 1

    # Fill diagonal to 1 to connect the group to itself
    np.fill_diagonal(efm, 1)

    # Normalize the edge feature matrix
    degmat = np.diag(np.sum(efm, axis=1))
    degmat_m0p5 = fractional_matrix_power(degmat, -0.5)
    efm = np.matmul(np.matmul(degmat_m0p5, efm), degmat_m0p5)

    # Padding the edge feature matrix to the maximum atom 25
    padding = max_atom - len(group_index_dict)
    efm = np.pad(efm, ((0, padding), (0, padding)), "constant", constant_values=0.0)

    return nfm, efm


def GCGCN(nfm: np.ndarray, efm: np.ndarray, param_list: list) -> float:
    """Predict the property of the molecule.

    Parameters
    ----------
    nfm : numpy.ndarray
        Node feature matrix.
    efm : numpy.ndarray
        Edge feature matrix.
    param_list : list
        Parameter list.

    Returns
    -------
    float
        Predicted property.
    """
    # GC-GCN layer
    x = np.dot(efm, nfm)
    x = np.dot(x, param_list[0]) + param_list[1]
    x = np.where(x > 0, x, 0)  # ReLU

    # GC-GCN layer
    x = np.dot(efm, x)
    x = np.dot(x, param_list[2]) + param_list[3]
    x = np.where(x > 0, x, 0)  # ReLU

    # Node-wise summation
    x = x.reshape(-1)
    x = np.dot(x, param_list[4])

    # Dense layer
    x = np.dot(x, param_list[5]) + param_list[6]
    x = np.where(x > 0, x, 0)

    # Dense layer
    x = np.dot(x, param_list[7]) + param_list[8]
    x = np.where(x > 0, x, 0)

    # Dense layer
    x = np.dot(x, param_list[9]) + param_list[10]

    return x


def predict_HFORM(SMILES: str) -> Tuple[float, float]:
    """
    Predict the heat of formation of the molecule.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.

    Returns
    -------
    Tuple[float, float]
        Predicted heat of formation and its uncertainty. (kJ/mol)
    """
    nfm, efm = _get_input_matrices(SMILES)

    result = []
    for i in range(10):
        with open(opj(FILE_DIR, "parameters", f"HFORM{i}.json"), "rb") as file:
            param_list = json.load(file)

        result.append(5400 * GCGCN(nfm, efm, param_list) - 4000)

    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_HFUS(SMILES: str) -> Tuple[float, float]:
    """
    Predict the heat of fusion of the molecule.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.

    Returns
    -------
    Tuple[float, float]
        Predicted heat of fusion and its uncertainty. (kJ/mol)
    """
    nfm, efm = _get_input_matrices(SMILES)

    result = []
    for i in range(10):
        with open(opj(FILE_DIR, "parameters", f"HFUS{i}.json"), "rb") as file:
            param_list = json.load(file)

        result.append(100 * GCGCN(nfm, efm, param_list))

    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_PC(SMILES: str) -> Tuple[float, float]:
    """
    Predict the critical pressure of the molecule.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.

    Returns
    -------
    Tuple[float, float]
        Predicted critical pressure and its uncertainty. (kPa)
    """
    nfm, efm = _get_input_matrices(SMILES)

    result = []
    for i in range(10):
        with open(opj(FILE_DIR, "parameters", f"PC{i}.json"), "rb") as file:
            param_list = json.load(file)

        result.append(13000 * GCGCN(nfm, efm, param_list))

    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_TC(SMILES: str) -> Tuple[float, float]:
    """
    Predict the critical temperature of the molecule.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.

    Returns
    -------
    Tuple[float, float]
        Predicted critical temperature and its uncertainty. (K)
    """
    nfm, efm = _get_input_matrices(SMILES)

    result = []
    for i in range(10):
        with open(opj(FILE_DIR, "parameters", f"TC{i}.json"), "rb") as file:
            param_list = json.load(file)

        result.append(2000 * GCGCN(nfm, efm, param_list))

    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_TBN(SMILES: str) -> Tuple[float, float]:
    """
    Predict the boiling point of the molecule.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.

    Returns
    -------
    Tuple[float, float]
        Predicted boiling point and its uncertainty. (K)
    """
    nfm, efm = _get_input_matrices(SMILES)

    result = []
    for i in range(10):
        with open(opj(FILE_DIR, "parameters", f"TBN{i}.json"), "rb") as file:
            param_list = json.load(file)

        result.append(1620 * GCGCN(nfm, efm, param_list))

    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_TF(SMILES: str) -> Tuple[float, float]:
    """
    Predict the flash point of the molecule.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.

    Returns
    -------
    Tuple[float, float]
        Predicted flash point and its uncertainty. (K)
    """
    nfm, efm = _get_input_matrices(SMILES)

    result = []
    for i in range(10):
        with open(opj(FILE_DIR, "parameters", f"TF{i}.json"), "rb") as file:
            param_list = json.load(file)

        result.append(550 * GCGCN(nfm, efm, param_list))

    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_TMN(SMILES: str) -> Tuple[float, float]:
    """
    Predict the melting point of the molecule.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.

    Returns
    -------
    Tuple[float, float]
        Predicted melting point and its uncertainty. (K)
    """
    nfm, efm = _get_input_matrices(SMILES)

    result = []
    for i in range(10):
        with open(opj(FILE_DIR, "parameters", f"TMN{i}.json"), "rb") as file:
            param_list = json.load(file)

        result.append(720 * GCGCN(nfm, efm, param_list))

    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_VC(SMILES: str) -> Tuple[float, float]:
    """
    Predict the critical volume of the molecule.

    Parameters
    ----------
    SMILES : str
        SMILES string of the molecule.

    Returns
    -------
    Tuple[float, float]
        Predicted critical volume and its uncertainty. (L/mol)
    """
    nfm, efm = _get_input_matrices(SMILES)

    result = []
    for i in range(10):
        with open(opj(FILE_DIR, "parameters", f"VC{i}.json"), "rb") as file:
            param_list = json.load(file)

        result.append(1.2 * GCGCN(nfm, efm, param_list))

    val = np.average(result)
    unc = np.std(result)

    return val, unc
