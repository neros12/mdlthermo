import os
import warnings
import logging
from typing import Tuple
from os.path import join as opj

# 로깅 수준 설정
warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import numpy as np
import tensorflow as tf
from pathlib import Path
from keras import models
from scipy.linalg import fractional_matrix_power
from rdkit.Chem import Descriptors, MolFromSmarts, MolFromSmiles

tf.get_logger().setLevel(logging.ERROR)
DIR_PATH = Path(__file__).parent

# Load Models
"""
TODO EXEC를 사용하지 않고 직접 입력하기...
"""
for i in range(10):
    exec(
        f"""
HFORM_{i} = models.load_model(opj(DIR_PATH, "models","HFORM","Model {i} Max=1400.h5"), compile=False)
HFUS_{i} = models.load_model(opj(DIR_PATH, "models","HFUS","Model {i} Max=100.h5"), compile=False)
PC_{i} = models.load_model(opj(DIR_PATH, "models","PC","Model {i} Max=13000.h5"), compile=False)
TBN_{i} = models.load_model(opj(DIR_PATH, "models","TBN","Model {i} Max=1620.h5"), compile=False)
TC_{i} = models.load_model(opj(DIR_PATH, "models","TC","Model {i} Max=2000.h5"), compile=False)
TF_{i} = models.load_model(opj(DIR_PATH, "models","TF","Model {i} Max=550.h5"), compile=False)
TMN_{i} = models.load_model(opj(DIR_PATH, "models","TMN","Model {i} Max=720.h5"), compile=False)
VC_{i} = models.load_model(opj(DIR_PATH, "models","VC","Model {i} Max=1.2.h5"), compile=False)
        """
    )

SMARTS = [
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


def fragement_group(SMILES: str) -> dict:
    mol_SMILES = MolFromSmiles(SMILES)
    Group_Index = {}
    List_HA = np.arange(0, Descriptors.HeavyAtomCount(mol_SMILES), 1)
    # -NO2
    if mol_SMILES.HasSubstructMatch(MolFromSmarts("[NX3v4;!R](=O)[OX1]")) == True:
        Candi = mol_SMILES.GetSubstructMatches(MolFromSmarts("[NX3v4;!R](=O)[OX1]"))
        for a in range(len(Candi)):
            if set(Candi[a]) & set(List_HA) == set(Candi[a]):
                Group_Index[len(Group_Index)] = [Candi[a], 46]
                List_HA = np.setdiff1d(List_HA, Candi[a])
    # -COOH
    if mol_SMILES.HasSubstructMatch(MolFromSmarts("[CH0v4;!R](=O)[OX2H1]")) == True:
        Candi = mol_SMILES.GetSubstructMatches(MolFromSmarts("[CH0v4;!R](=O)[OX2H1]"))
        for a in range(len(Candi)):
            if set(Candi[a]) & set(List_HA) == set(Candi[a]):
                Group_Index[len(Group_Index)] = [Candi[a], 47]
                List_HA = np.setdiff1d(List_HA, Candi[a])
    # -COO-
    if mol_SMILES.HasSubstructMatch(MolFromSmarts("[CH0v4;!R](=O)[OX2H0]")) == True:
        Candi = mol_SMILES.GetSubstructMatches(MolFromSmarts("[CH0v4;!R](=O)[OX2H0]"))
        for a in range(len(Candi)):
            if set(Candi[a]) & set(List_HA) == set(Candi[a]):
                Group_Index[len(Group_Index)] = [Candi[a], 48]
                List_HA = np.setdiff1d(List_HA, Candi[a])
    # Others
    for index, Subgroups in enumerate(SMARTS):
        if mol_SMILES.HasSubstructMatch(MolFromSmarts(Subgroups)) == True:
            Candi = mol_SMILES.GetSubstructMatches(MolFromSmarts(Subgroups))
            for a in range(len(Candi)):
                if set(Candi[a]) & set(List_HA) == set(Candi[a]):
                    Group_Index[len(Group_Index)] = [Candi[a], index]
                    List_HA = np.setdiff1d(List_HA, Candi[a])

    return Group_Index


def convert_gcgcn_input(SMILES: str, max_atom=25) -> Tuple[np.ndarray, np.ndarray]:
    try:
        mol = MolFromSmiles(SMILES)
        if mol == None:
            raise ValueError

        group_index_dict = fragement_group(SMILES)

        atom_to_group = {}
        for i in group_index_dict:
            for j in group_index_dict[i][0]:
                atom_to_group[j] = i

        # Node Feature Matrix
        nfm = np.zeros((max_atom, len(SMARTS) + 3))
        for group_index in group_index_dict:
            nfm[group_index, group_index_dict[group_index][1]] = 1

        # Edge Feature Matrix
        efm = np.zeros((len(group_index_dict), len(group_index_dict)))
        for group_index in group_index_dict:
            atom_idxs = group_index_dict[group_index][0]
            for atom_idx in atom_idxs:
                target_atom = mol.GetAtomWithIdx(atom_idx)
                for neigbor_atom in target_atom.GetNeighbors():
                    neigbor_atom_idx = neigbor_atom.GetIdx()
                    if neigbor_atom_idx not in atom_idxs:
                        efm[group_index, atom_to_group[neigbor_atom_idx]] += 1

        np.fill_diagonal(efm, 1)
        degmat = np.diag(np.sum(efm, axis=1))
        degmat_m0p5 = fractional_matrix_power(degmat, -0.5)
        efm = np.matmul(np.matmul(degmat_m0p5, efm), degmat_m0p5)
        padding = max_atom - len(group_index_dict)
        efm = np.pad(efm, ((0, padding), (0, padding)), "constant", constant_values=0.0)

        nfm = nfm.reshape((1, max_atom, len(SMARTS) + 3))
        efm = efm.reshape((1, max_atom, max_atom))

        return nfm, efm

    except:
        raise ValueError("Wrong SMILES format")


def predict_HFORM(SMILES: str) -> Tuple[float, float]:
    nfm, efm = convert_gcgcn_input(SMILES)
    result = []
    for i in range(10):
        exec(f"result.append(HFORM_{i}.predict([nfm, efm])[0,0] * 1400)")
    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_HFUS(SMILES: str) -> Tuple[float, float]:
    nfm, efm = convert_gcgcn_input(SMILES)
    result = []
    for i in range(10):
        exec(f"result.append(HFUS_{i}.predict([nfm, efm])[0,0] * 100)")
    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_PC(SMILES: str) -> Tuple[float, float]:
    nfm, efm = convert_gcgcn_input(SMILES)
    result = []
    for i in range(10):
        exec(f"result.append(PC_{i}.predict([nfm, efm])[0,0] * 13000)")
    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_TC(SMILES: str) -> Tuple[float, float]:
    nfm, efm = convert_gcgcn_input(SMILES)
    result = []
    for i in range(10):
        exec(f"result.append(TC_{i}.predict([nfm, efm])[0,0] * 2000)")
    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_TBN(SMILES: str) -> Tuple[float, float]:
    nfm, efm = convert_gcgcn_input(SMILES)
    TBN = []
    for i in range(10):
        exec(f"TBN.append(TBN_{i}.predict([nfm, efm])[0,0] * 1620)")
    TBN_val = np.average(TBN)
    TBN_unc = np.std(TBN)

    return TBN_val, TBN_unc


def predict_TF(SMILES: str) -> Tuple[float, float]:
    nfm, efm = convert_gcgcn_input(SMILES)
    result = []
    for i in range(10):
        exec(f"result.append(TF_{i}.predict([nfm, efm])[0,0] * 550)")
    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_TMN(SMILES: str) -> Tuple[float, float]:
    nfm, efm = convert_gcgcn_input(SMILES)
    result = []
    for i in range(10):
        exec(f"result.append(TMN_{i}.predict([nfm, efm])[0,0] * 720)")
    val = np.average(result)
    unc = np.std(result)

    return val, unc


def predict_VC(SMILES: str) -> Tuple[float, float]:
    nfm, efm = convert_gcgcn_input(SMILES)
    result = []
    for i in range(10):
        exec(f"result.append(VC_{i}.predict([nfm, efm])[0,0] * 1.2)")
    val = np.average(result)
    unc = np.std(result)

    return val, unc
