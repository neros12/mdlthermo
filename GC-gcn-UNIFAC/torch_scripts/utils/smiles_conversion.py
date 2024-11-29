import json
import subprocess
from typing import Tuple, List
from os.path import join as opj
from pathlib import Path


import numpy as np


ROOT_DIR = Path(__file__).parent.parent.parent
RDKIT_ENV = opj(Path.home(), "anaconda3", "envs", "rdkit-env", "python.exe")


def get_feature_matrices(SMILES: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    단일 SMILES에 대한 특성 행렬을 반환

    Parameters
    ----------
    smiles (str): SMILES 문자열


    Returns
    -------
    Tuple[np.ndarray, np.ndarray]:
        - nfm (numpy.ndarray): 노드 특성 행렬 (25, 50)
        - efm (numpy.ndarray): 엣지 특성 행렬 (25, 25)
    """
    PYTHON_SCRIPT = opj(ROOT_DIR, "rdkit_scripts", "utils", "run_get_feature.py")

    run_result = subprocess.run(
        [RDKIT_ENV, PYTHON_SCRIPT, SMILES],
        capture_output=True,
        text=True,
    )
    json_result = json.loads(run_result.stdout)

    return np.array(json_result["nfm"]), np.array(json_result["efm"])


def get_batch_feature_matrices(list_SMILES: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    """
    다수의 SMILES에 대한 특성 행렬을 반환

    Parameters
    ----------
    list_smiles (List[str]): SMILES 문자열들의 리스트

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]:
        - nfms (numpy.ndarray): 여러 SMILES의 노드 특성 행렬 리스트 (*, 25, 50)
        - efms (numpy.ndarray): 여러 SMILES의 엣지 특성 행렬 리스트 (*, 25, 25)
    """
    nfms = []
    efms = []
    for SMILES in list_SMILES:
        nfm, efm = get_feature_matrices(SMILES)
        nfms.append(nfm)
        efms.append(efm)

    return np.array(nfms), np.array(efms)
