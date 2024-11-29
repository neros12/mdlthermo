import json
import random
from typing import Tuple, List, Optional
from os.path import join as opj
from pathlib import Path


import pandas as pd


ROOT_DIR = Path(__file__).resolve().parents[2]


def parse_data(dir: Optional[str] = None) -> Tuple[
    List[str],
    List[str],
    List[float],
    List[float],
    List[float],
    List[float],
]:
    if dir is None:
        dir = opj(ROOT_DIR, "data", "DDB_Binary_VLE06.csv")

    df = pd.read_csv(dir)
    list_SMILES1 = df["cmp1_SMILES"].tolist()
    list_SMILES2 = df["cmp2_SMILES"].tolist()
    list_T = df["dTexp"].tolist()
    list_P = df["dPexp"].tolist()
    list_x1 = df["dX1exp"].tolist()
    list_y1 = df["dY1exp"].tolist()

    return list_SMILES1, list_SMILES2, list_T, list_P, list_x1, list_y1


def shuffle_data(list_SMILES1, list_SMILES2, list_T, list_P, list_x1, list_y1):
    indices = list(range(len(list_SMILES1)))
    random.shuffle(indices)

    shuffled_SMILES1 = [list_SMILES1[i] for i in indices]
    shuffled_SMILES2 = [list_SMILES2[i] for i in indices]
    shuffled_T = [list_T[i] for i in indices]
    shuffled_P = [list_P[i] for i in indices]
    shuffled_x1 = [list_x1[i] for i in indices]
    shuffled_y1 = [list_y1[i] for i in indices]

    return (
        shuffled_SMILES1,
        shuffled_SMILES2,
        shuffled_T,
        shuffled_P,
        shuffled_x1,
        shuffled_y1,
    )


def convert_data(
    list_SMILES1: List[str],
    list_SMILES2: List[str],
    list_T: List[float],
    list_P: List[float],
    list_x1: List[float],
    list_y1: List[float],
) -> Tuple[list, list]:

    with open(opj(ROOT_DIR, "data", "SMILES_to_Input.json"), "r") as file:
        SMILES_to_Input = json.load(file)

    nfm1 = [SMILES_to_Input[smiles]["nfm"] for smiles in list_SMILES1]
    efm1 = [SMILES_to_Input[smiles]["efm"] for smiles in list_SMILES1]
    nfm2 = [SMILES_to_Input[smiles]["nfm"] for smiles in list_SMILES2]
    efm2 = [SMILES_to_Input[smiles]["efm"] for smiles in list_SMILES2]

    x_value = [
        nfm1,
        efm1,
        nfm2,
        efm2,
        list_T,
        list_P,
        list_x1,
    ]
    y_value = list_y1

    return x_value, y_value


def retrieve_data():
    a, b, c, d, e, f = parse_data()
    # a, b, c, d, e, f = shuffle_data(a, b, c, d, e, f)
    x_value, y_value = convert_data(a, b, c, d, e, f)

    return x_value, y_value
