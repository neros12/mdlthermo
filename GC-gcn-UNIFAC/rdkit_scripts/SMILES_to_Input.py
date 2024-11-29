import json
from os.path import join as opj
from pathlib import Path


import pandas as pd


import utils


ROOT_DIR = Path(__file__).parent.parent


def get_SMILES_to_Input():
    df = pd.read_csv(opj(ROOT_DIR, "data", "DDB_Binary_VLE05.csv"))
    list_SMILES1 = df["cmp1_SMILES"].tolist()
    list_SMILES2 = df["cmp2_SMILES"].tolist()
    list_SMILES = list(set([*list_SMILES1, *list_SMILES2]))

    dict_SMILES = {}
    i = 0
    for SMILES in list_SMILES:
        i += 1
        print(i, len(list_SMILES))

        try:
            nfm, efm = utils.get_feature_matrces(SMILES)
            dict_SMILES[SMILES] = {"nfm": nfm.tolist(), "efm": efm.tolist()}
        except:
            pass

    with open(opj(ROOT_DIR, "data", "SMILES_to_Input.json"), "w") as f:
        json.dump(dict_SMILES, f, indent=4)


# with open(opj(ROOT_DIR, "data", "SMILES_to_Input.json"), "r") as file:
#     data_dict = json.load(file)

# df = pd.read_csv(opj(ROOT_DIR, "data", "DDB_Binary_VLE05.csv"))

# df = df[df["cmp1_SMILES"].isin(data_dict.keys())].reset_index(drop=True)
# df = df[df["cmp2_SMILES"].isin(data_dict.keys())].reset_index(drop=True)

# df.to_csv(opj(ROOT_DIR, "data", "DDB_Binary_VLE06.csv"), index=False)
# df.to_excel(opj(ROOT_DIR, "data", "DDB_Binary_VLE06.xlsx"), index=False)
