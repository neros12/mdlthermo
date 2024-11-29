from os.path import join as opj
from pathlib import Path
import time
import numpy as np
import requests, json
import pandas as pd


def convert_to_number(rn: str) -> int:

    return int(rn.replace("-", ""))


df = pd.read_csv(opj(Path(__file__).parent.parent, "data", "DDB_Binary_VLE05.csv"))

l1 = df["cmp1_SMILES"].to_list()
l2 = df["cmp2_SMILES"].to_list()

SMILES_list = list(set((*l1, *l2)))
CASRN_list = []
_i = 0
for SMILES in SMILES_list:
    _i += 1
    for _j in range(5):
        time.sleep(np.random.rand(1)[0])
        try:
            request_result = requests.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{SMILES}/xrefs/RN/JSON"
            )
            json_result = json.loads(request_result.text)
            CASRNs = json_result["InformationList"]["Information"][0]["RN"]

            CASRN = min(CASRNs, key=convert_to_number)
            CASRN = convert_to_number(CASRN)
            print(SMILES, CASRN, _i, len(SMILES_list))
            CASRN_list.append(CASRN)

            break
        except:
            pass

        if _j == 4:
            print(SMILES, "NONE", _i, len(SMILES_list))
            CASRN_list.append("")


df = pd.DataFrame(data={"CASRN": CASRN_list, "SMILES": SMILES_list})
df.to_csv("CASRN_SMILES.csv", index=False)
df.to_excel("CASRN_SMILES.xlsx", index=False)
