# import pandas as pd
# import re

# # df = pd.read_excel("idCmp.xlsx", dtype=str)

# for row in df.itertuples():
#     name = row.name
#     if re.match(r"PLEASE", name):
#         target_idCmp = name.split()[4].split("_")[0]
#         print(f"CHECK : {target_idCmp} = {df.loc[int(target_idCmp) - 1, 'idCmp']}")
#         df.loc[row.Index] = [
#             row.idCmp,
#             df.loc[int(target_idCmp) - 1, "CASRN"],
#             df.loc[int(target_idCmp) - 1, "name"],
#             df.loc[int(target_idCmp) - 1, "formula"],
#             df.loc[int(target_idCmp) - 1, "SMILES"],
#         ]

# df.to_excel("idCmp02.xlsx", index=False)

# <<< END OF THE SCRIPT <<<

# import pandas as pd
# import requests
# import json


# def get_SMILES(name):
#     url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
#     request_results = requests.get(url)
#     request_results = json.loads(request_results.text)["PropertyTable"]["Properties"][0]

#     return request_results["IsomericSMILES"]


# df = pd.read_excel("idCmp02.xlsx", dtype=str)

# for row in df.itertuples():
#     SMILES = row.SMILES
#     name = row.name

#     if type(SMILES) == float and type(name) != float:
#         try:
#             SMILES_result = get_SMILES(name)
#             df.loc[row.Index, "SMILES"] = SMILES_result
#             print(row.name, SMILES_result)
#         except:
#             print(name)

# df.to_excel("idCmp03.xlsx", index=False)

# <<< END OF THE SCRIPT <<<

import pandas as pd
import requests
import json


def get_SMILES(casrn):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{casrn}/property/IsomericSMILES/JSON"
    request_results = requests.get(url)
    request_results = json.loads(request_results.text)["PropertyTable"]["Properties"][0]

    return request_results["IsomericSMILES"]


df = pd.read_excel("idCmp04.xlsx", dtype=str)

for row in df.itertuples():
    SMILES = row.SMILES
    CASRN = row.CASRN

    if type(SMILES) == float and type(CASRN) != float:
        try:
            SMILES_result = get_SMILES(CASRN)
            df.loc[row.Index, "SMILES"] = SMILES_result
            print(row.CASRN, SMILES_result)
        except:
            print(CASRN)

df.to_excel("idCmp05.xlsx", index=False)

# for row in df.itertuples():
#     SMILES = row.SMILES
#     if type(SMILES) != float:
#         mol = Chem.MolFromSmiles(SMILES)
#         formula = row.formula
#         rdkit_formula = rdMolDescriptors.CalcMolFormula(mol).upper()

#         if formula != rdkit_formula:
#             print(
#                 "idCmp: ",
#                 row.idCmp,
#                 "  formula :",
#                 row.formula,
#                 "  rdkit_formula :",
#                 rdkit_formula,
#             )
#             pass
