import pandas as pd
from rdkit import Chem

# df = pd.read_excel("components_list.xlsx")
# list_of_SMILES = df["SMILES"].to_list()
# list_of_InChI = []
# list_of_InChIKey = []
# for SMILES in list_of_SMILES:
#     mol = Chem.rdmolfiles.MolFromSmiles(SMILES)
#     if mol != None:
#         InChI = Chem.inchi.MolToInchi(mol)
#         InChIKey = Chem.inchi.MolToInchiKey(mol)

#         list_of_InChI.append(InChI)
#         list_of_InChIKey.append(InChIKey)
#     else:
#         list_of_InChI.append("")
#         list_of_InChIKey.append("")

# df["InChI"] = list_of_InChI
# df["InChIKey"] = list_of_InChIKey

# df.to_excel("components_list.xlsx", index=False)


# import csv

# df = pd.read_excel("component_list.xlsx")
# df.to_csv("component_list.csv", index=False, quoting=csv.QUOTE_ALL)
