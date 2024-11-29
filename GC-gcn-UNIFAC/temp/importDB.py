import pandas as pd
import requests
import json


df_idCmp = pd.read_excel("idCmp05.xlsx", dtype=str)

MIXT01_DATA_VLE = pd.read_excel(
    "DDB2007Original.xlsx", sheet_name="HC_MIXT01_DATA_VLE", dtype=str
)
MIXT01_ID_VLE = pd.read_excel(
    "DDB2007Original.xlsx", sheet_name="HC_MIXT01_ID_VLE", dtype=str
)
# >>> REFERENCE DB >>>
# MIXT01_REF_VLE = pd.read_excel("DDB2007Original.xlsx", sheet_name="HC_MIXT01_REF_VLE")
# <<< REFERENCE DB <<<


df = pd.DataFrame(
    columns=[
        "cmp1_name",
        "cmp2_name",
        "cmp1_SMILES",
        "cmp2_SMILES",
        "dTexp",
        "dPexp",
        "dX1exp",
        "dY1exp",
    ]
)


dict_idCmp = {}
for row in df_idCmp.itertuples():
    idCmp = row.idCmp
    CASRN = row.CASRN
    name = row.name
    formula = row.formula
    SMILES = row.SMILES

    dict_idCmp[idCmp] = {
        "CASRN": CASRN,
        "name": name,
        "formula": formula,
        "SMILES": SMILES,
    }


dict_idExpSet = {}
for row in MIXT01_ID_VLE.itertuples():
    idExpSet = row.idExpSet
    nComp = row.nComp
    idCmp1 = row.idCmp1
    idCmp2 = row.idCmp2

    if int(nComp) == 2:
        dict_idExpSet[idExpSet] = {
            "idCmp1": idCmp1,
            "idCmp2": idCmp2,
        }

i = -1
for row in MIXT01_DATA_VLE.itertuples():
    i += 1
    print("Writing Data : ", i, len(MIXT01_DATA_VLE))

    idExpSet = row.idExpSet
    dTexp = row.dTexp
    dPexp = row.dPexp
    dX1exp = row.dX1exp
    dY1exp = row.dY1exp

    if idExpSet in dict_idExpSet:
        idCmp1 = dict_idExpSet[idExpSet]["idCmp1"]
        idCmp2 = dict_idExpSet[idExpSet]["idCmp2"]
        name1 = dict_idCmp[idCmp1]["name"]
        name2 = dict_idCmp[idCmp2]["name"]
        SMILES1 = dict_idCmp[idCmp1]["SMILES"]
        SMILES2 = dict_idCmp[idCmp2]["SMILES"]

        df.loc[len(df)] = [name1, name2, SMILES1, SMILES2, dTexp, dPexp, dX1exp, dY1exp]

df.to_excel("DDB_Binary_VLE.xlsx")
