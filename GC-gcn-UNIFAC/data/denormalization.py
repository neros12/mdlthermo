import os
import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))
from vapor_pressure.equations.NIST_Vp import Wagner25

df = pd.read_excel("HC_MIXT01_DATA_VLE.xlsx")
idCmp_df = pd.read_excel("HC_PURE_PROP01.xlsx")
idExpSet = pd.read_excel("HC_MIXT01_ID_VLE01.xlsx")

idCmpToCompound = {}
idCmpToCASRN = {}
for row in idCmp_df.itertuples():
    idCmpToCompound[row.idCmp] = row.strCmpName1
    idCmpToCASRN[row.idCmp] = str(row.CASRN)

idExpSetToCompound = {}
for row in idExpSet.itertuples():
    idExpSetToCompound[row.idExpSet] = (row.idCmp1, row.idCmp2)

for row in df.itertuples():
    idExpSet = row.idExpSet
    idCmp1, idCmp2 = idExpSetToCompound[idExpSetToCompound]
    compound1 = idCmpToCompound[idCmp1]
    compound2 = idCmpToCompound[idCmp2]
    casrn1 = idCmpToCompound[idCmp1]
    casrn2 = idCmpToCompound[idCmp2]
    Psat1 = Wagner25()
