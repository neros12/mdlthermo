import re
from typing import List

import pandas as pd

df = pd.read_excel("DDB2007Original.xlsx", sheet_name="HC_MIXT01_ID_VLE")

expNames: List[str] = df["strExpName"].tolist()
list_cmp1 = []
list_cmp2 = []
list_nData = []
for expName in expNames:
    match = re.search(r"\d+ pts", expName)

    nData = match.group().split("pts")[0].strip()

    temp2 = re.sub(r"\d+ pts", "", expName).strip()
    temp3 = temp2.split("+")
    cmp1 = temp3[0].split("of")[1].strip()
    cmp2 = temp3[1].split("at")[0].strip()

    print("compound1:", cmp1, "compound2:", cmp2, "nData:", nData)

    list_cmp1.append(cmp1)
    list_cmp2.append(cmp2)
    list_nData.append(nData)

df["compound1"] = list_cmp1
df["compound2"] = list_cmp2

df.to_excel("HC_MIXT01_ID_VLE.xlsx", index=False)
