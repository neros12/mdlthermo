import re
from typing import List

import pandas as pd

df = pd.read_excel("HC_MIXT01_ID_VLE.xlsx")
list_id_cmp = pd.read_excel("HC_PURE_PROP01.xlsx")["idCmp"].to_list()

df = df[df["idCmp1"].isin(list_id_cmp)]
df.reset_index(drop=True, inplace=True)
df = df[df["idCmp2"].isin(list_id_cmp)]
df.reset_index(drop=True, inplace=True)

df.to_excel("HC_MIXT01_ID_VLE01.xlsx", index=False)
