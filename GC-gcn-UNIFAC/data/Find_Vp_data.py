import os
import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))
from vapor_pressure.equations.NIST_Vp import wagner25_coef

available_casrn = list(map(int, wagner25_coef.keys()))
df = pd.read_excel(os.path.join("data", "HC_PURE_PROP.xlsx"))
df_filtered = df[df["CASRN"].isin(available_casrn)]
df_filtered.reset_index(drop=True, inplace=True)

df_filtered.to_excel("HC_PURE_PROP01.xlsx", index=False)
