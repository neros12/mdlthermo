from os.path import join as opj
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


ROOT_DIR = Path(__file__).resolve().parents[2]
data_path = opj(ROOT_DIR, "data", "DDB_Binary_VLE06.csv")

df = pd.read_csv(data_path)
list_data_set_ID = list(set(df["data_set_ID"].tolist()))

validated_ID = []
for data_set_ID in list_data_set_ID:
    parsed_df = df[df["data_set_ID"] == 1]

    if len(parsed_df) < 5:
        continue

    dx1 = parsed_df["dX1exp"].tolist()
    dy1 = parsed_df["dY1exp"].tolist()
    dT = parsed_df["dTexp"].tolist()
    dP = parsed_df["dPexp"].tolist()

    pass
