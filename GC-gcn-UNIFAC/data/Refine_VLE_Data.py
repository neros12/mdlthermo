import pandas as pd

df = pd.read_excel("DDB2007Original.xlsx", sheet_name="HC_MIXT01_DATA_VLE")
list_id = pd.read_excel("HC_MIXT01_ID_VLE01.xlsx")

list_id_set = list_id["idExpSet"].tolist()
df_filtered = df[df["idExpSet"].isin(list_id_set)]
df_filtered.reset_index(drop=True, inplace=True)

df_filtered.to_excel("HC_MIXT01_DATA_VLE.xlsx", index=False)
