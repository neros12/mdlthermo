import pandas as pd

df = pd.read_excel("DDB2007Original.xlsx", sheet_name="HC_PURE00_PROP")
list_id = pd.read_excel("HC_MIXT01_ID_VLE.xlsx")

list_id_set = list(set([*list_id["idCmp1"].tolist(), *list_id["idCmp2"].tolist()]))
df_filtered = df[df["idCmp"].isin(list_id_set)]
df_filtered.reset_index(drop=True, inplace=True)

df_filtered.to_excel("HC_PURE_PROP.xlsx", index=False)
