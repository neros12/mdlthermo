import pandas as pd

compounds = pd.read_excel("KUSRD_data/KUSRD01_compounds.xlsx")
boiling_points = pd.read_excel("KUSRD_data/KUSRD06_critical_pressures.xlsx")
boiling_points = boiling_points[boiling_points["grade"] > 2]

id_to_name = dict(
    zip(
        compounds["compound_id"].to_list(),
        compounds["name"].to_list(),
    )
)
id_to_SMILES = dict(
    zip(
        compounds["compound_id"].to_list(),
        compounds["smiles"].to_list(),
    )
)

from modules.GCGCN import modules as GCGCN_modules

df = pd.DataFrame(
    columns=[
        "compound_id",
        "name",
        "SMILES",
        "exp_value",
        "JR",
        "MG",
        "GC-GCN",
    ]
)

for row in boiling_points.itertuples():
    compound_id = row.compound_id
    name = id_to_name[compound_id]
    SMILES = id_to_SMILES[compound_id]
    try:
        pred_val, pred_unc = GCGCN_modules.predict_PC(SMILES)
        df.loc[len(df)] = [compound_id, name, SMILES, row.value, "", "", pred_val]
    except:
        pass

df.to_excel("GCGCN_PC_result.xlsx")
