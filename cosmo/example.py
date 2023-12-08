from cosmosac import CosmoSac
from cosmosac_ml import CosmoSacGcm, CosmoSacMgcn, CosmoSacGcgcn
from util import get_comp_list, get_file

# =============================================================================
# Calculating activity coefficients
# =============================================================================

# Importing COSMO-SAC models
model = {
    "orig": CosmoSac(2010),
    "gcm": CosmoSacGcm(),
    "mgcn": CosmoSacMgcn(),
    "gcgcn": CosmoSacGcgcn(),
}

# Finding molecules from database
database = get_comp_list()
hexane = (database[database["SMILES"] == "CCCCCC"].index)[0]
ethanol = (database[database["SMILES"] == "CCO"].index)[0]

# Calculating activity coeffieicnts
gamma = dict.fromkeys(model.keys(), [])

for name in model:
    m = model[name]

    if name == "orig":
        m.add_comp(get_file(hexane), name="hexane")
        m.add_comp(get_file(ethanol), name="ethanol")
    else:
        m.add_comp("CCCCCC", name="hexane")
        m.add_comp("CCO", name="ethanol")
    m.x = [0.5, 0.5]
    m.T = 373.15  # K

    gamma[name] = m.gam()

# Example
gcgcn_model = CosmoSacGcgcn()
gcgcn_model.add_comp("CCCCCC")
