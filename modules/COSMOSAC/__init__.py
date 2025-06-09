# from .modules import COSMOSAC
from modules import COSMOSAC, get_COSMO_file_dir


def cal_COSMOSAC(SMILES1, SMILES2, x1, x2, T, predict=False):
    if predict:
        COSMOSAC_module = COSMOSAC(version=2010, predict=True)
    else:
        COSMOSAC_module = COSMOSAC(version=2019, predict=False)

    COSMOSAC_module.add_comp(SMILES=SMILES1)
    COSMOSAC_module.add_comp(SMILES=SMILES2)
    COSMOSAC_module.x = [x1, x2]
    COSMOSAC_module.T = T

    return COSMOSAC_module.gam()


COSMOSAC_module = COSMOSAC(version=2010, predict=False)

dic = {}

COSMOSAC_module.add_comp(SMILES="O")
dic["O"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="C")
dic["C"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="CC")
dic["CC"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="C=C")
dic["C=C"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="C#C")
dic["C#C"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="CO")
dic["CO"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="C=O")
dic["C=O"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="CN")
dic["CN"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="C#N")
dic["C#N"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="FF")
dic["FF"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="BrBr")
dic["BrBr"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="ClCl")
dic["ClCl"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

COSMOSAC_module.add_comp(SMILES="O=C=O")
dic["O=C=O"] = {
    "area": COSMOSAC_module.A[0],
    "volume": COSMOSAC_module.V[0],
    "sigma_profiles": COSMOSAC_module.psigA.tolist(),
}
COSMOSAC_module.del_comp()

import json

with open("data.json", "w", encoding="utf-8") as f:
    json.dump(dic, f, indent=4)
