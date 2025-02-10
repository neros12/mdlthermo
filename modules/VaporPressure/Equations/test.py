import json
from os.path import join as opj
from pathlib import Path


from NIST_Vp import wagner25_coef

FILE_DIR = Path(__file__).parent

with open(opj(FILE_DIR, "InChIKey_to_CASRN.json"), "r") as f:
    InChIKey_to_CASRN = json.load(f)


CASRN_to_InChIKey = {value: key for key, value in InChIKey_to_CASRN.items()}

parameters = {}
for CASRN in wagner25_coef:
    try:
        InChIKey = CASRN_to_InChIKey[int(CASRN)]
        parameters[InChIKey] = wagner25_coef[CASRN]
    except:
        pass


with open(
    opj(FILE_DIR, "wagner25_parameters.json"), "w", encoding="utf-8"
) as json_file:
    json.dump(parameters, json_file, indent=4, ensure_ascii=False)
