import json
from os.path import join as opj
from pathlib import Path

from rdkit import Chem

from . import wagner_equation


FILE_DIR = Path(__file__).parent
with open(opj(FILE_DIR, "InChIKey_to_CASRN.json"), "r") as f:
    InChIKey_to_CASRN = json.load(f)


def SMILES_to_CASRN(SMILES: str) -> str:
    """
    Convert a SMILES string to a CAS Registry Number (CASRN) using InChIKey as an intermediate.

    Args:
        SMILES (str): The SMILES representation of the chemical compound.

    Returns:
        str: The corresponding CASRN.

    Raises:
        ValueError: If the molecule's InChIKey is not found in the mapping.
    """
    mol = Chem.MolFromSmiles(SMILES)
    InChIKey = Chem.inchi.MolToInchiKey(mol)

    if InChIKey in InChIKey_to_CASRN:
        return str(InChIKey_to_CASRN[InChIKey])
    else:
        raise ValueError("Unsupported Molecule!")


def cal_Wagner25(CASRN: str, T: float, check_temperature=True) -> float:
    """
    Calculate the reduced vapor pressure using the Wagner equation for a given CASRN and temperature.

    Args:
        CASRN (str): CAS Registry Number of the compound.
        T (float): Temperature in Kelvin.
        check_T (bool): Whether to check if the temperature is within the supported range.

    Returns:
        float: The reduced vapor pressure.

    Raises:
        Exception: If the temperature is outside the supported range (when check_T is True).
    """
    parameters = wagner_equation.wagner25_coef[CASRN]
    Tc = parameters["Tc"]
    lnPr = parameters["lnPr"]
    A = parameters["A"]
    B = parameters["B"]
    C = parameters["C"]
    D = parameters["D"]
    Tmin = parameters["Tmin"]
    Tmax = parameters["Tmax"]

    if check_temperature:
        if T < float(Tmin) or T > float(Tmax):
            raise Exception(
                f"""
                The supported temperature range has been exceeded!
                The supporting temperature range is
                {Tmin} K  ~   {Tmax} K            
                """
            )

    return wagner_equation.Wagner25(T, Tc, lnPr, A, B, C, D)


def get_temperature_range(CASRN1: str, CASRN2: str) -> tuple[float, float]:
    """
    Get the overlapping supported temperature range between two compounds.

    Args:
        CASRN1 (str): CASRN of the first compound.
        CASRN2 (str): CASRN of the second compound.

    Returns:
        tuple[float, float]: The (Tmin, Tmax) tuple representing the shared valid temperature range in Kelvin.
    """
    CASRN1_NIST = wagner_equation.wagner25_coef[CASRN1]
    CASRN2_NIST = wagner_equation.wagner25_coef[CASRN2]
    Tmin = max(CASRN1_NIST["Tmin"], CASRN2_NIST["Tmin"])
    Tmax = min(CASRN1_NIST["Tmax"], CASRN2_NIST["Tmax"])

    return Tmin, Tmax
