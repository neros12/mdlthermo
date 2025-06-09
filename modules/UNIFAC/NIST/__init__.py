from typing import List

from .modules import auto_fragmentation, calculate_gamma


def calculate_binary_gamma(
    SMILES1: str,
    SMILES2: str,
    x1: float,
    x2: float,
    T: float,
) -> List[float]:
    """
    Calculate activity coefficients for a binary mixture using the NIST-UNIFAC model.

    This function retrieves UNIFAC groups for two components
    defined by their SMILES strings, then computes their activity coefficients (γ)
    at a given temperature and composition.

    Parameters
    ----------
    SMILES1 : str
        SMILES string of the first component.
    SMILES2 : str
        SMILES string of the second component.
    x1 : float
        Mole fraction of the first component.
    x2 : float
        Mole fraction of the second component.
    T : float
        Temperature in Kelvin.

    Returns
    -------
    gamma : List[float]
        A list containing activity coefficients [γ1, γ2] for the two components.
    """
    cp1 = auto_fragmentation(SMILES1)
    cp2 = auto_fragmentation(SMILES2)
    gamma = calculate_gamma([cp1, cp2], [x1, x2], T)

    return gamma
