from numpy import ndarray
from typing import List, TypedDict

from cheker import have_COSMO, is_only_consist_with
from COSMO_calculation import (
    get_file_dir_from_SMILES,
    retrieve_sigma_profile,
    calculate_sigma_profile,
)


class ChemicalProfile(TypedDict):
    area: float
    volume: float
    sigma_profiles: ndarray


class ChemicalProfiles(TypedDict):
    version: int
    data: List[ChemicalProfile]


def retrieve_chemical_profiles(molecules: List[str]) -> ChemicalProfiles:
    """
    Retrieves or calculates the COSMO-based chemical profiles for a list of molecules.

    For each molecule (represented by its SMILES string), the function attempts to check
    whether precomputed COSMO data is available:

    - If COSMO data is available, it retrieves the sigma profile from a file.
    - If COSMO data is not available, it calculates the sigma profile from scratch.
    - If any molecule does not meet the requirements (e.g., contains unsupported elements),
      the function raises a ValueError.

    If even one molecule lacks COSMO data, the entire profile version is downgraded to 2014;
    otherwise, the version remains 2020.

    Args:
        molecules (List[str]): A list of SMILES strings representing the molecules.

    Returns:
        ChemicalProfiles: A dictionary containing the COSMO version used and
        a list of chemical profile data, each including area, volume, and sigma profile.

    Raises:
        ValueError: If any molecule does not meet the basic compositional requirements.
    """

    chemical_profiles = {"version": 2020, "data": []}

    for SMILES in molecules:
        check_COSMO = have_COSMO(SMILES)
        if not check_COSMO:
            chemical_profiles["version"] = 2014
            if not is_only_consist_with(SMILES):
                raise ValueError("Unsupported Molecule!")

    for SMILES in molecules:
        check_COSMO = have_COSMO(SMILES)
        if check_COSMO:
            COSMO_dir = get_file_dir_from_SMILES(SMILES)
            chemical_profiles["data"].append(
                retrieve_sigma_profile(
                    COSMO_dir,
                    chemical_profiles["version"],
                )
            )
        else:
            chemical_profiles["data"].append(calculate_sigma_profile(SMILES))

    return chemical_profiles
