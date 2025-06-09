from typing import List

from cheker import have_COSMO, is_only_consist_with
from COSMO_calculation import (
    get_file_dir_from_SMILES,
    retrieve_sigma_profile,
    calculate_sigma_profile,
)


def retrieve_chemical_profiles(molecules: List[str]):
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
