from .calculation_module import calculate_gamma
from .retrieve_comso import retrieve_chemical_profiles


def calculate_binary_gamma(SMILES1, SMILES2, x1, x2, T):
    chemical_profiles = retrieve_chemical_profiles([SMILES1, SMILES2])
    gamma = calculate_gamma(chemical_profiles, [x1, x2], T)

    return gamma
