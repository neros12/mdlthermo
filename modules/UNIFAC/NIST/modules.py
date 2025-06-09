from typing import List

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

from . import parameters


def auto_fragmentation(SMILES: str) -> List[int]:
    """
    Fragment a molecule into UNIFAC subgroups based on SMARTS matching.

    This function takes a SMILES string, converts it into an RDKit molecule,
    and attempts to match predefined SMARTS patterns to identify functional
    subgroups according to the UNIFAC model.

    Parameters
    ----------
    SMILES : str
        The SMILES representation of a molecule.

    Returns
    -------
    List[int]
        A list of subgroup IDs matched from the molecule.

    Raises
    ------
    ValueError
        If the SMILES string is invalid or the molecule contains atoms not covered
        by the current UNIFAC subgroup definitions.
    """
    try:
        mol = Chem.MolFromSmiles(SMILES)
    except:
        raise ValueError("Invalid SMILES format")

    if mol == None:
        raise ValueError("Invalid SMILES format")

    num_atoms = Descriptors.HeavyAtomCount(mol)
    list_atoms = np.arange(num_atoms)
    subgroups = []

    for item in parameters.subgroup_SMARTS:
        subgroup_ID = item[0]
        SMARTS = item[1]

        # These UNIFAC groups are not supported in this version.
        if SMARTS == "EX":
            continue

        if mol.HasSubstructMatch(Chem.MolFromSmarts(SMARTS)) == True:
            matched_fragments = mol.GetSubstructMatches(Chem.MolFromSmarts(SMARTS))
            for matched_atoms in matched_fragments:
                if set(matched_atoms) & set(list_atoms) == set(matched_atoms):
                    subgroups.append(subgroup_ID)
                    list_atoms = np.setdiff1d(list_atoms, matched_atoms)

        if len(list_atoms) == 0:
            break

    if len(list_atoms) != 0:
        raise ValueError("Molecule is not only consist of UNIFAC groups")

    return subgroups


def calculate_gamma(
    components: List[List[int]],
    x: List[float],
    T: float,
) -> List[float]:
    """
    Calculate activity coefficients using the UNIFAC group contribution method.

    This function computes activity coefficients (γ) for a mixture of components
    based on their subgroup decomposition. It includes both combinatorial and
    residual contributions and checks the validity of temperature and group interactions.

    Parameters
    ----------
    components : List[List[int]]
        A list of components, where each component is a list of UNIFAC subgroup IDs.
    x : List[float]
        Mole fractions of each component in the mixture.
    T : float
        Temperature in Kelvin.

    Returns
    -------
    List[float]
        The activity coefficients for each component.

    Raises
    ------
    ValueError
        If the temperature is outside the valid range, or if required group interaction
        parameters are missing.
    """
    total_subgroups = []
    for component in components:
        total_subgroups.extend(component)
    total_subgroups = set(total_subgroups)
    total_maingroups = {
        parameters.sub_to_main[subgroup] for subgroup in total_subgroups
    }

    Tmax = []
    Tmin = []
    iter_maingroups = []
    for m in total_maingroups:
        for n in total_maingroups:
            if m == n:
                continue
            else:
                if (m, n) in parameters.maingroup_par.keys():
                    iter_maingroups.append((m, n))
                    Tmax.append(parameters.maingroup_par[(m, n)]["Tmax"])
                    Tmin.append(parameters.maingroup_par[(m, n)]["Tmin"])
                else:
                    raise ValueError(
                        f"There are no interaction parameters between maingroup {m} and maingroup {n}."
                    )

    if len(Tmax) != 0:
        Tmax = max(Tmax)
        Tmin = min(Tmin)
        if T > Tmax or T < Tmin:
            raise ValueError(
                "Current method is unsupported for the given temperature condition."
            )

    num_components = len(x)
    vk = []
    for component in components:
        _vk = {}
        for k in total_subgroups:
            _vk[k] = component.count(k)
        vk.append(_vk)

    # Combinatorial Term
    r = []  # index = component
    q = []  # index = component
    for index in range(num_components):
        SumR, SumQ = 0.0, 0.0
        for k in components[index]:
            SumR += parameters.subgroup_par[k]["Ri"]
            SumQ += parameters.subgroup_par[k]["Qi"]
        r.append(SumR)
        q.append(SumQ)

    J = []  # index = component
    L = []  # index = component
    for i in range(num_components):
        rx = 0.0
        qx = 0.0
        for j in range(num_components):
            rx += r[j] ** (3 / 4) * x[j]
            qx += q[j] * x[j]
        J.append(r[i] ** (3 / 4) / rx)
        L.append(q[i] / qx)

    ln_rc = []  # index = component
    for i in range(num_components):
        ln_rc.append(
            1 - J[i] + np.log(J[i]) - 5 * q[i] * (1 - J[i] / L[i] + np.log(J[i] / L[i]))
        )

    # Residual Term
    e = {}  # key = (subgroup, component)
    for k in total_subgroups:
        for i in range(num_components):
            e[(k, i)] = vk[i][k] * parameters.subgroup_par[k]["Qi"] / q[i]

    theta = {}  # key = subgroup
    for k in total_subgroups:
        xqe = 0
        xq = 0
        for i in range(num_components):
            xqe += x[i] * q[i] * e[(k, i)]
            xq += x[i] * q[i]
        theta[k] = xqe / xq

    u = {}  # key = (subgroup, subgroup)
    for m in total_subgroups:
        for k in total_subgroups:
            i = parameters.sub_to_main[m]
            j = parameters.sub_to_main[k]
            if i == j:
                u[(m, k)] = 0.0
            else:
                u[(m, k)] = (
                    parameters.maingroup_par[(i, j)]["a1"]
                    + parameters.maingroup_par[(i, j)]["a2"] * T
                    + parameters.maingroup_par[(i, j)]["a3"] / 1000 * T * T
                )

    tow = {}  # key = (subgroup, subgroup)
    for m in total_subgroups:
        for k in total_subgroups:
            tow[(m, k)] = np.exp(-u[(m, k)] / T)

    betta = {}  # key = (component, subgroup)
    for i in range(num_components):
        for k in total_subgroups:
            sumbetta = 0
            for m in total_subgroups:
                sumbetta += e[(m, i)] * tow[(m, k)]
            betta[(i, k)] = sumbetta

    s = {}  # key = subgroup
    for k in total_subgroups:
        sums = 0
        for m in total_subgroups:
            sums += theta[m] * tow[(m, k)]
        s[k] = sums

    ln_rr = []  # index = component
    for i in range(num_components):
        sum_residual = 0
        for k in total_subgroups:
            sum_residual += theta[k] * betta[(i, k)] / s[k] - e[(k, i)] * np.log(
                betta[(i, k)] / s[k]
            )
        ln_rr.append(q[i] * (1 - sum_residual))

    ln_r = []  # key = component
    for i in range(num_components):
        ln_r.append(ln_rc[i] + ln_rr[i])

    activity_coefficient = np.exp(ln_r)

    return activity_coefficient.tolist()
