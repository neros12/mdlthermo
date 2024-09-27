import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
import os, sys

RDLogger.DisableLog("rdApp.*")
sys.path.append(os.path.dirname(__file__))
import parameters


def auto_fragmentation(SMILES: str) -> list[int]:
    mol = Chem.MolFromSmiles(SMILES)

    if mol == None:
        raise ValueError("Unsupported or invalid SMILES format.")

    num_atoms = Descriptors.HeavyAtomCount(mol)
    list_atoms = np.arange(num_atoms)
    subgroups = []

    for item in parameters.subgroup_SMARTS:
        subgroup_ID = item[0]
        SMARTS = item[1]

        if SMARTS == "EX":
            # These UNIFAC groups are not supported in this version.
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
        raise ValueError("Molecule is not only consist of UNIFAC subgroups.")

    return subgroups


def cal_activity_coefficient(
    SMILES1: str, SMILES2: str, x1: float, x2: float, T: float
) -> list[float, float]:
    ####################
    # Check Conditions #
    ####################
    c1 = []
    try:
        c1 = auto_fragmentation(SMILES1)
    except Exception as ex:
        raise ValueError(f"Fragmentation of component 1 has been failed: {ex}")
    c2 = []
    try:
        c2 = auto_fragmentation(SMILES2)
    except Exception as ex:
        raise ValueError(f"Fragmentation of component 2 has been failed: {ex}")

    total_subgroups = {*c1, *c2}
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

    ####################
    # Calcuation Start #
    ####################
    components = [c1, c2]
    x = [x1, x2]
    vk1 = {}
    vk2 = {}
    for k in total_subgroups:
        vk1[k] = c1.count(k)
        vk2[k] = c2.count(k)
    vk = [vk1, vk2]

    # Combinatorial Term
    r = []  # index = component
    q = []  # index = component
    for index in [0, 1]:
        SumR, SumQ = 0.0, 0.0
        for k in components[index]:
            SumR += parameters.subgroup_par[k]["Ri"]
            SumQ += parameters.subgroup_par[k]["Qi"]
        r.append(SumR)
        q.append(SumQ)

    J = []  # index = component
    L = []  # index = component
    for i in [0, 1]:
        rx = 0.0
        qx = 0.0
        for j in [0, 1]:
            rx += r[j] ** (3 / 4) * x[j]
            qx += q[j] * x[j]
        J.append(r[i] ** (3 / 4) / rx)
        L.append(q[i] / qx)

    ln_rc = []  # index = component
    for i in [0, 1]:
        ln_rc.append(
            1 - J[i] + np.log(J[i]) - 5 * q[i] * (1 - J[i] / L[i] + np.log(J[i] / L[i]))
        )

    # Residual Term
    e = {}  # key = (subgroup, component)
    for k in total_subgroups:
        for i in [0, 1]:
            e[(k, i)] = vk[i][k] * parameters.subgroup_par[k]["Qi"] / q[i]

    theta = {}  # key = subgroup
    for k in total_subgroups:
        xqe = 0
        xq = 0
        for i in [0, 1]:
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
    for i in [0, 1]:
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
    for i in [0, 1]:
        sum_residual = 0
        for k in total_subgroups:
            sum_residual += theta[k] * betta[(i, k)] / s[k] - e[(k, i)] * np.log(
                betta[(i, k)] / s[k]
            )
        ln_rr.append(q[i] * (1 - sum_residual))

    ln_r = []  # key = component
    for i in [0, 1]:
        ln_r.append(ln_rc[i] + ln_rr[i])

    activity_coefficient = np.exp(ln_r)

    return activity_coefficient
