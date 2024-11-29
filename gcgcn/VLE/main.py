import sys
from os.path import join as opj

import numpy as np
from pathlib import Path

ROOT_FOLDER = Path(__file__).parent.parent
if ROOT_FOLDER not in sys.path:
    sys.path.append(str(ROOT_FOLDER))

import cosmosac
import vapor_pressure


def cal_isotherm_VLE(SMILES1, SMILES2, T):
    COMSOSAC_models = cosmosac.load_models()
    Psat_Models = vapor_pressure.load_models()

    X1 = np.arange(0, 1.05, 0.05)
    Y1 = []
    Ps = []
    for x1 in X1:
        x2 = 1 - x1

        r1, r2 = cosmosac.cal_gamma(SMILES1, SMILES2, x1, x2, T, COMSOSAC_models)
        Psat1 = vapor_pressure.cal_vapor_pressure(SMILES1, T, Psat_Models)
        Psat2 = vapor_pressure.cal_vapor_pressure(SMILES2, T, Psat_Models)

        P = x1 * r1 * Psat1 + x2 * r2 * Psat2
        Ps.append(P)
        Y1.append(x1 * r1 * Psat1 / P)

    Ps = (np.array(Ps) / 101.325).tolist()
    X1 = X1.tolist()

    return X1, Y1, Ps


def cal_isobaric_VLE(SMILES1, SMILES2, P=101.325):
    COMSOSAC_models = cosmosac.load_models()
    Psat_Models = vapor_pressure.load_models()

    X1 = np.arange(0, 1.05, 0.05)
    Y1 = []
    Ts = []
    for x1 in X1:
        x2 = 1 - x1
        Tl = 100.0
        Tu = 1000.0
        for i in range(500):
            Tm = (Tl + Tu) / 2
            r1, r2 = cosmosac.cal_gamma(SMILES1, SMILES2, x1, x2, Tm, COMSOSAC_models)
            Psat1 = vapor_pressure.cal_vapor_pressure(SMILES1, Tm, Psat_Models)
            Psat2 = vapor_pressure.cal_vapor_pressure(SMILES2, Tm, Psat_Models)
            check_point = x1 * r1 * Psat1 + x2 * r2 * Psat2 - P

            if abs(check_point) < 0.001:
                break
            elif check_point < 0:
                Tl = Tm
            elif check_point > 0:
                Tu = Tm

        Ts.append(Tm)
        Y1.append(x1 * r1 * Psat1 / P)

    X1 = X1.tolist()

    return X1, Y1, Ts


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    X1, Y1, Ts = cal_isobaric_VLE("c1ccccc1", "CCCCCCCC", P=101)

    print()
    print(X1)
    print()
    print(Y1)
    print()
    print(Ts)

    plt.plot(X1, Ts)
    plt.plot(Y1, Ts)
    plt.show()
