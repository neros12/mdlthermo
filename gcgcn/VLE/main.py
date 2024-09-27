import sys
from os.path import join as opj

import numpy as np
from pathlib import Path

ROOT_FOLDER = Path(__file__).parent.parent
if ROOT_FOLDER not in sys.path:
    sys.path.append(str(ROOT_FOLDER))

import cosmosac
import vapor_pressure


def Modifed_Raults_Law(T, others):
    SMILES1 = others["SMILES1"]
    SMILES2 = others["SMILES2"]
    x1 = others["x1"]
    x2 = others["x2"]
    P = others["P"]
    COMSOSAC_models = others["COMSOSAC_models"]
    Psat_Models = others["Psat_Models"]

    r12, r21 = cosmosac.cal_gamma(SMILES1, SMILES2, x1, x2, T, COMSOSAC_models)
    _left = x1 * r12 * vapor_pressure.cal_vapor_pressure(SMILES1, T, Psat_Models)
    _right = x2 * r21 * vapor_pressure.cal_vapor_pressure(SMILES2, T, Psat_Models)

    return _left + _right - P


def Bisection(
    xl: float,
    xu: float,
    func,
    max_iter: int = 100,
    tol: float = 0.0001,
    **kwargs,
) -> float:
    for _iter in range(max_iter):
        print(_iter)
        xr = (xl + xu) / 2

        if len(kwargs) == 0:
            ff = func(xl) * func(xr)
        else:
            ff = func(xl, kwargs) * func(xr, kwargs)

        if abs(ff) < tol:
            break
        elif _iter == max_iter - 1:
            xr = 0.0
        elif ff < 0:
            xu = xr
        elif ff > 0:
            xl = xr

    return xr


def main(SMILES1, SMILES2, T):
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

    return X1, Y1, Ps


if __name__ == "__main__":

    pass


# Xmin = 100.0
# Xmax = 1000.0
# Patm = 101.325
# Ts = []
# for x1 in X1:
#     x2 = 1 - x1

#     Tsol = Bisection(
#         Xmin,
#         Xmax,
#         Modifed_Raults_Law,
#         SMILES1=SMILES1,
#         SMILES2=SMILES2,
#         x1=x1,
#         x2=x2,
#         P=Patm,
#         COMSOSAC_models=COMSOSAC_models,
#         Psat_Models=Psat_Models,
#     )

#     Ts.append(Tsol)
#     Y1.append(
#         x1
#         * cosmosac.cal_gamma(SMILES1, SMILES2, x1, x2, Tsol, COMSOSAC_models)[0]
#         * vapor_pressure.cal_vapor_pressure(SMILES1, Tsol, Psat_Models)
#         / Patm
#     )
# Ts = np.array(Ts) - 273.15
# plt.plot(X1, Ts, Y1, Ts)
