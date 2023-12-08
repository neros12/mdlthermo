import os, sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
import vapor_pressure_equations as vp
import NIST_UNIFAC.modules.NIST_UNIFAC as unifac


def modifed_raults_law(T, others):
    c1 = others["c1"]
    c2 = others["c2"]
    x1 = others["x1"]
    x2 = others["x2"]
    activity_model = others["activity_model"]
    Psat = others["Psat"]

    r12, r21 = activity_model(c1["SMILES"], c2["SMILES"], x1, x2, T)
    _left = x1 * r12 * Psat(c1["CASRN"], T)
    _right = x2 * r21 * Psat(c2["CASRN"], T)

    print(_left, _right)

    return _left + _right - 101.325


def bisection(
    xl: float,
    xu: float,
    func,
    max_iter: int = 1000,
    tol: float = 0.001,
    **kwargs,
) -> float:
    for _iter in range(max_iter):
        xr = (xl + xu) / 2

        if len(kwargs) == 0:
            ff = func(xl) * func(xr)
        else:
            ff = func(xl, kwargs) * func(xr, kwargs)

        print(xl, xu, xr, func(xl, kwargs), func(xr, kwargs))

        if abs(ff) < tol:
            break
        elif _iter == max_iter - 1:
            xr = 0.0
        elif ff < 0:
            xu = xr
        elif ff > 0:
            xl = xr

    return xr


Water = {"CASRN": "7732185", "SMILES": "O"}
Toluene = {"CASRN": "108883", "SMILES": "CC1=CC=CC=C1"}

# a = vp.CalVaporPressure(Toluene["CASRN"], 373.15)
# b = unifac.cal_activity_coefficient(
#     Water["SMILES"], Toluene["SMILES"], 0.1, 0.2, 293.15
# )

P = 101.325
X1 = np.arange(0, 1.05, 0.05)
X2 = 1 - X1
Xmin = 273
Xmax = 350

a = bisection(
    Xmin,
    Xmax,
    modifed_raults_law,
    c1=Water,
    c2=Toluene,
    x1=0.1,
    x2=0.2,
    activity_model=unifac.cal_activity_coefficient,
    Psat=vp.CalVaporPressure,
)

print(a)
