import os, sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
import vapor_pressure.equations as vp
import NIST_UNIFAC.modules.NIST_UNIFAC as unifac


def Modifed_Raults_Law(T, others):
    c1 = others["c1"]
    c2 = others["c2"]
    x1 = others["x1"]
    x2 = others["x2"]
    activity_model = others["activity_model"]
    Psat = others["vapor_pressure_model"]

    r12, r21 = activity_model(c1["SMILES"], c2["SMILES"], x1, x2, T)
    _left = x1 * r12 * Psat(c1["CASRN"], T)
    _right = x2 * r21 * Psat(c2["CASRN"], T)

    return _left + _right - 101.325


def Bisection(
    xl: float,
    xu: float,
    func,
    max_iter: int = 1000,
    tol: float = 0.0001,
    **kwargs,
) -> float:
    for _iter in range(max_iter):
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


component_1 = {"CASRN": "7732185", "SMILES": "O"}
component_2 = {"CASRN": "71363", "SMILES": "CCCCO"}
activity_model = unifac.cal_activity_coefficient
vapor_pressure_model = vp.Cal_Vapor_Pressure
Patm = 101.325
X1 = np.arange(0, 1.01, 0.01)
Xmin, Xmax = vp.Get_Temp_Range("7732185", "71363")
Xmin = round(Xmin + 5)
Xmax = round(Xmax - 5)

Y1 = []
Ts = []
for x1 in X1:
    x2 = 1 - x1

    Tsol = Bisection(
        Xmin,
        Xmax,
        Modifed_Raults_Law,
        c1=component_1,
        c2=component_2,
        x1=x1,
        x2=x2,
        activity_model=activity_model,
        vapor_pressure_model=vapor_pressure_model,
    )

    Ts.append(Tsol)
    Y1.append(
        x1
        * activity_model(component_1["SMILES"], component_2["SMILES"], x1, x2, Tsol)[0]
        * vapor_pressure_model(component_1["CASRN"], Tsol)
        / Patm
    )

plt.plot(X1, Ts, Y1, Ts)
plt.show()
