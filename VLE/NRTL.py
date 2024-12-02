import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
import vapor_pressure.equations as vp


A12 = 2.3015
A21 = -0.5363
B12 = -879.7008
B21 = 1412.7316
C = 0.3


def NRTL(x1, x2, T, a12, a21, b12, b21, c):
    tow12 = a12 + b12 / (T + 273.15)
    tow21 = a21 + b21 / (T + 273.15)

    x = [x1, x2]
    tow = np.array([[0, tow12], [tow21, 0]])
    G = np.exp(-c * tow)

    gamma = np.zeros(2)
    for i in range(2):
        a = 0
        for j in range(2):
            a += x[j] * tow[j][i] * G[j][i]
        b = 0
        for k in range(2):
            b += x[j] * G[j][i]
        c = a / b
        d = 0
        for j in range(2):
            e = 0
            for k in range(2):
                e += x[k] * G[k][j]
            f = x[j] * G[i][j] / e

            g = 0
            for m in range(2):
                g += x[m] * tow[m][j] * G[m][j]
            h = 0
            for k in range(2):
                h += x[k] * G[k][j]

            d += f * (tow[i][j] - g / h)

        gamma[i] = np.exp(d + c)

    return gamma


def Modifed_Raults_Law(T, others):
    c1 = others["c1"]
    c2 = others["c2"]
    x1 = others["x1"]
    x2 = others["x2"]
    activity_model = others["activity_model"]
    Psat = others["vapor_pressure_model"]

    r12, r21 = activity_model(
        x1, x2, T, 7.79716106, -0.107372719, -3131.20054, 1540.90442, 0.3
    )
    _left = x1 * r12 * Psat(c1["CASRN"], T)
    _right = x2 * r21 * Psat(c2["CASRN"], T)

    return _left + _right - 101.325


def Bisection(
    xl: float,
    xu: float,
    func,
    max_iter: int = 10000,
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


# component_2 = {"CASRN": "7732185"}
# component_1 = {"CASRN": "108952"}
# vapor_pressure_model = vp.Cal_Vapor_Pressure
# Patm = 101.325
# X1 = np.arange(0.01, 1, 0.01)
# Xmin = 315
# Xmax = 640

# Y1 = []
# Ts = []
# for x1 in [0.88]:
#     x2 = 1 - x1

#     Tsol = Bisection(
#         Xmin,
#         Xmax,
#         Modifed_Raults_Law,
#         c1=component_1,
#         c2=component_2,
#         x1=x1,
#         x2=x2,
#         activity_model=NRTL,
#         vapor_pressure_model=vapor_pressure_model,
#     )

#     if Tsol == 0:
#         pass

#     Ts.append(Tsol)
#     Y1.append(
#         x1
#         * NRTL(x1, x2, Tsol, 7.79716106, -0.107372719, -3131.20054, 1540.90442, 0.3)
#         * vapor_pressure_model(component_1["CASRN"], Tsol)
#         / Patm
#     )

# plt.plot(X1, Ts, Y1, Ts)
# plt.show()

a, b = NRTL(0.5, 0.5, 300, 0, 0, 0, 0, 0)

print(a)
