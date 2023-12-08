import os, sys

sys.path.append((os.path.dirname(__file__)))
import DIPPR_Vp
import NIST_Vp

#####################################
# TOTAL NUMBER OF COMPONENTS: 10289 #
#####################################


def CalWagner25(CASRN: str, T: float) -> float:
    parameters = NIST_Vp.wagner25_coef[CASRN]
    Tc = parameters["Tc"]
    lnPr = parameters["lnPr"]
    A = parameters["A"]
    B = parameters["B"]
    C = parameters["C"]
    D = parameters["D"]
    Psat = NIST_Vp.Wagner25(T, Tc, lnPr, A, B, C, D)

    return Psat


def CalEquation101(CASRN: str, T: float) -> float:
    parameters = DIPPR_Vp.equation101_coef[CASRN]
    A = parameters["A"]
    B = parameters["B"]
    C = parameters["C"]
    D = parameters["D"]
    E = parameters["E"]
    Psat = DIPPR_Vp.Equation101(T, A, B, C, D, E)

    return Psat


def CalVaporPressure(CASRN: str, T: float) -> float:
    """
    Input
    ------
      CAS Registry Number\n
      Temperature (K)
    Return
    -------
      Psat (kPa)
    """
    if (
        CASRN in DIPPR_Vp.equation101_coef
        and T > DIPPR_Vp.equation101_coef[CASRN]["Tmin"]
        and T < DIPPR_Vp.equation101_coef[CASRN]["Tmax"]
    ):
        parameters = DIPPR_Vp.equation101_coef[CASRN]
        A = parameters["A"]
        B = parameters["B"]
        C = parameters["C"]
        D = parameters["D"]
        E = parameters["E"]
        Psat = DIPPR_Vp.Equation101(T, A, B, C, D, E)
    elif (
        CASRN in NIST_Vp.wagner25_coef
        and T > NIST_Vp.wagner25_coef[CASRN]["Tmin"]
        and T < NIST_Vp.wagner25_coef[CASRN]["Tmax"]
    ):
        parameters = NIST_Vp.wagner25_coef[CASRN]
        Tc = parameters["Tc"]
        lnPr = parameters["lnPr"]
        A = parameters["A"]
        B = parameters["B"]
        C = parameters["C"]
        D = parameters["D"]
        Psat = NIST_Vp.Wagner25(T, Tc, lnPr, A, B, C, D)
    else:
        Psat = 0.0

    return Psat


# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     import numpy as np

#     CASRN = "64175"
#     Ts = np.arange(100, 310, 10)
#     Y1s = np.array([])
#     Y2s = np.array([])
#     for T in Ts:
#         Y1 = CalWagner25(CASRN, T)
#         Y2 = CalEquation101(CASRN, T)
#         Y1s = np.append(Y1s, Y1)
#         Y2s = np.append(Y2s, Y2)

#     plt.plot(Ts, Y1s, Ts, Y2s)
#     plt.show()
