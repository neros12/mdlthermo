import DIPPR_Vp
import NIST_Vp


#####################################
# TOTAL NUMBER OF COMPONENTS: 10289 #
#####################################
def cal_Wagner25(CASRN: str, T: float) -> float:
    parameters = NIST_Vp.wagner25_coef[CASRN]
    Tc = parameters["Tc"]
    lnPr = parameters["lnPr"]
    A = parameters["A"]
    B = parameters["B"]
    C = parameters["C"]
    D = parameters["D"]
    Psat = NIST_Vp.Wagner25(T, Tc, lnPr, A, B, C, D)

    return Psat


def cal_Equation101(CASRN: str, T: float) -> float:
    parameters = DIPPR_Vp.equation101_coef[CASRN]
    A = parameters["A"]
    B = parameters["B"]
    C = parameters["C"]
    D = parameters["D"]
    E = parameters["E"]
    Psat = DIPPR_Vp.Equation101(T, A, B, C, D, E)

    return Psat


def get_temp_range(CASRN1: str, CASRN2: str) -> tuple[float, float]:
    CASRN1_DIPPR = DIPPR_Vp.equation101_coef[CASRN1]
    CASRN2_DIPPR = DIPPR_Vp.equation101_coef[CASRN2]
    CASRN1_NIST = NIST_Vp.wagner25_coef[CASRN1]
    CASRN2_NIST = NIST_Vp.wagner25_coef[CASRN2]

    Tmin = max(
        CASRN1_DIPPR["Tmin"],
        CASRN2_DIPPR["Tmin"],
        CASRN1_NIST["Tmin"],
        CASRN2_NIST["Tmin"],
    )
    Tmax = min(
        CASRN1_DIPPR["Tmax"],
        CASRN2_DIPPR["Tmax"],
        CASRN1_NIST["Tmax"],
        CASRN2_NIST["Tmax"],
    )

    return Tmin, Tmax


def cal_vapor_pressure(CASRN: str, T: float) -> float:
    """
    Input
    ------
      CASRN: CAS Registry Number\n
      T: Temperature (K)
    Return
    -------
      Psat(T) (kPa)
    """

    if CASRN not in DIPPR_Vp.equation101_coef and CASRN not in NIST_Vp.wagner25_coef:
        raise Exception(f"There is no coeffcient for CAS Registery Number: {CASRN}")

    if (
        CASRN
        in DIPPR_Vp.equation101_coef
        # and T > DIPPR_Vp.equation101_coef[CASRN]["Tmin"]
        # and T < DIPPR_Vp.equation101_coef[CASRN]["Tmax"]
    ):
        parameters = DIPPR_Vp.equation101_coef[CASRN]
        A = parameters["A"]
        B = parameters["B"]
        C = parameters["C"]
        D = parameters["D"]
        E = parameters["E"]
        Psat = DIPPR_Vp.Equation101(T, A, B, C, D, E)
    elif (
        CASRN
        in NIST_Vp.wagner25_coef
        # and T > NIST_Vp.wagner25_coef[CASRN]["Tmin"]
        # and T < NIST_Vp.wagner25_coef[CASRN]["Tmax"]
    ):
        parameters = NIST_Vp.wagner25_coef[CASRN]
        Tc = parameters["Tc"]
        lnPr = parameters["lnPr"]
        A = parameters["A"]
        B = parameters["B"]
        C = parameters["C"]
        D = parameters["D"]
        Psat = NIST_Vp.Wagner25(T, Tc, lnPr, A, B, C, D)

    # if (
    #     CASRN in DIPPR_Vp.equation101_coef
    #     and T > DIPPR_Vp.equation101_coef[CASRN]["Tmin"]
    #     and T < DIPPR_Vp.equation101_coef[CASRN]["Tmax"]
    # ):
    #     parameters = DIPPR_Vp.equation101_coef[CASRN]
    #     A = parameters["A"]
    #     B = parameters["B"]
    #     C = parameters["C"]
    #     D = parameters["D"]
    #     E = parameters["E"]
    #     Psat = DIPPR_Vp.Equation101(T, A, B, C, D, E)
    # elif (
    #     CASRN in NIST_Vp.wagner25_coef
    #     and T > NIST_Vp.wagner25_coef[CASRN]["Tmin"]
    #     and T < NIST_Vp.wagner25_coef[CASRN]["Tmax"]
    # ):
    #     parameters = NIST_Vp.wagner25_coef[CASRN]
    #     Tc = parameters["Tc"]
    #     lnPr = parameters["lnPr"]
    #     A = parameters["A"]
    #     B = parameters["B"]
    #     C = parameters["C"]
    #     D = parameters["D"]
    #     Psat = NIST_Vp.Wagner25(T, Tc, lnPr, A, B, C, D)
    # else:
    #     raise Exception(
    #         f"""
    #         The given temperature ({T} K) is outside the valid range of the equations for {CASRN}.
    #         DIPPR: {DIPPR_Vp.equation101_coef[CASRN]['Tmin']} to {DIPPR_Vp.equation101_coef[CASRN]['Tmax']}
    #         NIST: {NIST_Vp.wagner25_coef[CASRN]['Tmin']} to {NIST_Vp.wagner25_coef[CASRN]['Tmax']}
    #         """
    #     )

    return Psat
