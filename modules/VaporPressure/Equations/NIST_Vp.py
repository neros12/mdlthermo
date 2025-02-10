import numpy as np


def Wagner25(
    T: float,
    Tc: float,
    lnPr: float,
    A: float,
    B: float = 0.0,
    C: float = 0.0,
    D: float = 0.0,
) -> float:
    """
    Input
    ------
      Temperature (K)\n
      Critical Temperature (K)\n
      lnPr (kPa)
      A\n
      B\n
      C\n
      D\n

    Return
    -------
      Psat (kPa)
    """

    if B == "":
        B = 0.0
    if C == "":
        C = 0.0
    if D == "":
        D = 0.0

    tow = 1 - T / Tc
    lnP = Tc / T * (A * tow + B * tow**1.5 + C * tow**2.5 + D * tow**5) + lnPr

    return np.exp(lnP)
