from typing import Tuple

from .modules import auto_fragment, UNIFAC


def cal_binary_UNIFAC(
    SMILES1: str,
    SMILES2: str,
    x1: float,
    x2: float,
    T: float,
    check_temperature=False,
) -> Tuple[float, float]:
    component1 = auto_fragment(SMILES1)
    component2 = auto_fragment(SMILES2)

    return UNIFAC(
        [component1, component2],
        [x1, x2],
        T,
        check_temperature=check_temperature,
    )
