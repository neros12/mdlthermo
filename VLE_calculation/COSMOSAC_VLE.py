import os
import sys
import numpy as np

import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
from vapor_pressure import equations as vp
from cosmo import cosmosac_ml


def main(
    SMILES_1="CCCCCC",
    SMILES_2="CC(C)O",
    CASRN_1="110543",
    CASRN_2="67630",
    P=101.315,
    max_iter=100,
    tol=1e-5,
):
    x1 = np.arange(0.01, 1, 0.01)
    y1 = np.arange(0.01, 1, 0.01)
    Tx = []
    Ty = []

    T_min, T_max = vp.Get_Temp_Range(CASRN_1, CASRN_2)
    T_med = (T_min + T_max) * 0.5

    cosmosac_class = cosmosac_ml.CosmoSacGcgcn()
    cosmosac_class.add_comp(SMILES_1)
    cosmosac_class.add_comp(SMILES_2)

    # Find molecule_1 boiling point
    Tb_1 = T_med
    for _iter in range(max_iter):
        check_point = P - vp.Cal_Vapor_Pressure(CASRN_1, Tb_1)
        if abs(check_point) < tol:
            break
        else:
            Tb_1 = Tb_1 - tol * check_point / (
                (P - vp.Cal_Vapor_Pressure(CASRN_1, Tb_1 + tol)) - check_point
            )

    # Find molecule_2 boiling point
    Tb_2 = T_med
    for _iter in range(max_iter):
        check_point = P - vp.Cal_Vapor_Pressure(CASRN_2, Tb_2)
        if abs(check_point) < tol:
            break
        else:
            Tb_2 = Tb_2 - tol * check_point / (
                (P - vp.Cal_Vapor_Pressure(CASRN_2, Tb_2 + tol)) - check_point
            )

    gamma_list = []
    for xi in x1:
        # Find Tx with Bisection Method
        Tl = max(50, min(Tb_1 - 100, Tb_2 - 100))
        Tu = min(Tb_1 + 100, Tb_2 + 100)
        cosmosac_class.x = [xi, 1 - xi]
        for _iter in range(max_iter):
            if _iter == max_iter - 1:
                Tx.append(0)
            Tr = (Tl + Tu) * 0.5

            cosmosac_class.T = Tl
            ac12, ac21 = cosmosac_class.gam()
            LT = (
                xi * vp.Cal_Vapor_Pressure(CASRN_1, Tl) * ac12
                + (1 - xi) * ac21 * vp.Cal_Vapor_Pressure(CASRN_2, Tl)
                - P
            )
            cosmosac_class.T = Tr
            ac12, ac21 = cosmosac_class.gam()
            RT = (
                xi * vp.Cal_Vapor_Pressure(CASRN_1, Tr) * ac12
                + (1 - xi) * ac21 * vp.Cal_Vapor_Pressure(CASRN_2, Tr)
                - P
            )

            check_point = LT * RT
            if abs(check_point) < tol:
                Tx.append(Tr)
                gamma_list.append([ac12, ac21])
                break
            elif check_point < 0:
                Tu = Tr
            else:
                Tl = Tr

    for yi in y1:
        # Find Ty with Newton-Rapson Method
        xn = 0.5
        Tn = T_med
        for _iter in range(max_iter):
            cosmosac_class.x = [xn, 1 - xn]
            cosmosac_class.T = Tn
            ac12, ac21 = cosmosac_class.gam()

            f_xn = yi * P - xn * ac12 * vp.Cal_Vapor_Pressure(CASRN_1, Tn)
            g_xn = (1 - yi) * P - (1 - xn) * ac21 * vp.Cal_Vapor_Pressure(CASRN_2, Tn)

            if abs(f_xn) + abs(g_xn) < tol:
                Ty.append(O_xn[1])
                break

            P_xn = np.array([xn, Tn]).T
            F_xn = np.array([f_xn, g_xn]).T

            cosmosac_class.x = [xn + tol, 1 - xn - tol]
            cosmosac_class.T = Tn
            ac12_dx, ac21_dx = cosmosac_class.gam()
            cosmosac_class.x = [xn, 1 - xn]
            cosmosac_class.T = Tn + tol
            ac12_dT, ac21_dT = cosmosac_class.gam()

            df_dx = (
                yi * P
                - (xn + tol) * ac12_dx * vp.Cal_Vapor_Pressure(CASRN_1, Tn)
                - f_xn
            ) / tol
            dg_dx = (
                (1 - yi) * P
                - (1 - xn - tol) * ac21_dx * vp.Cal_Vapor_Pressure(CASRN_2, Tn)
                - g_xn
            ) / tol
            df_dT = (
                yi * P - xn * ac12_dT * vp.Cal_Vapor_Pressure(CASRN_1, Tn + tol) - f_xn
            ) / tol
            dg_dT = (
                (1 - yi) * P
                - (1 - xn) * ac21_dT * vp.Cal_Vapor_Pressure(CASRN_2, Tn + tol)
                - g_xn
            ) / tol
            J = np.array(
                [
                    [df_dx, df_dT],
                    [dg_dx, dg_dT],
                ],
            )

            J_inv = np.linalg.pinv(J)
            O_xn = P_xn - np.matmul(J_inv, F_xn)
            xn = O_xn[0]
            Tn = O_xn[1]

            if xn < 0:
                xn = 0.01
            elif xn > 1:
                xn = 0.99

    print(gamma_list)
    plt.plot([0, *x1, 1], [Tb_2, *Tx, Tb_1])
    plt.plot([0, *x1, 1], [Tb_2, *Ty, Tb_1])
    plt.xlim(0, 1)
    plt.ylim(min([Tb_2, *Tx, *Ty, Tb_1]) - 2, max([Tb_2, *Tx, *Ty, Tb_1]) + 2)
    plt.show()


if __name__ == "__main__":
    main()
