import numpy as np
import matplotlib.pyplot as plt


# component 0 = Naphthalene
# component 1 = 1-Chloro-4-nitrobenzene
# Properties
Tm0 = 353.91
Tm1 = 356.77
Hfus0 = 18.811 * 1000
Hfus1 = 15.03 * 1000
R = 8.314


def NRTL(x1: float, T: float) -> dict[int, float]:
    # Binary NRTL
    num_components = 2
    x = {0: x1, 1: 1 - x1}

    # NRTL Parameters
    a = {(0, 1): 0.425142, (1, 0): -0.471558}
    b = {(0, 1): -0.00746743, (1, 0): -0.312964}
    c = 0.5

    tow = {
        (0, 0): 0,
        (0, 1): a[(0, 1)] + b[(0, 1)] / T,
        (1, 0): a[(1, 0)] + b[(1, 0)] / T,
        (1, 1): 0,
    }
    G = {
        (0, 0): 1,
        (0, 1): np.exp(-c * tow[(0, 1)]),
        (1, 0): np.exp(-c * tow[(1, 0)]),
        (1, 1): 1,
    }

    gamma = {}
    for i in range(num_components):
        par_1 = 0
        for j in range(num_components):
            par_1 += x[j] * tow[(j, i)] * G[(j, i)]

        par_2 = 0
        for k in range(num_components):
            par_2 += x[k] * G[(k, i)]

        par_3 = 0
        for j in range(num_components):
            par_4 = 0
            for j in range(num_components):
                par_4 += x[k] * G[(k, j)]

            par_5 = 0
            for m in range(num_components):
                par_5 += x[m] * tow[(m, j)] * G[(m, j)]

            par_6 = 0
            for k in range(num_components):
                par_6 += x[k] * G[(k, j)]

            par_3 += (x[j] * G[(i, j)]) / par_4 * (tow[(i, j)] - par_5 / par_6)

        gamma[i] = np.exp(par_1 / par_2 + par_3)

    return gamma


cal_x0 = np.arange(0.01, 1, 0.01)
cal_T0 = []
cal_T1 = []

# Bisection Method
Tmin = 100
Tmax = 400
tol = 1e-10
for _x in cal_x0:
    for iter in range(5000):
        if iter == 0:
            Tl = Tmin
            Tu = Tmax

        Tr = (Tl + Tu) / 2

        f_Tl = np.log(_x * NRTL(_x, Tl)[0]) + Hfus0 / R * (1 / Tl - 1 / Tm0)
        f_Tr = np.log(_x * NRTL(_x, Tr)[0]) + Hfus0 / R * (1 / Tr - 1 / Tm0)
        f_Tu = np.log(_x * NRTL(_x, Tu)[0]) + Hfus0 / R * (1 / Tu - 1 / Tm0)

        if abs(f_Tl * f_Tr) < tol:
            break
        elif f_Tl * f_Tr < 0:
            Tu = Tr
        elif f_Tl * f_Tr > 0:
            Tl = Tr

    cal_T0.append(Tr)

for _x in cal_x0:
    for iter in range(5000):
        if iter == 0:
            Tl = Tmin
            Tu = Tmax

        Tr = (Tl + Tu) / 2

        f_Tl = np.log((1 - _x) * NRTL((1 - _x), Tl)[1]) + Hfus1 / R * (1 / Tl - 1 / Tm0)
        f_Tr = np.log((1 - _x) * NRTL((1 - _x), Tr)[1]) + Hfus1 / R * (1 / Tr - 1 / Tm0)
        f_Tu = np.log((1 - _x) * NRTL((1 - _x), Tu)[1]) + Hfus1 / R * (1 / Tu - 1 / Tm0)

        if abs(f_Tl * f_Tr) < tol:
            break
        elif f_Tl * f_Tr < 0:
            Tu = Tr
        elif f_Tl * f_Tr > 0:
            Tl = Tr

    cal_T1.append(Tr)

# Experimental Results
exp_x = 1 - np.array([0, 0.1664, 0.2447, 0.3705, 0.5397, 0.6288, 0.7551, 0.8334, 1])
exp_T = np.array(
    [353.91, 341.93, 336.47, 324.65, 317.39, 327.07, 338.87, 345.18, 356.77]
)

cal_T = []
for i in range(len(cal_T0)):
    cal_T.append(max(cal_T0[i], cal_T1[i]))
cal_T = np.array(cal_T)

plt.scatter(exp_x, exp_T, s=10, c="dimgray", zorder=2)
plt.plot(cal_x0, cal_T, c="lightgray", zorder=1)
plt.xlim(0, 1)
plt.ylim(int(min(cal_T)), int(max(cal_T) + 10))
plt.tick_params(axis="y", direction="in")
plt.tick_params(axis="x", direction="in")
plt.xlabel("Naphthalene  mol%")
plt.ylabel("Temperature  (K)")
plt.show()
