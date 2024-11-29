import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from scipy.stats import t as t_distribution


def Z_Matrix(
    x: np.ndarray, par: np.ndarray, model, alpha: float, number_of_data: int
) -> np.ndarray:
    result = np.zeros((number_of_data, len(par)))
    for i in range(len(par)):
        delta = np.array(par)
        delta[i] += alpha
        result[:, i] = (model(x, delta) - model(x, par)) / alpha

    return result


def Square_Error(
    x: np.ndarray[float], y: np.ndarray[float], par: np.ndarray[float], model
) -> float:
    return sum((y - model(x, par)) ** 2)


def Gauss_Newton_Method(
    data_x: np.ndarray[float],
    data_y: np.ndarray[float],
    fitting_model,
    initial_parameters: np.ndarray[float],
    number_of_data: int,
    cost_function=None,
    custom_cost_function=False,
    max_iteration=1000,
    tolerance=1e-5,
    alpha=1e-5,
) -> np.ndarray:
    """
    fitting model must type must be like\n
    ####################################


    def fitting_model(
        data_x: ndarray[float],
        parameters: ndarray[float],
    )-> float:

        return __result__

    ###################################

    returns:
    --------
        parameter : ndarray[float], uncertainty: ndarray[float]
    """

    data_x = np.array(data_x, dtype="float64")
    data_y = np.array(data_y, dtype="float64")
    parameters = np.array(initial_parameters, dtype="float64")
    degree_of_freedom = number_of_data - len(parameters)
    old_square_error = 0
    new_square_error = 0
    for _iteration in range(max_iteration + 1):
        # Check if it failed at convergence
        if _iteration == max_iteration:
            raise ValueError("Failed at Convergence ")

        Z = Z_Matrix(data_x, parameters, fitting_model, alpha, number_of_data)
        ZZ = np.matmul(Z.T, Z)
        inv_ZZ = np.linalg.pinv(ZZ)
        D = data_y - fitting_model(data_x, parameters)
        ZD = np.matmul(Z.T, D)
        parameters += np.matmul(inv_ZZ, ZD)

        # Check if it reached to object
        if custom_cost_function:
            new_square_error = cost_function(data_x, data_y, parameters, fitting_model)
        else:
            new_square_error = Square_Error(data_x, data_y, parameters, fitting_model)

        if abs(old_square_error - new_square_error) / old_square_error < tolerance:
            variance = new_square_error / degree_of_freedom
            std_par = np.sqrt(np.diagonal(inv_ZZ) * variance)
            t_confidence = t_distribution.interval(0.95, degree_of_freedom)[0]
            uncertainty = abs(std_par * t_confidence)

            break
        else:
            old_square_error = new_square_error

    return parameters, uncertainty


# np.set_printoptions(linewidth=np.inf)
# df = pd.read_csv("DDB_Binary_VLE.csv")


def test_function2(x, p):

    return (
        p[0] * x**6
        + p[1] * x**5
        + p[2] * x**4
        + p[3] * x**3
        + p[4] * x**2
        + p[5] * x
        + p[6]
    )


def test_function(x, a, b, c, d, e, f, g):

    return a * x**6 + b * x**5 + c * x**4 + d * x**3 + e * x**2 + f * x + g


def is_well_fitted(par, fun, X, Y):
    relatve_error = abs(Y - fun(X, *par)) / abs(Y)
    std_relative_error = np.std(relatve_error)

    if std_relative_error > 0.008:
        return False
    else:
        return True


def is_reliable_data(X, Y):
    try:
        popt, pcov = curve_fit(test_function, X, Y)

        return is_well_fitted(popt, test_function, X, Y)
    except:
        return False


# def find_data(df: pd.DataFrame, component1: str, component2: str) -> pd.DataFrame:
#     df1 = df[df["cmp1_name"] == component1]
#     df1 = df1[df1["cmp2_name"] == component2]
#     # df2 = df[df["cmp1_name"] == component2]
#     # df2 = df2[df2["cmp2_name"] == component1]

#     # return pd.concat([df1, df2])

#     return df1


# dict_data = {}
# data_set = 1
# list_data_set_ID = []
# for row in df.itertuples():
#     cmp1_name = row.cmp1_name
#     cmp2_name = row.cmp2_name
#     dTexp = row.dTexp
#     dPexp = row.dPexp
#     dX1exp = row.dX1exp
#     dY1exp = row.dY1exp

#     data_key = (cmp1_name, cmp2_name, dPexp)
#     if data_key not in dict_data:
#         list_data_set_ID.append(data_set)
#         dict_data[data_key] = data_set
#         data_set += 1
#     else:
#         list_data_set_ID.append(dict_data[data_key])

# df["data_set_ID"] = list_data_set_ID


# df.to_csv("DDB_Binary_VLE02.csv", index=False)
# df.to_excel("DDB_Binary_VLE02.xlsx", index=False)


# df = pd.read_csv("DDB_Binary_VLE02.csv")

# num_data_set_id = 183690
# dict_is_reliable_data = {}
# for set_id in range(1, num_data_set_id + 1):
#     indexed_df = df[df["data_set_ID"] == set_id]
#     if len(indexed_df) > 9:
#         x1 = indexed_df["dX1exp"].to_numpy()
#         y1 = indexed_df["dY1exp"].to_numpy()
#         temperature = indexed_df["dTexp"].to_numpy()

#         if is_reliable_data(x1, temperature) and is_reliable_data(y1, temperature):
#             dict_is_reliable_data[set_id] = True
#         else:
#             dict_is_reliable_data[set_id] = False
#     else:
#         dict_is_reliable_data[set_id] = True

# list_is_reliable = []
# for row in df.itertuples():
#     list_is_reliable.append(dict_is_reliable_data[row.data_set_ID])

# df["is_reliable"] = list_is_reliable

# df.to_csv("DDB_Binary_VLE03.csv", index=False)
# df.to_excel("DDB_Binary_VLE03.xlsx", index=False)


df = pd.read_csv("DDB_Binary_VLE04.csv")
df_not_reliable = df[df["is_reliable"] == False]

list_of_ids = list(set(df_not_reliable["data_set_ID"].to_list()))
list_of_ids.sort()

for id in list_of_ids:
    sample_df = df_not_reliable[df_not_reliable["data_set_ID"] == id]

    X_data = sample_df[["dX1exp", "dTexp"]].values.tolist()
    Y_data = sample_df[["dY1exp", "dTexp"]].values.tolist()
    X_data.sort(key=lambda x: x[0])
    Y_data.sort(key=lambda x: x[0])
    X_data = np.array(X_data)
    Y_data = np.array(Y_data)

    plt.figure(id)
    plt.scatter(X_data[:, 0], X_data[:, 1], s=5)
    plt.scatter(Y_data[:, 0], Y_data[:, 1], s=5)
    plt.show()

# df = pd.read_csv("DDB_Binary_VLE04.csv")
# df_not_reliable = df[df["is_reliable"] == False]
# list_of_ids = set(df_not_reliable["data_set_ID"].to_list())


# dict_is_reliable_data = {}
# for set_id in list_of_ids:
#     indexed_df = df[df["data_set_ID"] == set_id]
#     x1 = indexed_df["dX1exp"].to_numpy()
#     y1 = indexed_df["dY1exp"].to_numpy()
#     temperature = indexed_df["dTexp"].to_numpy()

#     if is_reliable_data(x1, temperature) and is_reliable_data(y1, temperature):
#         dict_is_reliable_data[set_id] = True
#     else:
#         dict_is_reliable_data[set_id] = False

# for row in df.itertuples():
#     print(row.Index)
#     if not row.is_reliable:
#         if dict_is_reliable_data[row.data_set_ID]:
#             df.loc[row.Index, "is_reliable"] = True

# df.to_excel("DDB_Binary_VLE04.xlsx", index=False)
