from numpy import array, tile, exp, zeros, matmul, ndarray
from numpy import sum as npsum


def NRTL(X: ndarray, A: ndarray, B: ndarray, T: float):
    """
    Nc = Number of Component

    X : shape(Nc)
    A : shape(Nc, Nc)
    B : shape(Nc, Nc)
    """

    Nc = len(X)
    alpha = 0.3

    tau = A + B / T
    G = exp(-alpha * tau)

    temp0 = matmul(X, tau * G) / matmul(X, G)
    temp1 = tile(X.reshape(Nc, 1), (1, Nc)) * G.T
    temp2 = tile(npsum(tile(X, (Nc, 1)) * G.T, axis=1, keepdims=True), (1, Nc))
    temp3 = tile(npsum(tile(X, (Nc, 1)) * tau.T * G.T, axis=1, keepdims=True), (1, Nc))
    temp4 = tile(npsum(tile(X, (Nc, 1)) * G.T, axis=1, keepdims=True), (1, Nc))

    lnGamma = temp0 + npsum(temp1 / temp2 * (tau.T - temp3 / temp4), 0)

    return exp(lnGamma)


def Antonine(A, B, C, T):

    return exp(A - B / (T + C))


def FittingModel(x_data, a12, a21, b12, b21, A1, B1, C1, A2, B2, C2):
    x1, T = x_data
    x2 = 1 - x1

    r1, r2 = NRTL(
        array([x1, x2]),
        array([[0, a12], [a21, 0]]),
        array([[0, b12], [b21, 0]]),
        T,
    )
    Psat1 = Antonine(A1, B1, C1, T)
    Psat2 = Antonine(A2, B2, C2, T)

    return x1 * r1 * Psat1 + x2 * r2 * Psat2
