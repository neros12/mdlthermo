import numpy as np
from typing import List, TypedDict

_q0 = 79.53  # area normalization parameter [Å**2]
_r0 = 66.69  # volume normalization parameter [Å**3]
_z = 10  # coordination number
_R = 1.987204258e-3  # gas constant [kcal/K/mol]
_aeff = 7.25  # effective area [Å**2], number of sigma profiles,


def _cal_DW(T, version):
    """
    Calculate the exchange energy.

    The exchange energy has the values for each charge density combinations
    and sigma profile type combinations, therefore having the shape of
    (num_sp, num_sp, 51, 51).

    Parameters
    ----------
    T : float
        The system temperature.

    Returns
    -------
    DW : numpy.ndarray of shape=(num_sp, num_sp, 51, 51)
        The exchange energy.
    """
    if version == 2020:
        num_sp = 4
        _chb = np.array(
            [
                [0, 0, 0, 0],
                [0, 4013.78, 3016.43, 3020.18],
                [0, 3016.43, 932.31, 1872.84],
                [0, 3020.18, 1872.84, 2225.67],
            ]
        )  # hydrogen bonding parameter [kcal*Å^4/mol/e^2],
    else:
        num_sp = 3
        _chb = np.array(
            [
                [0, 0, 0],
                [0, 4013.78, 3016.43],
                [0, 3016.43, 932.31],
            ]
        )  # hydrogen bonding parameter [kcal*Å^4/mol/e^2],

    _AES = 6525.69  # electrostatic constant A [kcal*ang**4/mol/e**2]
    _BES = 1.4859e8  # electrostatic constant B [kcal*Å**4*K**2/mol/e**2]
    _cES = _AES + _BES / T / T  # electrostatic parameter [kcal*Å^4/mol/e^2]

    # Initialize parameters
    sig = np.linspace(-0.025, 0.025, 51)
    sigT = sig.reshape(-1, 1)
    DW = np.zeros((num_sp, num_sp, 51, 51))

    # Calculate exchange energy for each pair of sigma profile types
    for i in range(num_sp):
        for j in range(i + 1):
            mask = (sig * sigT) < 0
            chb_part = np.where(mask, _chb[i, j] * (sig - sigT) ** 2, 0)

            # Calculate total exchange energy
            DW[i, j] = DW[j, i] = _cES * (sig + sigT) ** 2 - chb_part

    return DW


def cal_ln_gam_comb(A, V, x):
    """
    Calculate log of combinatory activity coefficients.

    Parameters
    ----------
    None.

    Returns
    -------
    ln_gam_comb : numpy.ndarray of shape=(num_comp,)
        Combinatory activity coefficients of components.
    """
    # calculate normalized areas and volumes
    q = A / _q0
    r = V / _r0
    L = (_z / 2) * (r - q) - (r - 1)

    theta = q / np.sum(x * q)
    phi = r / np.sum(x * r)

    # calcualte combinatory activity coefficients
    ln_gam_comb = (
        np.log(phi) + _z * q * np.log(theta / phi) / 2 + L - phi * np.sum(x * L)
    )

    return ln_gam_comb


def cal_ln_gam_res(A, psigA, x, T, version):
    """
    Calculate residual activity coefficients.

    Parameters
    ----------
    None.

    Returns
    -------
    ln_gam_res : numpy.ndarray of shape=(num_comp,)
        Residual activity coefficients of components.
    """
    # calculate intermediate terms
    psig = np.einsum("itm,i->itm", psigA, 1 / A)
    psig_mix = np.einsum("i,itm->tm", x, psigA) / np.sum(x * A)

    exp_DW = np.exp(-_cal_DW(T, version) / _R / T)

    A_plus = np.einsum("stmn,isn->istmn", exp_DW, psig)  # A^(+)
    A_plus_mix = np.einsum("stmn,sn->stmn", exp_DW, psig_mix)  # A^(+)_mix

    # calculate the segment activity coefficients
    Gam = np.ones(np.shape(psig))
    Gam_mix = np.ones(np.shape(psig_mix))
    diff = 1

    for _ in range(500):
        Gam_old = np.array(Gam)
        Gam_mix_old = np.array(Gam_mix)

        # Update Gam element-wise
        for i in range(Gam.shape[0]):
            for t in range(Gam.shape[1]):
                for m in range(Gam.shape[2]):
                    Gam[i, t, m] = 1 / np.einsum(
                        "sn,sn->", A_plus[i, :, t, m, :], Gam[i, :, :]
                    )

        # Update Gam_mix element-wise
        for t in range(Gam_mix.shape[0]):
            for m in range(Gam_mix.shape[1]):
                Gam_mix[t, m] = 1 / np.einsum(
                    "sn,sn->", A_plus_mix[:, t, m, :], Gam_mix[:, :]
                )

        # check convergence
        diff = np.sum((Gam - Gam_old) ** 2)
        diff_mix = np.sum((Gam_mix - Gam_mix_old) ** 2)

        if diff <= 1e-6 and diff_mix <= 1e-6:
            break
    else:
        raise Exception("Converge failed")

    # calculate residual activity coefficients
    Gam_part = np.log(Gam_mix) - np.log(Gam)
    ln_gam_res = np.einsum("itm,itm->i", psigA, Gam_part) / _aeff

    return ln_gam_res


def cal_ln_gam_dsp(x, ek, dnatr):
    """
    Calculate dispersive activity coefficients.

    Parameters
    ----------
    None.

    Returns
    -------
    ln_gam_dsp : numpy.ndarray of shape=(num_comp,)
        Dispersive activity coefficients of components.
    """
    num_mol = len(x)
    if None in ek or np.isnan(ek).any():

        return np.zeros(num_mol)

    ekT = ek.reshape(-1, 1)

    # check if dispersion activity coefficients are applicable
    if None in ek or None in dnatr:
        ln_gam_dsp = np.array([0] * num_mol)

        return ln_gam_dsp

    # calculate interaction parameters
    w = np.ones((num_mol, num_mol)) * 0.27027
    wpair = [
        {"WATER", "HBOA"},
        {"COOH", "NHB"},
        {"COOH", "HBDA"},
        {"WATER", "COOH"},
    ]
    for i in range(num_mol):
        for j in range(i):
            if {dnatr[i], dnatr[j]} in wpair:
                w[i][j] = w[j][i] = -0.27027

    A = w * (0.5 * (ek + ekT) - np.sqrt(ek * ekT))  # not area

    # calculate dispersive activity coefficients
    ln_gam_dsp = np.zeros(num_mol)
    for i in range(num_mol):
        for j in range(num_mol):
            if i != j:
                ln_gam_dsp[i] = ln_gam_dsp[i] + x[j] * A[i, j]
            if j > i:
                ln_gam_dsp[i] = ln_gam_dsp[i] - x[i] * x[j] * A[i, j]

    return ln_gam_dsp


class ChemicalProfile(TypedDict):
    area: float
    volume: float
    sigma_profiles: np.ndarray


class ChemicalProfiles(TypedDict):
    version: int
    data: List[ChemicalProfile]


def calculate_gamma(
    chemical_profiles: ChemicalProfiles,
    x: List[float],
    T: float,
) -> List[float]:
    """
    Calculate COSMO-SAC activity coefficients for a mixture of components.

    This function computes the activity coefficients of components in a mixture
    based on the COSMO-SAC model. It considers combinatorial, residual, and dispersion
    contributions using sigma profiles and molecular properties.

    Parameters
    ----------
    chemical_profiles : ChemicalProfiles
        Dictionary containing 'version' and a list of component profiles. Each profile
        includes area, volume, sigma_profiles, ek (electrostatic energy), and natr
        (atom type information).
    x : List[float]
        Mole fraction of each component in the mixture.
    T : float
        Temperature in Kelvin.

    Returns
    -------
    gam : List[float]
        A list of activity coefficients for each component in the same order as `x`.
    """
    version = chemical_profiles["version"]
    if version == 2020:
        num_sp = 4
    else:
        num_sp = 3

    areas = np.array([])
    volumes = np.array([])
    psigA = np.array([]).reshape(0, num_sp, 51)
    eks = np.array([])
    natrs = []
    for chemical_profile in chemical_profiles["data"]:
        areas = np.append(areas, chemical_profile["area"])
        volumes = np.append(volumes, chemical_profile["volume"])
        psigA = np.vstack((psigA, chemical_profile["sigma_profiles"]))
        eks = np.append(eks, chemical_profile["ek"])
        natrs.append(chemical_profile["natr"])

    ln_gam_comb = cal_ln_gam_comb(areas, volumes, x)
    ln_gam_res = cal_ln_gam_res(areas, psigA, x, T, version)
    ln_gam_dsp = cal_ln_gam_dsp(x, eks, natrs)

    ln_gam = ln_gam_comb + ln_gam_res + ln_gam_dsp
    gam: np.ndarray = np.exp(ln_gam)

    return gam.tolist()
