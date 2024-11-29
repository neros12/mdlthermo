# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:15:45 2024.

@author: Beom Chan Ryu

COSMO-SAC ML
"""
import numpy as np
import pickle
from cosmosac import CosmoSac
from rdkit import Chem
from scipy.linalg import fractional_matrix_power


list_smarts = [
    "[CX4H3]",
    "[CX3H2v4]",
    "[CX2H1v4]",
    "[!R;CX4H2]",
    "[!R;CX4H]",
    "[!R;CX4H0]",
    "[!R;CX3H1v4]",
    "[!R;CX3H0v4]",
    "[!R;CX2H0;$([CX2H0](=*)(=*))]",
    "[!R;CX2H0;$([CX2H0](#*))]",
    "[R;CX4H2]",
    "[R;CX4H]",
    "[R;CX4H0]",
    "[R;CX3H1v4]",
    "[R;CX3H0v4]",
    "[R;CX2H0;$([CX2H0](=*)(=*))]",
    "[R;CX2H0;$([CX2H0](#*))]",
    "[R;cX4h2]",
    "[R;cX4h]",
    "[R;cX3h1v4]",
    "[R;cX3h0v4]",
    "[FX1H0]",
    "[ClX1H0]",
    "[BrX1H0]",
    "[IX1H0]",
    "[OX2H1]",
    "[!R;OX2H0]",
    "[R;OX2H0]",
    "[R;oX2h0]",
    "[OX1H0v2]",
    "[NX3H2v3]",
    "[NX3H1v3;!R]",
    "[NX3H1v3;R]",
    "[nX3h1v3;R]",
    "[NX3H0v3;!R]",
    "[NX3H0v3;R]",
    "[nX3h0v3;R]",
    "[NX2H0v3;!R]",
    "[NX2H0v3;R]",
    "[nX2h0v3;R]",
    "[NX1H0v3]",
    "[SX2H1v2]",
    "[SX2H0v2;!R]",
    "[SX2H0v2;R]",
    "[sX2h0v2;R]",
    "[SX1H0v2]",
    "[OX1H0v1]",
    "[NX3H0v4]",
]

with open("sig_param.pickle", "rb") as f:
    sig_param = pickle.load(f)
with open("nhb_param.pickle", "rb") as f:
    nhb_param = pickle.load(f)
with open("oh_param.pickle", "rb") as f:
    oh_param = pickle.load(f)
with open("ot_param.pickle", "rb") as f:
    ot_param = pickle.load(f)
with open("vol_param.pickle", "rb") as f:
    vol_param = pickle.load(f)


def get_gcgcn_input(smiles):
    """Generate node and edge feature matrices for GC-GCN.

    Parameters
    ----------
    smiles : str
        SMILES string of the compound.

    Returns
    -------
    nfm, efm : tuple of numpy.ndarray
        Node feature matrix (nfm) of shape (25, len(list_smarts)) and edge
        feature matrix (efm) of shape (25, 25).
    """
    mol = Chem.MolFromSmiles(smiles)

    # Get node feature matrix
    nfm = np.zeros((25, len(list_smarts)))

    for smarts_index, smarts in enumerate(list_smarts):
        pat = Chem.MolFromSmarts(smarts)
        nfm[mol.GetSubstructMatches(pat), smarts_index] = 1

    # Get edge feature matrix
    efm = Chem.GetAdjacencyMatrix(mol).astype("float64")
    np.fill_diagonal(efm, 1)

    diag = np.diag(np.sum(efm, axis=1))
    diag_half = fractional_matrix_power(diag, -0.5)
    efm = np.matmul(np.matmul(diag_half, efm), diag_half)

    # Padding edge feature matrix
    n_heavyatom = len(efm)
    pad = 25 - n_heavyatom
    efm = np.pad(efm, ((0, pad), (0, pad)), "constant", constant_values=0.0)

    return nfm, efm


def pred_gcgcn(efm, nfm, param_list):
    """Predict output using GC-GCN model with learned parameters.

    Parameters
    ----------
    efm : np.ndarray
        Edge feature matrix of shape (25, 25).
    nfm : np.ndarray
        Node feature matrix of shape (25, len(list_smarts)).
    param_list : list of np.ndarray
        List of parameters for the GC-GCN model layers.

    Returns
    -------
    x : numpy.ndarray
        Predicted values as a flat array.
    """
    x = np.dot(efm, nfm)
    x = np.dot(x, param_list[0]) + param_list[1]  # Graph convolution layer
    x = np.where(x > 0, x, 0)  # ReLU

    x = np.dot(efm, x)
    x = np.dot(x, param_list[2]) + param_list[3]
    x = np.where(x > 0, x, 0)

    x = x.reshape(-1)
    x = np.dot(x, np.tile(np.eye(256), (25, 1)))  # Node-wise summation

    x = np.dot(x, param_list[4]) + param_list[5]  # Dense layer
    x = np.where(x > 0, x, 0)

    x = np.dot(x, param_list[6]) + param_list[7]
    x = np.where(x > 0, x, 0)

    x = np.dot(x, param_list[8]) + param_list[9]
    x = np.where(x > 0, x, 0)

    x = np.dot(x, param_list[10]) + param_list[11]
    return x


class CosmoSacGcgcn(CosmoSac):
    """COSMO-SAC model using GC-GCN for sigma and volume prediction."""

    def add_comp(self, smiles, name=None):
        """Add a component with GC-GCN predictions to COSMO-SAC.

        Parameters
        ----------
        smiles : str
            SMILES string of the component.
        name : str, optional
            Name of the component.

        Returns
        -------
        None
        """
        nfm, efm = get_gcgcn_input(smiles)

        volume = 562 * pred_gcgcn(efm, nfm, vol_param)[0]

        sigma_profiles = np.zeros((3, 51))
        sigma_profiles[0] = 145 * pred_gcgcn(efm, nfm, nhb_param)
        sigma_profiles[1] = 7 * pred_gcgcn(efm, nfm, oh_param)
        sigma_profiles[2] = 16 * pred_gcgcn(efm, nfm, ot_param)
        sigma_profiles = np.where(sigma_profiles < 0, 0, sigma_profiles)

        # area = np.sum(145*pred_gcgcn(efm, nfm, sig_param))
        area = np.sum(sigma_profiles)

        self.A.append(area)
        self.V.append(volume)
        self.psigA.append(sigma_profiles)
        self.name.append(name)

    def get_totsig(self, smiles):
        """Compute total sigma profile for a compound.

        Parameters
        ----------
        smiles : str
            SMILES string of the compound.

        Returns
        -------
        psigA : np.ndarray
            Predicted sigma profile as a 1D array.
        """
        nfm, efm = get_gcgcn_input(smiles)

        psigA = 145 * pred_gcgcn(efm, nfm, sig_param)
        psigA = np.where(psigA < 0, 0, psigA)

        return psigA

    def gam(self):
        """Calculate COSMO-SAC activity coefficient (gamma).

        Returns
        -------
        gam : np.ndarray
            Activity coefficients for each component.
        """
        ln_gam_comb = self.ln_gam_comb()
        ln_gam_res = self.ln_gam_res()

        ln_gam = ln_gam_comb + ln_gam_res
        gam = np.exp(ln_gam)
        return gam


if __name__ == "__main__":
    # Define COSMO-SAC ML model
    model = CosmoSacGcgcn(version=2010)

    # Add components
    model.add_comp("CCO", name="ethanol")
    model.add_comp("c1ccccc1", name="benzene")

    # Add system conditions
    model.x = [0.5, 0.5]
    model.T = 298.15

    # Calculate activity coefficients
    gamma = model.gam()

    # Print results
    print(f"Components: {model.name}")
    print(f"Mole fractions: {model.x}")
    print(f"System temperature: {model.T} K")
    print(f"Activity coefficients: {gamma}")
