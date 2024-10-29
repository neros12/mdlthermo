import os
import warnings
import logging
from os.path import join as opj

# 로깅 수준 설정
# warnings.filterwarnings("ignore")
# os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import keras
import numpy as np
import tensorflow as tf
from pathlib import Path
from rdkit import Chem
from keras import initializers
from keras.layers import Input, Dense, LayerNormalization, LeakyReLU
from keras.models import Model, load_model
from scipy.linalg import fractional_matrix_power

# tf.get_logger().setLevel(logging.ERROR)
DIR_PATH = Path(__file__).parent

_descript_group = [
    "[CX4v4]",
    "[CX3v4;$([CX3v4](=*))]",
    "[CX2v4;$([CX2v4](=*)(=*))]",
    "[CX2v4;$([CX2v4](#*))]",
    "[c]",
    "[NX3v3]",
    "[NX2v3;$([NX2v3](=*))]",
    "[NX1v3;$([NX1v3](#*))]",
    "[N+v4]",
    "[N-v2]",
    "[n+0]",
    "[OX2v2]",
    "[OX1v2;$([OX1v2](=*))]",
    "[O+v3]",
    "[O-v1]",
    "[o]",
    "[Sv2]",
    "[Sv4]",
    "[Sv6]",
    "[s]",
    "[Pv3]",
    "[Pv5]",
    "[p]",
    "[FX1v1]",
    "[ClX1v1]",
    "[BrX1v1]",
    "[IX1v1]",
]

_num_feat = len(_descript_group) + 6  # 5 Hs and aliphatic rings


class ReduceSumLayer(keras.layers.Layer):
    def call(self, x):
        return tf.reduce_sum(x, axis=1)


class Antoine(keras.layers.Layer):
    """
    Custom layer to perform the Antoine equation.

    ln(p) = A - B/(T + C)

    Parameters
    ----------
    None
    """

    def __init__(self, **kwargs):
        super(Antoine, self).__init__(**kwargs)
        self.num_param = 3

    def call(self, inputs):
        """
        Perform the Antoine equation calculation.

        Parameters
        ----------
        inputs: List of input tensors [x, t]
            - x: Input tensor of parameters with shape (?, 3).
                The parameters are A, B, and C.
            - t: Input tensor of temperature with shape (?, 1).

        Returns
        -------
        result : tensorflow.Tensor with shape (?, 1)
            Output tensor of log vapor pressure with shape (?, 1).
        """
        x, t = inputs
        a, b, c = tf.split(x, num_or_size_splits=3, axis=1)  # shape=(?, 1)

        result = a - b / (t + c)
        return result


class Wagner(keras.layers.Layer):
    """
    Custom layer to perform the Wagner equation.

    ln(pr) = (a*tau + b*tau**1.5 + c*tau**2.5 + d*tau**5)/tau
    tau = 1 - Tr

    Note that the Tc and ln(pc) is acting as parameters, including a, b, c,
    and d. Also, tc is added by 1500, to prevent tau from being less then zero.

    Parameters
    ----------
    None
    """

    def __init__(self, **kwargs):
        super(Wagner, self).__init__(**kwargs)
        self.num_param = 6

    def call(self, inputs):
        """
        Perform the Wagner equation calculation.

        Parameters
        ----------
        inputs: List of input tensors [x, t]
            - x: Input tensor of parameters with shape (?, 4).
                Parameters are a, b, c, and d.
            - t: Input tensor of temperature with shape (?, 1).

        Returns
        -------
        result : tensorflow.Tensor with shape (?, 1)
            Output tensor of log vapor pressure with shape (?, 1).
        """
        x, t = inputs
        a, b, c, d, tcs, lnpc = tf.split(x, num_or_size_splits=6, axis=1)
        # tcs : scaled critical temperature [kK]
        tc = tcs + 1.5  # 충분히 큰 값으로 tau 값을 0 이상으로 한다

        tau = tf.maximum(0.0, 1 - 0.001 * t / tc)  # T가 Tc 이상이어도 tau는 0
        result = (a * tau + b * tau**1.5 + c * tau**2.5 + d * tau**5) / (1 - tau) + lnpc
        return result


class KingAlNajjar(keras.layers.Layer):
    """
    Custom layer to perform the King and Al-Najaar's equation.

    d/dT(T**2*d/dT(ln(p))) = Tr*(3*b/4/tau**0.5 + 3.75*c*tau**0.5 +
                                 20*d*tau**3)
    tau = 1 - T/Tc
    This custom layer calculates ln(p).

    Parameters
    ----------
    None
    """

    def __init__(self, **kwargs):
        super(KingAlNajjar, self).__init__(**kwargs)
        self.num_param = 6

    def call(self, inputs):
        """
        Perform the King-Al-Najaar equation calculation.

        Parameters
        ----------
        inputs: List of input tensors [x, tau]
            - x: Input tensor of parameters with shape (?, 4).
                Parameters are b, c, d, and e.
            - t: Input tensor of temperature with shape (?, 1).

        Returns
        -------
        result : tensorflow.Tensor with shape (?, 1)
            Output tensor of log vapor pressure with shape (?, 1).
        """
        x, t = inputs
        b, c, d, e, tcs, lnpc = tf.split(x, num_or_size_splits=6, axis=1)
        # tcs : scaled critical temperature
        tc = tcs + 1.5  # [kK] 충분히 큰 값으로 tau 값을 0 이상으로 한다.

        tau = tf.maximum(0.0, 1 - 0.001 * t / tc)  # T가 Tc 이상이어도 tau는 0
        result = (
            (-b - c) * tau**0.5
            - d * tau
            - c * tau**1.5
            - d * (tau**2 + tau**4 + tau**6)
            + (-2 * b + c + d + e) * tau**0.5 / (1 - tau)
            + 3 * b * tf.math.atanh(tau**0.5)
        ) + lnpc
        return result


class GraphConvolution(keras.layers.Layer):
    """
    Graph Convolutional Layer.

    This layer performs graph convolution on the input features using the
    adjacency matrix and weight matrix.

    Parameters
    ----------
    units: int
        Number of output units.
    """

    def __init__(
        self,
        units,
        activation=None,
        bias_initializer="zeros",
        kernel_initializer="glorot_uniform",
        **kwargs,
    ):
        super(GraphConvolution, self).__init__(**kwargs)

        self.units = units
        self.activation = keras.activations.get(activation)
        self.bias_initializer = keras.initializers.get(bias_initializer)
        self.kernel_initializer = keras.initializers.get(kernel_initializer)

    def get_config(self):
        """
        Override the arguments in __init__.

        Parameters
        ----------
        None

        Return
        ------
        congfig : pass
            pass
        """
        config = super().get_config()
        config.update(
            {
                "units": self.units,
                "activation": self.activation,
                "bias_initializer": self.bias_initializer,
                "kernel_initializer": self.kernel_initializer,
            }
        )
        return config

    def build(self, input_shape):
        """
        Build paremeters.

        Parameters
        ----------
        input_shape : int
            pass

        Notes
        -----
        The weight 'w' has the shape of (input_dim, units).
        The bias 'b' has the shape of (units,).
        """
        input_dim = input_shape[0][-1]

        self.w = self.add_weight(
            shape=(input_dim, self.units),
            initializer=self.kernel_initializer,
            trainable=True,
            name="W",
        )
        self.b = self.add_weight(
            shape=(self.units,),
            initializer=self.bias_initializer,
            trainable=True,
            name="b",
        )

    def call(self, inputs):
        """
        Perform Graph Convolution calculation.

        This layer performs graph convolution on the input features using the
        adjacency matrix and weight matrix.

        Parameters
        ----------
        inputs : list of [x, a]
            - x: tensorflow.Tensor
                Input feature tensor of shape
                (batch_size, num_nodes, input_dim).
            - a: tensorflow.Tensor
                Normalized adjacency matrix tensor of shape
                (batch_size, num_nodes, num_nodes).

        Returns
        -------
        h: tensorflow.Tensor
            Output tensor after graph convolution of shape
            (batch_size, num_nodes, units).
        """
        # Assuming inputs shape: (batch_size, num_nodes, input_dim)
        # Assuming adjacency_matrix shape: (batch_size, num_nodes, num_nodes)
        x, a = inputs

        # Perform graph convolution: H^(i+1) = A H^(i) W + b
        h = tf.matmul(a, tf.matmul(x, self.w)) + self.b

        if self.activation is not None:
            h = self.activation(h)

        return h


def is_frag(smiles, max_atom):
    """
    Find if the molecule is fragmentable.

    Parameters
    ----------
    smiles : str
        The SMILES representation of the molecule.
    max_atom : int
        The maximum number of atoms in the molecule.

    Returns
    -------
    bool
        If True, the molecule is fragmentable. Else, it is not fragmentable.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol == None:
        raise ValueError()

    # If the number of heavy atoms in the molecule is larger than max_atom,
    if mol.GetNumAtoms() > max_atom:
        raise ValueError()

    # 분자의 기초 그룹을 찾는다.
    num_find = 0
    for smarts in _descript_group:
        atom_id = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
        atom_id = np.array(atom_id).flatten()

        num_find += atom_id.size

    # 분자의 기초 그룹 합이 전체 heavy atom 수와 안맞으면 fragmentation 실패.
    if not num_find == mol.GetNumHeavyAtoms():
        raise ValueError("fragmentation 실패")

    return True


def get_h(mol, max_atom):
    """
    Get the feature matrix of the molecule.

    The index of the feature matrix corresponds to each group of the group
    contribution method. That is, the feature matrix is an integer encoding in
    which the number of groups constituting a molecule is counted.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The molecule read by rdkit.
    max_atom : int
        The maximum number of atoms in the molecule.

    Returns
    -------
    h : np.ndarray of shape (max_atom, num_feat)
        The (node) feature matrix of the molecule. If the number of heavy atoms
        in the molecule is larger than max_atom, the array of -1 is returned.
    """
    h = np.zeros((max_atom, _num_feat))

    # If the number of heavy atoms in the molecule is larger than max_atom,
    if mol.GetNumAtoms() > max_atom:
        h = h - 1  # Set -1 to the feature matrix.
        return h

    # 분자의 기초 그룹을 찾는다.
    for smarts_id, smarts in enumerate(_descript_group):
        atom_id = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
        atom_id = np.array(atom_id).flatten()

        if atom_id.size > 0:  # If atom_id is not empty,
            h[np.array(atom_id), smarts_id] = 1

    # 분자의 기초 그룹 합이 전체 heavy atom 수와 안맞으면 fragmentation 실패.
    if not np.sum(h) == mol.GetNumHeavyAtoms():
        raise ValueError()

    # 원자의 수소 개수와 고리 여부를 찾는다.
    for atom in mol.GetAtoms():
        atom_id = atom.GetIdx()

        # 수소 개수를 구한다.
        h[atom_id, _num_feat - 7 + atom.GetTotalNumHs()] = 1

        # Aliphatic 고리 구조 여부를 구한다.
        if atom.IsInRing() and not atom.GetIsAromatic():
            h[atom_id, _num_feat - 1] = 1

    return h


def get_a(mol, max_atom):
    """
    Get the normalized adjacency matrix of the molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The molecule read by rdkit.
    max_atom : int
        The maximum number of atoms in the model input.

    Returns
    -------
    a_norm : np.ndarray of shape (max_atom, max_atom)
        Normalized adjacency matrix of the molecule.
    """
    # If the number of heavy atoms in the molecule is larger than max_atom,
    if mol.GetNumAtoms() > max_atom:
        a_norm = np.zeros((max_atom, max_atom)) - 1
        return a_norm

    # Calculate the adjacency matrix A.
    a = Chem.GetAdjacencyMatrix(mol).astype("float64")

    # Fill the diagonal of A with 1 for self-loops.
    np.fill_diagonal(a, 1)

    # Calculate the D^(-0.5) matrix for normalization.
    d = np.diag(np.sum(a, axis=1))
    d_sqrt_inv = fractional_matrix_power(d, -0.5)

    # Compute the normalized adjacency matrix A^(~) = D^(-0.5) A D^(-0.5).
    a_norm = np.matmul(np.matmul(d_sqrt_inv, a), d_sqrt_inv)

    # Pad the matrix with zeros to match the size of max_atom.
    pad = max_atom - len(a)

    a_norm = np.pad(a_norm, ((0, pad), (0, pad)), "constant", constant_values=0.0)
    return a_norm


def get_input(smiles, max_atom):
    """
    Get the matrices of the molecule.

    Parameters
    ----------
    smiles : str
        The SMILES representation of the molecule.
    max_atom : int
        The maximum number of atoms in the molecule.

    Returns
    -------
    Tuple of the matrices (h, a)
        - h : np.ndarray of shape (max_atom, num_feat)
            The (node) feature matrix of the molecule.
        - a : np.ndarray of shape (max_atom, max_atom)
            Normalized adjacency matrix of the molecule.

    See Also
    --------
    get_h, get_a
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError()
    else:
        h = get_h(mol, max_atom)
        a = get_a(mol, max_atom)

    return h, a


def load_models():
    model1 = load_model(
        opj(DIR_PATH, "models", "31.h5"),
        compile=False,
        custom_objects={
            "GraphConvolution": GraphConvolution,
            "Antoine": Antoine,
            "GlorotUniform": initializers.GlorotUniform,
        },
    )
    model2 = load_model(
        opj(DIR_PATH, "models", "32.h5"),
        compile=False,
        custom_objects={
            "GraphConvolution": GraphConvolution,
            "Wagner": Wagner,
            "GlorotNormal": initializers.GlorotNormal,
        },
    )
    model3 = load_model(
        opj(DIR_PATH, "models", "33.h5"),
        compile=False,
        custom_objects={
            "GraphConvolution": GraphConvolution,
            "KingAlNajjar": KingAlNajjar,
            "GlorotNormal": initializers.GlorotNormal,
        },
    )

    return model1, model2, model3


def cal_vapor_pressure(SMILES, T, models=None):
    try:
        T = float(T)
    except:
        raise ValueError()

    if models is None:
        (
            model1,
            model2,
            model3,
        ) = load_models()
    else:
        model1, model2, model3 = models

    h1, a1 = get_input(SMILES, 25)

    h = h1.reshape(1, 25, 33)  # shape=(batch_number, 25 ,33)
    a = a1.reshape(1, 25, 25)  # shape=(batch_number, 25, 25)
    t = np.array([T])

    # # Predict the vapor pressure
    ln_vp1 = model1.predict([h, a, t], verbose=0)[0][0]  # shape=(batch_number, 1)
    ln_vp2 = model2.predict([h, a, t], verbose=0)[0][0]
    ln_vp3 = model3.predict([h, a, t], verbose=0)[0][0]

    ln_vp = 0.249305638 * ln_vp1 + 0.234145692 * ln_vp2 + 0.51654867 * ln_vp3
    vp = np.exp(ln_vp * np.log(10)) / 1000  # [kPa]

    return vp
