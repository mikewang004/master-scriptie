import numpy as np 
import scipy as sp 
from numba import jit, njit
from tqdm import tqdm
import pandas as pd


def calc_nematic_tensor_2(array):
    """Calculation for the nematic tensor of a local box. NB this is not used yet in the analysis."""
    array_length = array.shape[0]
    array = array / np.linalg.norm(array, axis = 1, keepdims = True)
    Q = np.zeros([3,3])
    outer = (np.einsum('ni,nj->ij', array, array)) / array_length
    #Q =  np.mean(outer  - (1/3) * np.eye(3), axis = 0) # According to Sommer/Luo Sep 2010
    Q = 1.5 * outer - 0.5* np.eye(3) # Sara 2015
    #order_param = np.sqrt(1.5 * np.trace(Q**2)) #Sommer/Luo 2010
    labda, ev = np.linalg.eigh(Q)
    max_labda = np.max(labda)
    max_ev = ev[:, np.argmax(labda)]
    order_param = max_labda #Sara 2015
    return max_labda, max_ev, labda, ev


@njit
def norm_keep_dims_1d(A):
    n, m = A.shape
    result = np.empty((n, 1), dtype = A.dtype)
    for i in range(n):
        s = 0.0
        for j in range(m):
            s += A[i, j] * A[i, j]
        result[i, 0] = np.sqrt(s)
    return result




def calc_nematic_tensor_3(array):
    """Calculation for the nematic tensor of a local box. NB this is not used yet in the analysis."""
    array_length = array.shape[0]
    array = array / norm_keep_dims_1d(array)
    Q = np.zeros((3,3))
    outer = (array.T @ array) / array_length

    Q = 1.5 * outer - 0.5* np.eye(3) 
    labda, ev = np.linalg.eigh(Q)
    max_labda = np.max(labda)
    max_ev = ev[:, np.argmax(labda)]
    return max_labda, max_ev, labda, ev


def calc_nematic_tensor_pandas(block: pd.DataFrame, bond_vectors) -> pd.Series:
    array = bond_vectors.loc[block["atom_id"]].to_numpy()[:, 1:4]
    array_length = array.shape[0]
    if array_length == 0:
        return pd.Series(
            {"cryst_bool": 0.0, "x_ev": 0.0, "y_ev": 0.0, "z_ev": 0.0}
        )
    array = array / norm_keep_dims_1d(array)
    Q = np.zeros((3,3))
    outer = (array.T @ array) / array_length

    Q = 1.5 * outer - 0.5* np.eye(3) 
    labda, ev = np.linalg.eigh(Q)
    idx = np.argmax(labda)  # or np.argmax(np.abs(labda)), just be consistent
    max_labda = labda[idx]
    max_ev = ev[:, idx]

    if max_ev[2] < 0:
        max_ev = -max_ev
    return pd.Series({"cryst_bool": max_labda, "x_ev": max_ev[0], "y_ev":max_ev[1], "z_ev": max_ev[2]})



def nematic_vector_loop(data: pd.DataFrame, bond_vectors):
    data_indexed = data.reset_index().rename({"index": "atom_id"}).set_index(["xid", "yid", "zid"]).sort_index()
    g = data_indexed.groupby(level=["xid", "yid", "zid"])
    df_cryst = g.apply(calc_nematic_tensor_pandas, bond_vectors)
    df_cryst = df_cryst.reset_index()
    df_cryst = df_cryst.sort_values(by=["zid", "xid", "yid"]).reset_index(drop=True)
    #df_cryst.to_csv("fast_cryst.txt", sep = " ", mode = "w")
    return df_cryst




def nematic_vector_loop_2(bond_vectors):
    """Uses bond vector cell alignment to calculate crystallinity per cell"""
    bond_vectors_indexed = bond_vectors.reset_index().rename({"index": "atom_id"}).set_index(["xid", "yid", "zid"]).sort_index()
    g = bond_vectors_indexed.groupby(level=["xid", "yid", "zid"])
    df_cryst = g.apply(calc_nematic_tensor_pandas, bond_vectors)
    df_cryst = df_cryst.reset_index()
    df_cryst = df_cryst.sort_values(by=["zid", "xid", "yid"]).reset_index(drop=True)
    return df_cryst

def compute_Q(block: pd.DataFrame) -> pd.Series:
    vecs = block[['bx', 'by', 'bz']].to_numpy()
    if len(vecs) == 0:
        return pd.Series({f"Q{a}{b}": 0.0 for a in range(3) for b in range(3)})

    norms = np.linalg.norm(vecs, axis=1, keepdims=True)
    vecs_unit = vecs / norms
    N = vecs_unit.shape[0]
    outer = (vecs_unit.T @ vecs_unit) / N
    Q = 1.5 * outer - 0.5 * np.eye(3)

    # symmetric eigen-decomposition
    vals, vecs_ev = np.linalg.eigh(Q)               # vals ascending

    # largest eigenvalue & eigenvector
    idx = np.argmax(vals)
    #idx = np.argmax(np.abs(vals))
    S = vals[idx]                                   # order parameter
    director = vecs_ev[:, idx]                      # eigenvector

    # enforce a consistent sign (optional; e.g. nz >= 0)
    if director[2] < 0:
        director = -director

    return pd.Series({'cryst_bool': S, 'x_ev': director[0], 'y_ev': director[1], 'z_ev': director[2]})

def orderparameter(block: pd.DataFrame) -> pd.Series:
    """
    Python translation of your C code:
    - Build second moment matrix A = <u u^T>
    - Build nematic tensor Q = (3/2) A - (1/2) I
    - Diagonalize Q
    - Return nematic order parameter S and director

    Parameters
    ----------
    um : array_like, shape (npar, 3)
        Orientation vectors (ideally unit vectors).

    Returns
    -------
    A : (3, 3) ndarray
        Second moment matrix <u u^T>.
    Q : (3, 3) ndarray
        Nematic tensor.
    S : float
        Nematic order parameter (largest eigenvalue of Q in magnitude).
    director : (3,) ndarray
        Corresponding eigenvector (director).
    lambdas : (3,) ndarray
        Eigenvalues of Q (ascending order).
    evecs : (3, 3) ndarray
        Eigenvectors of Q as columns; evecs[:, k] corresponds to lambdas[k].
    ss : float
        sqrt(2/3 * sum(lambdas^2)), as in your C code.
    """

    um = block[['bx', 'by', 'bz']].to_numpy()
    if len(um) == 0:
        return pd.Series({f"Q{a}{b}": 0.0 for a in range(3) for b in range(3)})
    um = np.asarray(um, dtype=float)
    npar = um.shape[0]
    if npar == 0:
        raise ValueError("Need at least one vector")

    # --- 1) Build second moment matrix A = <u u^T> ---
    # This is the explicit analog of your a11, a12, ..., a33 sums.
    ux = um[:, 0]
    uy = um[:, 1]
    uz = um[:, 2]

    a11 = np.mean(ux * ux)
    a12 = np.mean(ux * uy)
    a13 = np.mean(ux * uz)
    a22 = np.mean(uy * uy)
    a23 = np.mean(uy * uz)
    a33 = np.mean(uz * uz)

    A = np.array([[a11, a12, a13],
                  [a12, a22, a23],
                  [a13, a23, a33]])

    # --- 2) Nematic tensor Q = (3/2) A - (1/2) I ---
    Q = 1.5 * A - 0.5 * np.eye(3)

    # --- 3) Eigenvalues/eigenvectors of Q (symmetric) ---
    lambdas, evecs = np.linalg.eigh(Q)  # eigenvalues in ascending order

    # --- 4) Pick nematic order parameter S and director ---
    idx = np.argmax(np.abs(lambdas))   # largest |lambda|
    S = np.abs(lambdas[idx])
    director = evecs[:, idx].copy()

    # Optional: fix sign convention for director (e.g. z >= 0)
    if director[2] < 0:
        director = -director

    # --- 5) ss = sqrt(2/3 * sum(lambda_i^2)) ---
    #ss = np.sqrt(2.0 / 3.0 * np.sum(lambdas**2))

    return pd.Series({'cryst_bool': S, 'x_ev': director[0], 'y_ev': director[1], 'z_ev': director[2]})