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
    max_labda = np.max(np.abs(labda)) #Try also without the np.abs()
    max_ev = ev[:, np.argmax(labda)]
    #max_ev = np.abs(ev[:, np.argmax(labda)])
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
    bond_vectors_indexed = bond_vectors.reset_index().rename({"index": "atom_id"}).set_index(["xid", "yid", "zid"]).sort_index()
    g = bond_vectors_indexed.groupby(level=["xid", "yid", "zid"])
    df_cryst = g.apply(calc_nematic_tensor_pandas, bond_vectors)
    df_cryst = df_cryst.reset_index()
    df_cryst = df_cryst.sort_values(by=["zid", "xid", "yid"]).reset_index(drop=True)
    return df_cryst