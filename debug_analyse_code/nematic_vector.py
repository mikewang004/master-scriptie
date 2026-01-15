import numpy as np 
import scipy as sp 
from numba import jit

@jit
def calc_nematic_tensor_3(array):
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
    return max_labda, max_ev


def nematic_vector_loop(cryst_array, data, bond_vectors):
    for t in range(0, cryst_array.shape[0]):
    #for t in range(0, 10):
        # Filter out subset from data
        mask = (data[:, 5] == cryst_array[t,0]) & (data[:, 6] == cryst_array[t,1]) & (data[:, 7] == cryst_array[t,2])
        subset = data[mask]
        if subset.size > 0:
            subset_bond_vectors = bond_vectors.loc[subset[:, 0]]
            order_param, order_ev = calc_nematic_tensor_3(subset_bond_vectors.iloc[:, 1:4].to_numpy())
            cryst_array[t, 3] = order_param
            cryst_array[t, 4:7] = order_ev
    return cryst_array