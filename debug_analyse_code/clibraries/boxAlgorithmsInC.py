import numpy as np
import ctypes
import os




def box_algos_lib(rel_path):
    #lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'./libboxAlgorithmsInC.so'))
    lib = ctypes.CDLL(os.path.join(rel_path, "libboxAlgorithmsInC.so"))


    # Define the function signature
    lib.find_nearest_value.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # const double nearest_values[]
        ctypes.c_size_t,                   # size_t a_size
        ctypes.POINTER(ctypes.c_double),  # const double data[]
        ctypes.c_size_t                    # size_t n_size
    ]


    lib.inner_products_per_polymer.argtypes = [
        np.ctypeslib.ndpointer(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_int
    ]

    lib.inner_products_columnwise_array.argtypes = [
        np.ctypeslib.ndpointer(ctypes.c_double), #Input array 1 of size [rows x cols]
        np.ctypeslib.ndpointer(ctypes.c_double), #Input array 2 of size [rows x cols]
        ctypes.c_int,
        ctypes.c_int
    ]
    lib.bond_bond_correlation.argtypes = [
        np.ctypeslib.ndpointer(ctypes.c_double), 
        np.ctypeslib.ndpointer(ctypes.c_double), 
        ctypes.c_int,
        ctypes.c_int
    ]
    lib.find_nearest_value.restype = ctypes.POINTER(ctypes.c_int)  # int* return type
    lib.inner_products_per_polymer.restype = ctypes.POINTER(ctypes.c_double)
    lib.inner_products_columnwise_array.restype = ctypes.POINTER(ctypes.c_double)
    lib.bond_bond_correlation.restype = None
    #lib.hoshen_kopelman_crystallisation.restype = ctypes.c_double
    return lib

def hoshen_kopelman_lib(rel_path):
    lib = ctypes.CDLL(os.path.join(rel_path, "libhoshenKopelmanInC.so"))
    lib.hoshen_kopelman_crystallisation.argtypes = [
        ctypes.POINTER(ctypes.c_double), #Input cryst_array [box-id, cryst_value, xev, yev, zev]
        ctypes.c_int, #no. rows
        ctypes.c_int, #no. cols
        ctypes.c_int, #nridges
        ctypes.c_float, #cutoff for crystallisation yes/no
        ctypes.c_float, #cutoff for ndot product 
        ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_int))), #Return array of size nridges**#
        ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_int))), #Return array of size nridges**#
        #(ctypes.POINTER(ctypes.c_int)) #Return array containing cluster indexes and size
        np.ctypeslib.ndpointer(ctypes.c_int, flags = "C_CONTIGUOUS")
    ]
    lib.hoshen_kopelman_crystallisation.restype = None
    return lib


def nematic_vector_lib(rel_path):
    lib = ctypes.CDLL(os.path.join(rel_path, "libnematicVector.c"))


here = os.path.dirname(os.path.abspath(__file__))  # directory of this .py file
box_algos_lib = box_algos_lib(here)
hk_lib = hoshen_kopelman_lib(here)