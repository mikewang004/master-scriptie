import numpy as np
import ctypes
#import os

def box_algos_lib_init():
    #lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'./libboxAlgorithmsInC.so'))
    lib = ctypes.CDLL('./libboxAlgorithmsInC.so')


    # Define the function signature
    lib.find_nearest_value.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # const double nearest_values[]
        ctypes.c_size_t,                   # size_t a_size
        ctypes.POINTER(ctypes.c_double),  # const double data[]
        ctypes.c_size_t                    # size_t n_size
    ]
    
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
    lib.find_nearest_value.restype = ctypes.POINTER(ctypes.c_int)  # int* return type
    lib.hoshen_kopelman_crystallisation.restype = None
    lib.inner_products_per_polymer.restype = ctypes.POINTER(ctypes.c_double)
    lib.inner_products_columnwise_array.restype = ctypes.POINTER(ctypes.c_double)
    #lib.hoshen_kopelman_crystallisation.restype = ctypes.c_double
    return lib



box_algos_lib = box_algos_lib_init()