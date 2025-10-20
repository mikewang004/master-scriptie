import ctypes
import numpy as np
import os
# Load the shared library
#lib = ctypes.CDLL('./find_box_id.so')
lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "find_box_id.so"))
# Define the function signature
lib.find_nearest_value.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # const double nearest_values[]
    ctypes.c_size_t,                   # size_t a_size
    ctypes.POINTER(ctypes.c_double),  # const double data[]
    ctypes.c_size_t                    # size_t n_size
]
lib.find_nearest_value.restype = ctypes.POINTER(ctypes.c_int)  # int* return type

# # Helper function to free memory (if you create one in C)
# lib.free_results.argtypes = [ctypes.POINTER(ctypes.c_int)]
# lib.free_results.restype = None


def find_box_id(nearest_values, data):
    """
    nearest_values: 1D numpy array of doubles
    data: 1D numpy array of doubles
    returns: 1D numpy array of integers
    """
    a_size = len(nearest_values)
    n_size = len(data)
    
    # Convert numpy arrays to ctypes pointers
    nearest_values_ptr = nearest_values.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    data_ptr = data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    # Call C function
    result_ptr = lib.find_nearest_array(nearest_values_ptr, a_size, data_ptr, n_size)
    
    # Convert result pointer to numpy array
    # This creates a view without copying data
    result_array = np.ctypeslib.as_array(result_ptr, shape=(n_size,))
    
    # If you need to take ownership and free later, make a copy:
    result_copy = result_array.copy()
    
    # Free the C-allocated memory
    # lib.free_results(result_ptr)
    
    return result_copy