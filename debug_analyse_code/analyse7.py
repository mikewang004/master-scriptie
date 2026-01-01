import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import re 
import pandas as pd 
from tqdm import tqdm
import ctypes
import functools
import sys

def gaussian(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0)**2 / (2 * sigma**2))

# Load the shared library
def c_lib_init():
    lib = ctypes.CDLL('./boxAlgorithmsInC.so')

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



def filter_out_subset(data, combination):
    subset = data[(data['xid'] == combination[0]) & (data['yid'] == combination[1]) & (data['zid'] == combination[2])]
    return subset

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
    result_ptr = lib.find_nearest_value(nearest_values_ptr, a_size, data_ptr, n_size)
    
    # Convert result pointer to numpy array
    # This creates a view without copying data
    result_array = np.ctypeslib.as_array(result_ptr, shape=(n_size,))
    
    # If you need to take ownership and free later, make a copy:
    result_copy = result_array.copy()
    
    # Free the C-allocated memory
    #lib.free_results(result_ptr)
    
    return result_copy

def calc_nematic_tensor(array):
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

def get_time_temp_from_slurm(file_to_path):
    n_columns = 2
    with open(file_to_path, "r") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if "Per MPI rank memory allocation (min/avg/max)" in line:
            start_index = i+2
            break
    collected_rows = []
    for line in lines[start_index:]:
        if "error: *** JOB" in line:
            print(line)
            break
        row = line.strip().split()
        collected_rows.append(row[:n_columns])
    return np.array(collected_rows, dtype = float)

class atom_coords:

    """Read in files and do basic data preprocessing"""

    def __init__(self, file_to_path, nridges = 33):
        #Read in data 
        self.nridges = nridges
        self.datapd = self.prepare_position_data(file_to_path)
        self.volume, self.maxlength, self.minlength, self.total_volume_length = self.get_volume_box(file_to_path)
        self.current_timestep = self.get_timestep_from_file_name(file_to_path) 
        #Calculate box properties 
        self.n_atoms = len(self.datapd.index)
        self.no_polymers = self.datapd["mol_id"].max()
        self.polymer_length = self.n_atoms/self.no_polymers
        self.combinations = self.generate_box_list()
        self.bond_vectors = self.calculate_bond_vectors()
        self.datapd, self.local_volume = self.assign_center_of_mass(nridges = nridges)
        self.results = results()


    def prepare_position_data(self, file_to_path):
        """Reads position data (starts from line 9) in given a file name"""
        datapd = pd.read_csv(file_to_path, sep = " ", header = None, skiprows = 9)
        datapd.columns = ["atom_id", "mol_id", "xu", "yu", "zu"]
        datapd = datapd.sort_values("atom_id")
        datapd = datapd.set_index("atom_id")

        return datapd

    def get_volume_box(self, file_to_path):
        """Calculates volume of the total (!) simulation box."""
        datapd_first_rows = (pd.read_csv(file_to_path, sep = " ", header = None, skiprows = 5, nrows = 3))
        minlength = float(datapd_first_rows.iloc[0, 0])
        maxlength = float(datapd_first_rows.iloc[0, 1])
        total_volume_length = maxlength - minlength
        datapd_first_rows = np.abs(datapd_first_rows)
        axes_length = datapd_first_rows.sum(axis = 1)
        volume = (axes_length[0] * axes_length[1] * axes_length[2])
        return volume, maxlength, minlength, total_volume_length

    def get_timestep_from_file_name(self, file_to_path):
        """Only reads timestep from the file name."""
        pattern = r"_(\d+)\.txt$"
        timestep = re.search(pattern, file_to_path).group(1)
        return int(timestep)


    def calculate_bond_vectors(self):
        """Returns np.diff in positions except for each 100th particle."""
        # Filter out every nth row
        data = self.datapd 
        bond_vectors_array = np.diff(data.iloc[:, 1:4], axis = 0)
        polymer_mask = np.ones(len(bond_vectors_array), dtype = bool)
        polymer_mask[99::100] = False
        bond_vectors_array = bond_vectors_array[polymer_mask]
        bond_vectors = data.iloc[:, :5].copy()
        bond_vectors = bond_vectors[bond_vectors.index % 100 != 0]
        bond_vectors.iloc[:, 1:4] = bond_vectors_array
        bond_vectors = bond_vectors.rename(columns = {"xu" : "bondx", "yu" : "bondy", "zu" : "bondz"})
        return bond_vectors

    def wrap_coordinates(self, data):
        """Converts coordinates to wrapped coordinates. Input can be float or array (float)"""
        return (data - self.minlength) % self.total_volume_length + self.minlength

    def assign_center_of_mass(self, nridges = 33):
        """Loops over all polymers to assign center of mass 
        Also assignes a box id to each polymer group
        Takes about 110 seconds over dataset 720k big"""

        data = self.datapd
        length_array = np.linspace(self.minlength, self.maxlength, nridges+1)
        midpoint_ridges = ((length_array + np.roll(length_array, 1))/2)[1:] #Serves as box id 
        box_length =  (self.total_volume_length)/nridges
        df_com = pd.DataFrame(np.zeros([data.shape[0], 3]), columns = ["xid", "yid", "zid"], index = data.index)
        #print(self.total_volume_length, self.minlength, self.maxlength)
        local_volume = box_length**3
        wrapped_coordinates = self.wrap_coordinates(data.iloc[:, 1:5])
        print("test here below")
        print(wrapped_coordinates)
        df_com.iloc[:, 0] = find_box_id(midpoint_ridges, wrapped_coordinates.iloc[:, 0].to_numpy())
        df_com.iloc[:, 1] = find_box_id(midpoint_ridges, wrapped_coordinates.iloc[:, 1].to_numpy())
        df_com.iloc[:, 2] = find_box_id(midpoint_ridges, wrapped_coordinates.iloc[:, 2].to_numpy())

            
        data = pd.concat([data, df_com], axis=1)
        #data = data.iloc[:-1, :]
        #data.to_csv("data_test_ctypes.txt", sep = " ", mode = "w")
        return data, local_volume


    def generate_box_list(self, nridges = 33):
        numbers = np.arange(0, nridges)  # 0 to 32 inclusive
        combinations = np.array(np.meshgrid(numbers, numbers, numbers)).T.reshape(-1, 3)
        return combinations

    def create_new_polymer_df(self, column_names):
        """Creates an empty dataframe per mol_id with len(column_names) columns. column_names must be a list"""
        new_df = pd.DataFrame(np.zeros([self.no_polymers, len(column_names)]), columns = column_names)
        new_df.index.name = "mol_id"
        new_df.index = new_df.index + 1
        return new_df


    def get_nematic_vector_4(self, nridges = 33, save_ev = False, save_string = None):
        data = self.datapd
        # Prepare masks of all possible combinations 
        data = data[data.index % 100 != 0] # Filter out all last monomers as they do not have a bond vector per definiton
        df_cryst = pd.DataFrame(np.zeros([self.combinations.shape[0], 7]), columns = ["xid", "yid", "zid", "cryst_bool", "x_ev", "y_ev", "z_ev"])
        df_labdas = pd.DataFrame(np.zeros([self.combinations.shape[0], 3]), columns = ["labda_1", "labda_2", "labda_3"])
        df_cryst.iloc[:, :3] = self.combinations
        #for t in tqdm(range(0, len(self.combinations))):
        for t in range(0, len(self.combinations)):
            combination = self.combinations[t]
            #print(combination)
            #subset = data[(data['xid'] == combination[0]) & (data['yid'] == combination[1]) & (data['zid'] == combination[2])]
            subset = filter_out_subset(data, combination)
            if subset.empty == False:
                # Get index molecules 
                indexes = subset.index
                #print(indexes)
                subset_bond_vectors = self.bond_vectors.loc[indexes]
                order_param, order_ev, labda, ev = calc_nematic_tensor_2(subset_bond_vectors.iloc[:, 1:4])
                df_cryst.iloc[t,3] = order_param
                df_cryst.iloc[t,4:7] = order_ev
                df_labdas.iloc[t, :] = labda
        self.df_labdas = df_labdas
        self.fraction_crystallinity = fraction_crystallinity(df_cryst.iloc[:,3])
        self.df_cryst = df_cryst
        print(self.df_labdas)
        #print(self.fraction_crystallinity)
        if save_ev == True:
            df_labdas.to_csv("10e5_T1_debug_labdas.txt", sep = " ", mode = "w")
            df_cryst.to_csv("%s" %save_string, sep = " ", mode = "w")

        return self.fraction_crystallinity


    def merge_boxes(self, ndot_cutoff = 0.97, nridges = 33, cryst_cutoff = 0.8):
        max_labels = 837
        cryst_array = self.df_cryst.iloc[:, 1:].to_numpy()
        rows, cols = cryst_array.shape[0], cryst_array.shape[1]
        cryst_array_c = cryst_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        label_matrix_np = np.zeros([nridges, nridges, nridges], dtype = int)
        label_matrix = (ctypes.POINTER(ctypes.POINTER(ctypes.c_int)) * nridges)()
        new_label_matrix = (ctypes.POINTER(ctypes.POINTER(ctypes.c_int)) * nridges)()
        labels = np.zeros(max_labels, dtype = np.int32)
        for i in range(nridges):
            label_matrix[i] = (ctypes.POINTER(ctypes.c_int) * nridges)()
            new_label_matrix[i] = (ctypes.POINTER(ctypes.c_int) * nridges)()
            for j in range(nridges):
                label_matrix[i][j] = (ctypes.c_int * nridges)()
                new_label_matrix[i][j] = (ctypes.c_int * nridges)()
                for k in range(nridges):
                    label_matrix[i][j][k] = label_matrix_np[i, j, k]
                    new_label_matrix[i][j][k] = label_matrix_np[i, j, k]
        lib.hoshen_kopelman_crystallisation(cryst_array_c, rows, cols, nridges, cryst_cutoff, ndot_cutoff, label_matrix, new_label_matrix, labels)

        #Return label matrix to c 

        for i in range(nridges):
            for j in range(nridges):
                for k in range(nridges):
                    label_matrix_np[i,j,k] = label_matrix[j][i][k]



        label_matrix = label_matrix_np
        size = nridges
        labels2 = labels[labels != 0]


        # Loop over matrix to set all points not connected to a cluster to 0 

        for i in range(0, nridges):
            for j in range(0, nridges):
                for k in range(0, nridges):
                    if label_matrix[i,j,k] >0:
                        if label_matrix[(i -1)% nridges, j,k] ==0 and label_matrix[(i +1)% nridges, j,k] == 0:
                            if label_matrix[i, (j - 1) % nridges, k]==0 and label_matrix[i, (j + 1) %nridges, k] == 0:
                                if label_matrix[i,j,(k -1) %nridges]==0 and label_matrix[i,j, (k + 1) %nridges] == 0:
                                    label_matrix[i,j,k] = 0
                        if label_matrix[i,j,k] > 0:
                            print("position %i %i %i, current label %i, neighbours in x-direction %i and %i" %(i,j,k,label_matrix[i,j,k], label_matrix[(i -1)% nridges, j,k], label_matrix[(i +1)% nridges, j,k]))
                            print("neighbours in y-direction %i and %i" %(label_matrix[i, (j - 1) %nridges, k], label_matrix[i, (j + 1) %nridges, k]))
                            print("neighbours in z-direction %i and %i" %(label_matrix[i, j, (k - 1)%nridges], label_matrix[i, j, (k + 1)%nridges]))


        unique_values, counts = np.unique(label_matrix, return_counts=True)
        for value, count in zip(unique_values, counts):
            print(f"  {value}: {count}")
        print(np.sum(counts[1:]))





class polymer_properties():
    def __init__(self, atom_coords):
        self.atom_coords = atom_coords
        self.results = results()

    def create_new_polymer_df(self, column_names):
        """Creates an empty dataframe per mol_id with len(column_names) columns. column_names must be a list"""
        new_df = pd.DataFrame(np.zeros([self.atom_coords.no_polymers, len(column_names)]), columns = column_names)
        new_df.index.name = "mol_id"
        new_df.index = new_df.index + 1
        return new_df

    def end_to_end_distance(self, show_plot = False, save_plot_string = None):
        # Calculates end-to-end distance of each polymer 
        #print(self.datapd)
        df_end_end_length = self.create_new_polymer_df(["end_end_length"])
        for i in range(0, self.atom_coords.no_polymers):
            # First calculate end to end distance 
            # defined as r_n - r_i for each position 
            subset = self.bond_vectors[(self.bond_vectors["mol_id"] == i+1)]
            dist = np.sum(subset.iloc[:, 1:4])
            df_end_end_length.iloc[i] = (dist.iloc[0]* dist.iloc[0] + dist.iloc[1] * dist.iloc[1] + dist.iloc[2] * dist.iloc[2])
        end_to_end_length = (np.sum(df_end_end_length.to_numpy())/self.no_polymers)
        end_end_distance_normalised = np.sqrt(df_end_end_length.iloc[:])
        self.results.end_to_end_length = end_to_end_length
        self.mean_squared_end_to_end = np.sqrt(end_to_end_length)
        print("mean end-to-end is %f" %self.mean_squared_end_to_end)

    def gyration_radius(self, nridges = 33, show_plot = False):
        """Should confirm to Saras 2018 paper eq. 3"""
        df_gyration_radius = self.create_new_polymer_df(["comx", "comy", "comz", "gyration_radius"])
        counter = 0
        for i in range(0, self.atom_coords.no_polymers):
            subset = self.atom_coords.datapd[(self.atom_coords.datapd["mol_id"] == i+1)].iloc[:, 1:4]

            com = np.mean(subset, axis = 0)

            df_gyration_radius.iloc[i, :3] = com
            # Shift system to have new center of mass as center 
            subset_com = subset - com
            gyration_radius_squared = np.sum((subset_com**2), axis = 1)
            #print(np.mean(gyration_radius_squared))
            df_gyration_radius.iloc[i, 3] = np.mean(gyration_radius_squared) 
        self.results.gyration_radius = df_gyration_radius
        self.mean_gyration_radius = np.mean(df_gyration_radius["gyration_radius"]) #Ensemble average
        print("mean gyration is %f" %np.sqrt(self.mean_gyration_radius))


    def bond_bond_correlation(self, show_plot = False):
        #df_mean_bond_bond_per_polymer = self.create_new_polymer_df(["cos_theta"])
        df_bond_per_position = pd.DataFrame(np.zeros([int(self.atom_coords.n_atoms/self.atom_coords.no_polymers)-2, self.atom_coords.no_polymers]))
        df_bond_per_position.index.name = "bead_position"
        df_bond_per_position.index = df_bond_per_position.index + 1
        for i in range(0, self.atom_coords.no_polymers):
            subset = self.atom_coords.bond_vectors[(self.atom_coords.bond_vectors["mol_id"] == i + 1)].to_numpy()[:, 1:4]
            #print(subset)
            #Normalise bond vectors 
            subset = subset/np.linalg.norm(subset, axis = 1, keepdims = True)
            #print(subset)
            bond_bond_array = lib.inner_products_per_polymer(subset, subset.shape[0], subset.shape[1])
            bond_bond_array = np.ctypeslib.as_array(bond_bond_array, shape=(subset.shape[0]-1,))
            #print(bond_bond_array)
            df_bond_per_position.iloc[:, i] = bond_bond_array
            #print(df_bond_per_position.iloc[:, i])
        cos_per_position = np.mean(df_bond_per_position, axis = 1)
        print(cos_per_position)




class results(object):
    def __init__(self):
        return

    # Include plotting functions here 


lib = c_lib_init()