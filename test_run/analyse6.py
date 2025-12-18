import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import re 
import pandas as pd 
from tqdm import tqdm
#from find_box_id import *
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


        
def find_nearest_array(nearest_values, data):
    """nearest value: 1d np array of size [a], data: 1d np array of size [n]"""
    results = np.zeros(data.size)
    for i in range(0, data.size):
        idx = np.abs(nearest_value - data).argmin()
        results[i] = idx
    return results

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

    """Used to read in files and analyse them"""

    def __init__(self, file_to_path, nridges = 33):
        #Note datapd has following columns: 
        self.datapd = self.prepare_position_data(file_to_path)
        self.n_atoms = len(self.datapd.index)
        self.no_polymers = self.datapd["mol_id"].max()
        self.polymer_length = self.n_atoms/self.no_polymers
        combinations = self.generate_box_list()
        self.bond_vectors = self.calculate_bond_vectors()
        self.get_volume_box(file_to_path)
        self.get_timestep_from_file_name(file_to_path)
        self.datapd = self.assign_center_of_mass(nridges = nridges)
        #self.box_atom_list = divide_into_box(self.datapd)
        #self.check_box_atom_list_exist()

    def read_in_temperature_scale(self, file_to_nano_slurm):
        pass
        # Should read temperature per timestep from the Slurm file 
        # Returns a 2d array with the temperature per timestep 

    def read_cryst(self, location):
        self.df_cryst = pd.read_csv(location, sep = " ", header = None, skiprows = 1)
        return 0;
        
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
        self.minlength = float(datapd_first_rows.iloc[0, 0])
        self.maxlength = float(datapd_first_rows.iloc[0, 1])
        datapd_first_rows = np.abs(datapd_first_rows)
        axes_length = datapd_first_rows.sum(axis = 1)
        self.volume = (axes_length[0] * axes_length[1] * axes_length[2])
        return self.volume

    def get_timestep_from_file_name(self, file_to_path):
        """Only reads timestep from the file name."""
        pattern = r"_(\d+)\.txt$"
        timestep = re.search(pattern, file_to_path).group(1)
        self.current_timestep = int(timestep)
        return timestep

    def check_box_atom_list_exist(self):
        """Subdivides box into [n] smaller boxes of a volume V"""
        try: 
            self.box_atom_list
        except AttributeError:
            self.box_atom_list, self.volume_box = self.divide_into_box()
        return 0;

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
        bond_vectors = bond_vectors.rename(columns = {"xu" : "x", "yu" : "y", "zu" : "z"})
        #np.savetxt("bond_vectors.txt", bond_vectors)
        return bond_vectors

    def wrap_coordinates(self, data):
        """Converts coordinates to wrapped coordinates. Input can be float or array (float)"""
        return (data - self.minlength) % self.total_volume_length + self.minlength

    def assign_center_of_mass(self, nridges = 33):
        """Loops over all polymers to assign center of mass 
        Also assignes a box id to each polymer group
        Takes about 110 seconds over dataset 720k big"""
        self.total_volume_length = self.maxlength - self.minlength
        data = self.datapd
        length_array = np.linspace(self.minlength, self.maxlength, nridges+1)
        self.midpoint_ridges = ((length_array + np.roll(length_array, 1))/2)[1:] #Serves as box id 
        box_length =  (self.total_volume_length)/nridges
        df_com = pd.DataFrame(np.zeros([data.shape[0], 3]), columns = ["xid", "yid", "zid"], index = data.index)
        print(self.total_volume_length, self.minlength, self.maxlength)
        self.local_volume = box_length**3
        # Wrap coordinates 
        #data.iloc[:, 2:5] = (data.iloc[:, 2:5] - self.minlength) % self.total_volume_length + self.minlength

        #Wrapped implementation below 
        ###############################################
        #data.iloc[:, 1:5] = self.wrap_coordinates(data.iloc[:, 1:5])
        # data.rename(columns={"x": "xu", "y": "yu", "z": "zu"})
        # df_com.iloc[:, 0] = find_box_id(self.midpoint_ridges, data.iloc[:, 1].to_numpy())
        # df_com.iloc[:, 1] = find_box_id(self.midpoint_ridges, data.iloc[:, 2].to_numpy())
        # df_com.iloc[:, 2] = find_box_id(self.midpoint_ridges, data.iloc[:, 3].to_numpy())
        # for i in range(0, data.shape[0]): #Below is an all-python approach 
        #     df_com.iloc[i, 0] = find_nearest(self.midpoint_ridges, data.iloc[i, 1])
        #     df_com.iloc[i, 1] = find_nearest(self.midpoint_ridges, data.iloc[i, 2])
        #     df_com.iloc[i, 2] = find_nearest(self.midpoint_ridges, data.iloc[i, 3])
        #########################################################
        #Implementation below preserves unwrapped original coordinates
        wrapped_coordinates = self.wrap_coordinates(data.iloc[:, 1:5])
        print("test here below")
        print(wrapped_coordinates)
        df_com.iloc[:, 0] = find_box_id(self.midpoint_ridges, wrapped_coordinates.iloc[:, 0].to_numpy())
        df_com.iloc[:, 1] = find_box_id(self.midpoint_ridges, wrapped_coordinates.iloc[:, 1].to_numpy())
        df_com.iloc[:, 2] = find_box_id(self.midpoint_ridges, wrapped_coordinates.iloc[:, 2].to_numpy())

            
        data = pd.concat([data, df_com], axis=1)
        #data = data.iloc[:-1, :]
        #data.to_csv("data_test_ctypes.txt", sep = " ", mode = "w")
        return data


    def generate_box_list(self, nridges = 33):
        numbers = np.arange(0, nridges)  # 0 to 32 inclusive
        self.combinations = np.array(np.meshgrid(numbers, numbers, numbers)).T.reshape(-1, 3)
        return self.combinations

    def create_new_polymer_df(self, column_names):
        """Creates an empty dataframe per mol_id with len(column_names) columns. column_names must be a list"""
        new_df = pd.DataFrame(np.zeros([self.no_polymers, len(column_names)]), columns = column_names)
        new_df.index.name = "mol_id"
        new_df.index = new_df.index + 1
        return new_df

    def end_to_end_distance(self, nridges = 33, show_plot = False, save_plot_string = None):
        # Calculates end-to-end distance of each polymer 
        #print(self.datapd)
        df_end_end_length = self.create_new_polymer_df(["end_end_length"])
        for i in range(0, self.no_polymers):
            # First calculate end to end distance 
            # defined as r_n - r_i for each position 
            subset = self.bond_vectors[(self.bond_vectors["mol_id"] == i+1)]
            first_element = subset.iloc[0, 1:4] #r_{i,1}
            last_element = subset.iloc[-1, 1:4] #r_{i,N}
            dist = np.sum(subset.iloc[:, 1:4])
            df_end_end_length.iloc[i] = (dist.iloc[0]* dist.iloc[0] + dist.iloc[1] * dist.iloc[1] + dist.iloc[2] * dist.iloc[2])
        end_to_end_length = (np.sum(df_end_end_length.to_numpy())/self.no_polymers)
        end_end_distance_normalised = np.sqrt(df_end_end_length.iloc[:])
        self.mean_squared_end_to_end = np.sqrt(end_to_end_length)
        print("mean end-to-end is %f" %self.mean_squared_end_to_end)

        if show_plot == True:
            # Fit gaussian 
            values, bins, _ = plt.hist(end_end_distance_normalised, bins = 200, density = True)
            mu, sigma = sp.stats.norm.fit(end_end_distance_normalised)
            gauss_fit = sp.stats.norm.pdf(bins, mu, sigma)
            plt.plot(bins, gauss_fit, label = r"Gaussian fit, $\mu = %.2f$, $\sigma = %.2f$" %(mu, sigma))
            plt.vlines((self.mean_squared_end_to_end), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = r"mean end-to-end length = %.2f" %self.mean_squared_end_to_end)
            plt.title("Distribution of end-to-end distance, PVA-100, simulation start")
            plt.xlabel(r"end-end distance ($R_e\sigma$)")
            #plt.ylabel("count")
            plt.legend()
            plt.savefig("end_to_end_distance_PVA-100_start.pdf")
            plt.show()
        return self.mean_squared_end_to_end;

    
    def distribution_bond_vectors(self, show_plot = False):
        for i in range(0, self.no_polymers):
            subset = self.bond_vectors[(self.bond_vectors["mol_id"] == i + 1)]
            

    # def end_to_end_distance_2(self, nridges = 33, show_plot = False):
    #     df_end_end_length = self.create_new_polymer_df(["end_end_length"])
    #     for i in range(0, 3):
    #     #for i in range(0, self.no_polymers):
    #         subset = self.datapd[(self.datapd["mol_id"] == i+1)].iloc[:, 1:4]
    #         first_element = subset.iloc[0, 1:4]
    #         dist = subset.iloc[1:, 1:4] - first_element
    #         dist_sum = np.sum(dist**2, axis = 1)
    #         mean_squared_end_to_end = np.sum(dist_sum)
    #         df_end_end_length.iloc[i] = mean_squared_end_to_end
    #     print(np.mean(df_end_end_length))
    

    def gyration_radius(self, nridges = 33, show_plot = False):
        """Should confirm to Saras 2018 paper eq. 3"""
        #TODO fix definition of rhs and gyration radius 
        # First calculate center of mass of each polymer 
        df_gyration_radius = self.create_new_polymer_df(["comx", "comy", "comz", "gyration_radius"])
        counter = 0
        for i in range(0, self.no_polymers):
            subset = self.datapd[(self.datapd["mol_id"] == i+1)].iloc[:, 1:4]

            com = np.mean(subset, axis = 0)

            df_gyration_radius.iloc[i, :3] = com
            # Shift system to have new center of mass as center 
            subset_com = subset - com
            gyration_radius_squared = np.sum((subset_com**2), axis = 1)
            #print(np.mean(gyration_radius_squared))
            df_gyration_radius.iloc[i, 3] = np.mean(gyration_radius_squared) 
        self.mean_gyration_radius = np.mean(df_gyration_radius["gyration_radius"]) #Ensemble average
        print("mean gyration is %f" %np.sqrt(self.mean_gyration_radius))
        if show_plot == True:
            values, bins, _ = plt.hist(np.sqrt(df_gyration_radius.iloc[:, -1]), bins = 200, density = True)
            #mu, sigma = sp.stats.norm.fit(np.sqrt(df_gyration_radius.iloc[:, -1]))
            #gauss_fit = sp.stats.norm.pdf(bins, mu, sigma)
            #plt.plot(bins, gauss_fit, label = r"Gaussian fit, $\mu = %.2f$, $\sigma = %.2f$" %(mu, sigma))
            plt.vlines(np.sqrt(self.mean_gyration_radius), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "mean gyration radius = %.4f" %np.sqrt(self.mean_gyration_radius))
            #plt.vlines(((np.mean(gyration_2[gyration_2 < 10.0]))), ymin = 0, ymax = np.max(values), linestyles =":", color = "red", label = "gaussian peak gyration radius = %.4f" %(np.mean(gyration_2[gyration_2 < 10.0])))
            plt.title("Distribution of the gyration radius, PVA-100, simulation start")
            plt.xlabel(r"gyration radius  ($R_g/\sigma$)")
            #plt.ylabel("probability")
            plt.legend()
            plt.savefig("gyration_radius_PVA-100_normalised_start.pdf")
            plt.show()

    def gyration_radius_debug(self, nridges = 33, show_plot = False):
        molid = 999
        for i in range(molid, molid +1):
            subset = self.datapd[(self.datapd["mol_id"] == i+1)].iloc[:, 1:4]
            print(subset)
            com = np.mean(subset, axis = 0)
            #print(com)
            #df_gyration_radius.iloc[i, :3] = com
            # Shift system to have new center of mass as center 
            subset_com = subset - com
            #subset_com = (subset - com) %self.total_volume_length
            print(subset_com)
            gyration_radius_squared = np.sum((subset_com**2), axis = 1)/self.polymer_length
            print(np.sqrt(gyration_radius_squared))
            

    def gyration_radius_rmsd(self, nridges = 33, show_plot = False):
        df_gyration_radius = self.create_new_polymer_df(["comx", "comy", "comz", "gyration_radius"])
        counter = 0
        for i in range(0, self.no_polymers):
        #for i in range(0, 50):
            subset = np.array(self.datapd[(self.datapd["mol_id"] == i+1)].iloc[:, 1:4])
            sqd = np.sum((np.abs(subset[:, np.newaxis] - subset))**2, axis = 2)
            upper_triangle = sqd[np.triu_indices(self.polymer_length, k=1)]

            #print(np.sum(upper_triangle)/(2 * self.polymer_length**2))
            gyration_radius_squared = np.sum(upper_triangle)/(2 * self.polymer_length**2)
            #print(i, np.mean(gyration_radius_squared) )
            df_gyration_radius.iloc[i, 3] = np.mean(gyration_radius_squared) 
        self.mean_gyration_radius = np.mean(df_gyration_radius["gyration_radius"]) #Ensemble average
        print("mean gyration is %f" %np.sqrt(self.mean_gyration_radius))
        gyration_2 = np.sqrt(df_gyration_radius.iloc[: , -1].to_numpy())

        if show_plot == True:
            values, bins, _ = plt.hist(np.sqrt(df_gyration_radius.iloc[:, -1]), bins = 50, density = True)
            plt.vlines(np.sqrt(self.mean_gyration_radius), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "mean gyration radius = %.4f" %np.sqrt(self.mean_gyration_radius))
            plt.vlines(((np.mean(gyration_2[gyration_2 < 10.0]))), ymin = 0, ymax = np.max(values), linestyles =":", color = "red", label = "gaussian peak gyration radius = %.4f" %(np.mean(gyration_2[gyration_2 < 10.0])))
            plt.title("Distribution of the gyration radius, PVA-100, simulation start")
            plt.xlabel("gyration radius")
            #plt.ylabel("probability")
            plt.legend()
            plt.savefig("gyration_radius_PVA-100_unnormalised_start.pdf")
            plt.show()



    def gyration_tensor(self, nridges = 33, show_plot = False):
        df_gyration_radius = self.create_new_polymer_df(["comx", "comy", "comz", "gyration_radius"])
        for i in range(0, self.no_polymers):
            subset = np.array(self.datapd[(self.datapd["mol_id"] == i+1)].iloc[:, 1:4])
            com = np.mean(subset, axis = 0)
            #print(com)
            df_gyration_radius.iloc[i, :3] = com
            # Shift system to have new center of mass as center 
            subset_com = subset - com
            # Normalise result 
            #subset_com = np.linalg.norm(subset_com, axis = 1, keepdims = True)
            gyration_tensor =  np.einsum('im,in->mn', subset_com,subset_com)/self.polymer_length
            #print(gyration_tensor)
            labda, ev = np.linalg.eig(gyration_tensor)
            #print(labda)
            gyration_radius_squared = np.sum(labda)
            df_gyration_radius.iloc[i, 3] = np.mean(gyration_radius_squared) 
        self.mean_gyration_radius = np.mean(df_gyration_radius["gyration_radius"]) #Ensemble average
        gyration_2 = np.sqrt(df_gyration_radius.iloc[: , -1].to_numpy())

        if show_plot == True:
            values, bins, _ = plt.hist(np.sqrt(df_gyration_radius.iloc[:, -1]), bins = 100, density = True)
            plt.vlines(np.sqrt(self.mean_gyration_radius), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "mean gyration radius = %.4f" %np.sqrt(self.mean_gyration_radius))
            plt.vlines(((np.mean(gyration_2[gyration_2 < 10.0]))), ymin = 0, ymax = np.max(values), linestyles =":", color = "red", label = "gaussian peak gyration radius = %.4f" %(np.mean(gyration_2[gyration_2 < 10.0])))
            plt.title("Distribution of the gyration radius, PVA-100, simulation start")
            plt.xlabel("gyration radius")
            #plt.ylabel("probability")
            plt.legend()
            plt.savefig("gyration_radius_PVA-100_unnormalised_start.pdf")
            plt.show()



    def bond_bond_correlation(self, show_plot = False):
        #df_mean_bond_bond_per_polymer = self.create_new_polymer_df(["cos_theta"])
        df_bond_per_position = pd.DataFrame(np.zeros([int(self.n_atoms/self.no_polymers)-2, self.no_polymers]))
        df_bond_per_position.index.name = "bead_position"
        df_bond_per_position.index = df_bond_per_position.index + 1
        for i in range(0, self.no_polymers):
            subset = self.bond_vectors[(self.bond_vectors["mol_id"] == i + 1)].to_numpy()[:, 1:4]
            print(subset)
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
        plt.scatter(cos_per_position.index, cos_per_position)
        plt.xlabel("n")
        plt.ylabel(r"cos\theta(n)")
        plt.title("Distribution of bond-bond correlations")
        plt.show()

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


    def get_distribution_eigenvalues(self, title, nridges = 33, savestring = None, readfile = None):
        #TODO modify this to include all eigenvalues instead of the highest one and normalise result
        #labdas = self.df_labdas.to_numpy().flatten()
        if readfile != None:
            self.df_labdas = pd.read_csv(readfile, sep = " ", header = None, skiprows = 1)
        labdas = self.df_labdas.to_numpy()[:, 1:]
        print(labdas)
        #Get largest labdas per row 
        max_indices = np.argmax(np.abs(labdas), axis = 1)

        max_labda = labdas[range(labdas.shape[0]), max_indices]
        values, bins, _ = plt.hist(max_labda, bins = 250, density = True)
        print(max_labda)
        plt.title(title)
        plt.xlabel("eigenvalues")
        plt.vlines(0.8, ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "cut-off crystalisation")
        plt.legend()
        if savestring == None:
            plt.show()
        else:
            plt.savefig("%s.pdf" %savestring)
        plt.close()
        #plt.show()

    def apply_nn_cutoff(self, ndot, ndot_cutoff = 0.97):
        if ndot >= cutoff:
            return 1
        else:
            return 0


    def calc_rdf(self):
        """Calculates radial distribution function"""

        data = self.datapd
        positions = data.iloc[:-1, 1:4].to_numpy()
        #dists = np.zeros([positions.shape[9**2])
        j = 0
        #for i in range(0, positions):
        current_position = positions[0, :]
        dists = np.abs(current_position - positions)
        print(dists)



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
        #print(labels[labels != 0])
        #print(labels)
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

                        # if label_matrix[(i -1)% nridges, j,k] or label_matrix[(i +1)% nridges, j,k] != 0:
                        #     print("position %i %i %i, current label %i, neighbours in x-direction %i and %i" %(i,j,k,label_matrix[i,j,k], label_matrix[(i -1)% nridges, j,k], label_matrix[(i +1)% nridges, j,k]))
                        # if label_matrix[i, (j - 1) % nridges, k] or label_matrix[i, (j + 1) %nridges, k] != 0:
                        #     print("position %i %i %i, current label %i, neighbours in y-direction %i and %i" %(i,j,k,label_matrix[i,j,k], label_matrix[i, (j - 1) %nridges, k], label_matrix[i, (j + 1) %nridges, k]))
                        # if label_matrix[i,j,(k -1) %nridges] or label_matrix[i,j, (k + 1) %nridges] != 0:
                        #     print("position %i %i %i, current label %i, neighbours in z-direction %i and %i" %(i,j,k,label_matrix[i,j,k], label_matrix[i, j, (k - 1)%nridges], label_matrix[i, j, (k + 1)%nridges]))


        unique_values, counts = np.unique(label_matrix, return_counts=True)
        #print(label_matrix[label_matrix != 0])
        for value, count in zip(unique_values, counts):
            print(f"  {value}: {count}")
        #print(np.unique(label_matrix))
        print(np.sum(counts[1:]))

        #print(label_matrix[:,:,0])
        #print(label_matrix[:,:,1])
        








            

    def get_density_dist(self, nridges = 33):
        """Uses new method with combinations to calculate local density (i.e. density per box)"""
        #print(self.combinations)
        local_densities = np.zeros([len(self.combinations), 4]) #First three columns reserved for combination, last for corresponding density
        local_densities[:, :3] = self.combinations
        data = self.assign_center_of_mass(nridges = nridges)
        for t in tqdm(range(0, len(self.combinations))):
            combination = self.combinations[t]
            subset = data[(data['xid'] == combination[0]) & (data['yid'] == combination[1]) & (data['zid'] == combination[2])]
            local_densities[t, 3] = len(subset.index)/self.local_volume
        return local_densities



def fraction_crystallinity(data, cutoff = 0.8):
    """Data 1d 1d list/array, defined as data > 0.8 -> crystallinity = 1;
    data <= 0.8 -> crystallinity = 0 as in Sommer/Luo Sep 2010"""

    mask = data > cutoff
    fraction = len(data[mask]) / len(data)
    return fraction

def plot_density_dist(atom_coord, title):
    """Returns histogram of local density per cube"""
    data = atom_coord.get_density_dist() #Returns [nx4] array with first 3 columns indicating box id while 4 contains local density
    plt.hist(data[:, 3], bins = 500)
    plt.title(title) #Include time and temperature in title 
    plt.xlabel("Local density / cube")
    plt.ylabel("count")
    #plt.xlim(2.0, 3.0)
    plt.savefig(title + ".pdf")
    plt.close()
    return 0;


def get_list_atom_coords(shared_name, n_samples, starttime = 0, endtime = 1e7):
    """Wrapper to extract the monomers distribution over the local boxes over multiple LAMPPS-files"""
    timespace = np.linspace(starttime, endtime, n_samples)
    #print(timespace)
    list_atom_coords = []
    for i in range(0, n_samples):
        list_atom_coords.append(atom_coords("%s_%i.txt" %(shared_name, timespace[i])))
    #print(list_atom_coords)
    return list_atom_coords


def get_list_different_tdot(shared_name, min_time, tdots):
    """tdots array, time int, shared_name string
    NB assumes that tdots is ordered from smallest total runtime to largest"""
    list_atom_coords = []
    for i in range(0, len(tdots)):
        current_tdot = "e-%i" %(tdots[i])
        time = min_time * 10**i
        list_atom_coords.append(atom_coords("%s_%s_time_%i.txt" %(shared_name, current_tdot, time)))
    return list_atom_coords

def get_crystallinity_tdots(shared_name, min_time, tdots):
    """tdots array, time int, shared_name string
    NB assumes that tdots is ordered from smallest total runtime to largest"""
    psic_array = np.zeros([2,len(tdots)])
    psic_array[0, :] = 10.0**(-tdots)
    list_atom_coords = get_list_different_tdot(shared_name, min_time, tdots)
    i = 0
    for atom_coords in list_atom_coords:
        psic_array[1, i] = atom_coords.get_nematic_vector_4()
        i = i + 1
    print(psic_array)
    return psic_array



def plot_volume_line(list_atom_coords, title,savestring = None, n_atoms = 720000, starttemp = 1.0, endtemp = 0.5):
    """Returns plot of volume as function of temperature"""
    list_volumes = []
    for atom_coords in list_atom_coords:
        list_volumes.append((atom_coords.volume)/n_atoms)
        print(atom_coords.volume/n_atoms)

    temps = np.linspace(starttemp, endtemp, num = len(list_atom_coords))

    plt.scatter(temps, list_volumes)
    plt.title(title)
    plt.xlabel("Temperature [unitless]")
    plt.ylabel("Volume per monomer")
    if savestring == None:
        pass
    else:
        plt.savefig("%s.pdf" %savestring)
        np.savetxt("%s.txt" %savestring)
    plt.show()

def plot_order_param(list_atom_coords, title,savestring = None, starttemp = 1.0, endtemp = 0.5, plot_time_instead = False, starttime = None, endtime = None, n_samples = None):
    """Returns plot of volume as function of temperature"""
    list_order_params = []
    for t in tqdm(range(0, len(list_atom_coords))):
        print(t)
        atom_coords = list_atom_coords[t]
    #for atom_coords in list_atom_coords:
        atom_coords.get_nematic_vector_4()
        list_order_params.append((atom_coords.fraction_crystallinity))
        print(atom_coords.fraction_crystallinity)

    list_order_params = np.array(list_order_params)
    if plot_time_instead == True:
        timespace = np.linspace(starttime, endtime, n_samples)
        plt.scatter(timespace, list_order_params)
        plt.xlabel("Simulation time")
    else:
        temps = np.linspace(starttemp, endtemp, num = len(list_atom_coords))
        plt.scatter(temps, list_order_params)
        plt.xlabel("Temperature [unitless]")
    plt.title(title)
    plt.ylabel("Fraction of crystallinity")
    if savestring == None:
        pass
    else:
        plt.savefig("%s.pdf" %savestring)
        #np.savetxt("%s.txt" %savestring, np.column_stack((temps, list_order_params)))
    plt.show()


def plot_multiple_dists_eigenvalues(list_atom_coords, starttemp = 1.0, endtemp = 0.5):
    temps = np.linspace(starttemp, endtemp, 21)
    for t in tqdm(range(0, len(list_atom_coords))):
        current_temp = temps[t]
        atom_coords = list_atom_coords[t]
        atom_coords.get_nematic_vector_4()
        title = r"Distribution of eigenvalues at $T = %s$, $\dot{T} = 10^{-5}$" %temps[t]
        savestring = "plots/eigenvalue_dist_pva_100_t_%s" %temps[t]
        atom_coords.get_distribution_eigenvalues(title = title, savestring = savestring)
    return 0;

lib = c_lib_init()





