import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import re 
import pandas as pd 
from tqdm import tqdm
import functools
import sys
from clibraries.boxAlgorithmsInC import box_algos_lib, hk_lib
import ctypes
import time
from numba import jit
from nematic_vector import calc_nematic_tensor_2, nematic_vector_loop, nematic_vector_loop_2
from matplotlib.colors import ListedColormap, BoundaryNorm

def fraction_crystallinity(data, cutoff = 0.8):
    """Data 1d 1d list/array, defined as data > 0.8 -> crystallinity = 1;
    data <= 0.8 -> crystallinity = 0 as in Sommer/Luo Sep 2010"""

    mask = data > cutoff
    fraction = len(data[mask]) / len(data)
    return fraction, len(data[mask])

def gaussian(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0)**2 / (2 * sigma**2))





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
    result_ptr = box_algos_lib.find_nearest_value(nearest_values_ptr, a_size, data_ptr, n_size)
    
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

    def __init__(self, path_to_file, nridges = 33):
        #Read in data 
        self.nridges = nridges
        self.path_to_file = path_to_file
        self.datapd = self.prepare_position_data(path_to_file)
        self.volume, self.maxlength, self.minlength, self.total_volume_length = self.get_volume_box(path_to_file)
        self.current_timestep = self.get_timestep_from_file_name(path_to_file) 
        #Calculate box properties 
        self.n_atoms = len(self.datapd.index)
        self.no_polymers = self.datapd["mol_id"].max()
        self.polymer_length = int(self.n_atoms/self.no_polymers)
        self.combinations = self.generate_box_list()
        self.bond_vectors = self.calculate_bond_vectors()
        self.datapd, self.local_volume, _ = self.assign_center_of_mass(nridges = nridges)
        self.results = results()


    def prepare_position_data(self, path_to_file):
        """Reads position data (starts from line 9) in given a file name"""
        datapd = pd.read_csv(path_to_file, sep = " ", header = None, skiprows = 9)
        datapd.columns = ["atom_id", "mol_id", "xu", "yu", "zu"]
        datapd = datapd.sort_values("atom_id")
        datapd = datapd.set_index("atom_id")

        return datapd

    def get_volume_box(self, path_to_file):
        """Calculates volume of the total (!) simulation box."""
        datapd_first_rows = (pd.read_csv(path_to_file, sep = " ", header = None, skiprows = 5, nrows = 3))
        minlength = float(datapd_first_rows.iloc[0, 0])
        maxlength = float(datapd_first_rows.iloc[0, 1])
        total_volume_length = maxlength - minlength
        datapd_first_rows = np.abs(datapd_first_rows)
        axes_length = datapd_first_rows.sum(axis = 1)
        volume = (axes_length[0] * axes_length[1] * axes_length[2])
        return volume, maxlength, minlength, total_volume_length

    def get_timestep_from_file_name(self, path_to_file):
        """Only reads timestep from the file name."""
        pattern = r"_(\d+)\.txt$"
        timestep = re.search(pattern, path_to_file).group(1)
        return int(timestep)


    def calculate_bond_vectors(self):
        """Returns np.diff in positions except for each 100th particle."""
        # Filter out every nth row
        data = self.datapd 
        midpoint_ridges = self.make_midpoint_ridges()
        bond_vectors_array = (np.diff(data.iloc[:, 1:4], axis = 0))
        polymer_mask = np.ones(len(bond_vectors_array), dtype = bool)
        polymer_mask[99::100] = False
        bond_vectors_array = bond_vectors_array[polymer_mask]
        bond_vectors = data.iloc[:, :5].copy()
        bond_vectors = bond_vectors[bond_vectors.index % 100 != 0]
        bond_vectors.iloc[:, 1:4] = bond_vectors_array
        bond_vectors = bond_vectors.rename(columns = {"xu" : "bondx", "yu" : "bondy", "zu" : "bondz"})
        #Append center of mass to bond_vectors

        wrapped_coordinates = self.wrap_coordinates(data.iloc[:, 1:5])
        #wrapped_coordinates = wrapped_coordinates[wrapped_coordinates.index % 100 != 0]
        midpoints = (wrapped_coordinates.rolling(window = 2).sum().dropna())/2
        midpoints.index = midpoints.index - 1
        midpoints = midpoints[midpoints.index % 100 != 0]
        midpoints.iloc[:, 0] = find_box_id(midpoint_ridges, midpoints.iloc[:, 0].to_numpy())
        midpoints.iloc[:, 1] = find_box_id(midpoint_ridges, midpoints.iloc[:, 1].to_numpy())
        midpoints.iloc[:, 2] = find_box_id(midpoint_ridges, midpoints.iloc[:, 2].to_numpy())
        bond_vectors_2 = pd.concat([bond_vectors, midpoints], axis =1).rename(columns = {"xu" : "xid", "yu" : "yid", "zu" : "zid"})
        #return bond_vectors
        return bond_vectors_2


    def wrap_coordinates(self, data):
        """Converts coordinates to wrapped coordinates. Input can be float or array (float)"""
        return (data - self.minlength) % self.total_volume_length + self.minlength

    def make_midpoint_ridges(self):
        length_array = np.linspace(self.minlength, self.maxlength, self.nridges+1)
        midpoint_ridges = ((length_array + np.roll(length_array, 1))/2)[1:] #Serves as box id 
        return midpoint_ridges

    def assign_center_of_mass(self, nridges = 33):
        """Loops over all polymers to assign center of mass 
        Also assignes a box id to each polymer group
        Takes about 110 seconds over dataset 720k big"""

        data = self.datapd
        midpoint_ridges = self.make_midpoint_ridges()
        box_length =  (self.total_volume_length)/nridges
        df_com = pd.DataFrame(np.zeros([data.shape[0], 3]), columns = ["xid", "yid", "zid"], index = data.index)
        local_volume = self.volume/self.nridges**3
        wrapped_coordinates = self.wrap_coordinates(data.iloc[:, 1:5])
        df_com.iloc[:, 0] = find_box_id(midpoint_ridges, wrapped_coordinates.iloc[:, 0].to_numpy())
        df_com.iloc[:, 1] = find_box_id(midpoint_ridges, wrapped_coordinates.iloc[:, 1].to_numpy())
        df_com.iloc[:, 2] = find_box_id(midpoint_ridges, wrapped_coordinates.iloc[:, 2].to_numpy())

    
        data = pd.concat([data, df_com], axis=1)
        return data, local_volume, df_com


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

    def get_nematic_vector_4(self, save_ev: bool = False, save_string = None):
        data = self.datapd
        # Prepare masks of all possible combinations 
        data = data[data.index % 100 != 0] # Filter out all last monomers as they do not have a bond vector per definiton
        df_cryst = pd.DataFrame(np.zeros([self.combinations.shape[0], 7]), columns = ["xid", "yid", "zid", "cryst_bool", "x_ev", "y_ev", "z_ev"])
        #df_labdas = pd.DataFrame(np.zeros([self.combinations.shape[0], 3]), columns = ["labda_1", "labda_2", "labda_3"])
        df_cryst.iloc[:, :3] = self.combinations
        for t in tqdm(range(0, len(self.combinations))):
            combination = self.combinations[t] #Selects a random box 
            #print(combination)
            #subset = data[(data['xid'] == combination[0]) & (data['yid'] == combination[1]) & (data['zid'] == combination[2])]
            subset = filter_out_subset(data, combination)
            if subset.empty == False:
                # Get index molecules 
                indexes = subset.index
                subset_bond_vectors = self.bond_vectors.loc[indexes]
                order_param, order_ev, labda, ev = calc_nematic_tensor_2(subset_bond_vectors.iloc[:, 1:4].to_numpy())
                df_cryst.iloc[t,3] = order_param
                df_cryst.iloc[t,4:7] = order_ev
                #df_labdas.iloc[t, :] = labda
        #self.df_labdas = df_labdas
        self.fraction_crystallinity, _ = fraction_crystallinity(df_cryst.iloc[:,3])
        self.df_cryst = df_cryst
        #print(self.df_labdas)
        #print(self.fraction_crystallinity)
        if save_ev == True:
            #df_labdas.to_csv("10e5_T1_debug_labdas.txt", sep = " ", mode = "w")
            df_cryst.to_csv("%s" %save_string, sep = " ", mode = "w")
        return self.fraction_crystallinity

    def get_nematic_vector_5(self, save_string = None, cryst_cutoff = 0.8):
        data = self.datapd
        data = data[data.index % 100 != 0] # Filter out all last monomers as they do not have a bond vector per definiton
        #self.df_cryst = nematic_vector_loop(data, self.bond_vectors)
        self.df_cryst = nematic_vector_loop_2(self.bond_vectors)
        self.fraction_crystallinity, _ = fraction_crystallinity(self.df_cryst.iloc[:,3], cutoff= cryst_cutoff)
        if isinstance(save_string, str):
            self.df_cryst.to_csv("%s" %save_string, sep = " ", mode = "w")
        return self.fraction_crystallinity



class polymer():
    def __init__(self, path_to_file):
        self.atom_coords = atom_coords(path_to_file)
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
            subset = self.atom_coords.bond_vectors[(self.atom_coords.bond_vectors["mol_id"] == i+1)]
            dist = np.sum(subset.iloc[:, 1:4])
            df_end_end_length.iloc[i] = (dist.iloc[0]* dist.iloc[0] + dist.iloc[1] * dist.iloc[1] + dist.iloc[2] * dist.iloc[2])
        end_to_end_length = (np.sum(df_end_end_length.to_numpy())/self.atom_coords.no_polymers)
        end_end_distance_normalised = np.sqrt(df_end_end_length.iloc[:])
        self.results.end_to_end_distribution = (end_end_distance_normalised)
        self.results.end_to_end_length = end_to_end_length
        self.results.mean_squared_end_to_end = np.sqrt(end_to_end_length)
        print("mean end-to-end is %f" %self.results.mean_squared_end_to_end)

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
        self.results.gyration_radius_distribution = np.sqrt(df_gyration_radius.iloc[:, 3])
        self.results.gyration_radius = df_gyration_radius
        self.results.mean_gyration_radius = np.mean(df_gyration_radius["gyration_radius"]) #Ensemble average
        print("mean gyration is %f" %np.sqrt(self.results.mean_gyration_radius))


    def bond_bond_correlation(self, show_plot = False):
        #For PVA-100: 7200 x99 x 99 upper triangular matrix, axis-0: j, axis-1: j + 1
        bond_bond_correlation_array = np.zeros([self.atom_coords.no_polymers, self.atom_coords.polymer_length-1, self.atom_coords.polymer_length-1])
        for i in range(0, self.atom_coords.no_polymers):
        #for i in range(0,1):
            subset = self.atom_coords.bond_vectors[(self.atom_coords.bond_vectors["mol_id"] == i + 1)].to_numpy()[:, 1:4]
            #Normalise bond vectors 
            subset = subset/np.linalg.norm(subset, axis = 1, keepdims = True)
            #print(subset)
            bond_bond_corr_array_polymer = np.zeros([self.atom_coords.polymer_length-1, self.atom_coords.polymer_length-1])
            box_algos_lib.bond_bond_correlation(subset, bond_bond_corr_array_polymer, subset.shape[0], subset.shape[1])

            bond_bond_correlation_array[i,:,:] = bond_bond_corr_array_polymer
        bond_correlation_average = np.mean(bond_bond_correlation_array, axis = 0)
        # print(bond_correlation_average)
        #np.savetxt("debug_bond_bond_average.txt", bond_correlation_average)
        #bond_correlation_average = np.loadtxt("debug_bond_bond_average.txt")
        diag_means = np.zeros(self.atom_coords.polymer_length-1)
        for i in range(0, self.atom_coords.polymer_length -1):
        #     #print(np.mean(np.diagonal(bond_correlation_average, offset = i)))
            diag_means[i] = np.mean(np.diagonal(bond_correlation_average, offset = i))
        self.results.bond_bond_correlation = diag_means
        self.atom_coords.positions = np.arange(1, 100)
        # plt.scatter(positions, (diag_means))
        # plt.xlabel("n")
        # plt.ylabel(r"cos\theta(n)")
        # plt.title("Distribution of bond-bond correlations")
        #plt.show()

    def bond_distribution(self):
        #self.atom_coords.get_nematic_vector_5()
        print(self.atom_coords.bond_vectors)
        counts = (
            self.atom_coords.bond_vectors
            .groupby(["xid", "yid", "zid"])
            .size()                      # count rows per combination
            .reset_index(name="count")   # put counts in a column
        )

        dist = counts["count"].value_counts().sort_index()
        print(dist)




    def get_density_dist_per_box(self, nridges = 33):
        """Uses new method with combinations to calculate local density (i.e. density per box)"""
        data = self.atom_coords.datapd
        data_indexed = data.reset_index().rename({"index": "atom_id"}).set_index(["xid", "yid", "zid"]).sort_index()
        #print(data_indexed)
        data_grouped = data_indexed.groupby(level=["xid", "yid", "zid"])
        no_points_per_box = data_grouped.size()
        print(no_points_per_box)
        local_densities = no_points_per_box / self.atom_coords.local_volume
        print(data.shape[0]/self.atom_coords.volume)
        #no_points_per_box.to_csv("local_densities.txt")
        self.results.avg_local_density = np.mean(local_densities)
        return local_densities

    def get_density_dist(self, nridges = None):
        if nridges is None:
            nridges = self.atom_coords.nridges
        _, _, boxes = self.atom_coords.assign_center_of_mass(nridges = nridges)
        boxes = boxes.set_index(["xid", "yid", "zid"]).sort_index()
        data_grouped = boxes.groupby(level=["xid", "yid", "zid"])
        no_points_per_box = data_grouped.size()
        return no_points_per_box/self.atom_coords.local_volume




    def read_cryst(self, location):
        self.df_cryst = pd.read_csv(location, sep = " ", header = None, skiprows = 1, index_col = False).iloc[:, 1:]
        self.df_cryst.columns = ["xid", "yid", "zid", "cryst_bool", "x_ev", "y_ev", "z_ev"]
        self.frac_cryst, _ = fraction_crystallinity(self.df_cryst.iloc[:, 3])
        return 0;

    def merge_boxes(self, ndot_cutoff = 0.97, nridges = 33, cryst_cutoff = 0.8, save = False, print_results: bool = False):
        max_labels = int(nridges**3)
        #max_labels = int(200)
        try:
            self.df_cryst
        except(AttributeError):
            print(self.atom_coords.path_to_file)
            self.read_cryst()
        cryst_array = self.df_cryst.to_numpy()
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
        hk_lib.hoshen_kopelman_crystallisation(cryst_array_c, rows, cols, nridges, cryst_cutoff, ndot_cutoff, label_matrix, new_label_matrix, labels)

        #Return label matrix to c 

        for i in range(nridges):
            for j in range(nridges):
                for k in range(nridges):
                    #label_matrix_np[i,j,k] = label_matrix[j][i][k]
                    label_matrix_np[i,j,k] = label_matrix[i][j][k]



        label_matrix = label_matrix_np
        size = nridges
        #print(labels[labels != 0])
        #print(labels)
        labels2 = labels[labels != 0]

        unique_values, counts = np.unique(label_matrix, return_counts=True) #Labels and how much each label occurs
        # print(unique_values, counts)
        # print(np.mean(counts[1:]))
        # for value, count in zip(unique_values, counts):
        #     if count == 1:
        #         #print(f"count 1 :   {value}: {count}")
        #         pass
        #     else:
        #         print(f"{value}: {count}")

        total_number_merged_clusters = counts[counts > 1]
        self.results.total_number_clusters = total_number_merged_clusters.size
        self.results.total_number_independent_clusters = unique_values.size-1
        self.results.mean_cluster_size = np.mean(counts[1:])
        self.results.total_number_crystalline_grid_elements = np.sum(counts[1:])
        self.results.fraction_crystallinity, no_crystalline_elements = fraction_crystallinity(self.df_cryst.iloc[:,3])
        if print_results == True:
            print("total number clusters w/ >= 2 elements: %i" %(self.results.total_number_clusters))
            print("total number independent crystalline domains: %i" %(self.results.total_number_independent_clusters))
            print("average cluster size crystalline domains: %f" %(self.results.mean_cluster_size))
            print("total number crystalline/all grid elements: %i/%i -> cryt_frac = %f" 
                %(self.results.total_number_crystalline_grid_elements,self.atom_coords.nridges**3, 
                    self.results.total_number_crystalline_grid_elements/self.atom_coords.nridges**3))
        #print("current crystallinity w/o h-k algo: %i/%i -> %f" %(no_crystalline_elements, self.atom_coords.nridges**3, fraction_of_crystallinity))
            print(" ")
        #print(label_matrix)

        if save is True:
            I, J, K = label_matrix.shape

            rows = []
            for k in range(K):        # k: slowest
                for i in range(I):    # i: middle
                    for j in range(J):# j: fastest
                        l = label_matrix[i, j, k]
                        rows.append((i, j, k, l))


            df = pd.DataFrame(rows, columns=["x", "y", "z", "label"])
            df.to_csv("hk_label_matrix.txt", sep = " ", index = True)

        #np.save("hk_label_matrix.npy", label_matrix)
        return label_matrix;

def bin_label_matrix(label_matrix, bins: np.array):
    """Function to bin the cluster size into a histogram"""
    unique_values, counts = np.unique(label_matrix, return_counts=True)

    print(unique_values, counts)
    hist, bin_edges = np.histogram(counts[1:], bins)



    
def check_labels(label_matrix, nridges, i,j,k):
    """Helper functuion to check if label_matrix is equal to one of its neighbours. If not, print current value of label_matrix and all neighbours"""
    if label_matrix[i,j,k] != label_matrix[(i -1)% nridges,j,k] or label_matrix[i,j,k] != label_matrix[(i+1)%nridges, j, k]:
        if label_matrix[i,j,k] != label_matrix[i,(j-1)% nridges,k] or label_matrix[i,j,k] != label_matrix[i, (j+1) % nridges, k]:
            if label_matrix[i,j,k] != label_matrix[i,j,(k-1) % nridges] or label_matrix[i,j,k] != label_matrix[i, j, (k+1) % nridges]:
                print("position %i %i %i, current label %i \nneighbours in -x %i and +x %i" %(i,j,k,label_matrix[i,j,k], label_matrix[(i -1)% nridges, j,k], label_matrix[(i +1)% nridges, j,k]))
                print("neighbours in -y %i, +y %i" %(label_matrix[i, (j - 1) %nridges, k], label_matrix[i, (j + 1) %nridges, k]))
                print("neighbours in -z direction %i, +z %i \n" %(label_matrix[i, j, (k - 1)%nridges], label_matrix[i, j, (k + 1)%nridges]))


    # if label_matrix[i,j,k] != label_matrix[(i-1)% nridges,j,k]:
    #     if label_matrix[i,j,k] != label_matrix[(i+1)%nridges, j, k]:  
    #         if label_matrix[i,j,k] != label_matrix[i,(j-1)% nridges,k]: 
    #             if label_matrix[i,j,k] != label_matrix[i, (j+1) % nridges, k]:
    #                 if label_matrix[i,j,k] != label_matrix[i,j,(k-1) % nridges]:
    #                     if label_matrix[i,j,k] != label_matrix[i, j, (k+1) % nridges]:

class results(object):
    def __init__(self):
        return

    # Include plotting functions here 

def plot_hk_matrix_2d(polymer):
    data = polymer.merge_boxes(print_results = True)
    Nx, Ny, Nz = data.shape

    # ---- GLOBAL COLOUR SCALE (same for all k) ----
    global_max = int(data.max())
    n_labels = max(global_max, 0)

    # Base colours from viridis
    base_colors = plt.cm.viridis(np.linspace(0, 1, n_labels + 1))
    # Make label 0 stand out (e.g. bright red)
    base_colors[0] = np.array([1.0, 0.0, 0.0, 1.0])  # RGBA

    cmap = ListedColormap(base_colors)
    bounds = np.arange(-0.5, n_labels + 1.5, 1)
    norm = BoundaryNorm(bounds, cmap.N)

    # ----------------------------------------------
    for k in range(Nz):
        slice_k = data[:, :, k]

        # Scale figure size with grid size so each cell is big enough
        cell_size = 0.3  # inches per cell (adjust up/down)
        fig_width = max(4, Ny * cell_size)
        fig_height = max(4, Nx * cell_size)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # use the same cmap/norm for every slice
        im = ax.imshow(slice_k, cmap=cmap, norm=norm, origin='lower')

        ax.set_title(
            f"Cluster labels, PVA-100, T = 0.5, Tdot = 10**(-5), "
            f"T_init = 1.0, z = {k}"
        )
        ax.set_xlabel("y")
        ax.set_ylabel("x")

        ax.set_xlim(-0.5, Ny - 0.5)
        ax.set_ylim(-0.5, Nx - 0.5)

        # --- Control which cells get labels ---
        label_only_nonzero = False   # set True if you only want val != 0
        step_x = 1                   # label every 'step_x' cells in x
        step_y = 1                   # label every 'step_y' cells in y

        for x in range(0, Nx, step_x):
            for y in range(0, Ny, step_y):
                val = slice_k[x, y]
                if label_only_nonzero and val == 0:
                    continue
                ax.text(
                    y, x, str(int(val)),
                    ha='center', va='center',
                    color='white' if val < global_max / 2 else 'black',
                    fontsize=6,
                )

        cbar = fig.colorbar(im, ax=ax, label='cluster label')
        cbar.set_ticks(np.arange(0, n_labels + 1))

        fig.tight_layout()
        plt.savefig(f"hk_debug/PVA_100_t=20_zlayer_{k}.pdf")
        # plt.show()
        plt.close()


def plot_hk_matrix(polymer):
    label_matrix = polymer.merge_boxes()


    Nx, Ny, Nz = label_matrix.shape

    # Create coordinate grids for x, y, z (index positions)
    x_idx, y_idx, z_idx = np.indices((Nx, Ny, Nz))

    # Flatten to 1D
    x = x_idx.ravel()
    y = y_idx.ravel()
    z = z_idx.ravel()
    l = label_matrix.ravel()      # values

    mask = l > 0
    x = x[mask]
    y = y[mask]
    z = z[mask]
    l = l[mask]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(x, y, z, c=l, cmap='viridis', s=10)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    fig.colorbar(sc, ax=ax, label='l')

    plt.show()



def main():
    """Testing function"""

    #first_timestep_e5 = polymer("../../data/pva-100/quick_quench/equil_t_085_tdot_e-3_run2_goodrun_time_0.txt")
    slow_quench_e5_time_100e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
    # slow_quench_e5_time_95e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_9500000.txt")
    # slow_quench_e5_time_90e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_9000000.txt")
    # slow_quench_e5_time_85e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_8500000.txt")
    # slow_quench_e5_time_80e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_8000000.txt")


    slow_quench_e5_time_100e5.atom_coords.get_nematic_vector_5("10e5_T05_last_timestep_boxes_ev.txt")
    # slow_quench_e5_time_95e5.atom_coords.get_nematic_vector_5("10e5_T05_timestep_boxes_ev_time_95e5.txt")
    # slow_quench_e5_time_90e5.atom_coords.get_nematic_vector_5("10e5_T05_timestep_boxes_ev_time_90e5.txt")
    # slow_quench_e5_time_85e5.atom_coords.get_nematic_vector_5("10e5_T05_timestep_boxes_ev_time_85e5.txt")
    # slow_quench_e5_time_80e5.atom_coords.get_nematic_vector_5("10e5_T05_timestep_boxes_ev_time_80e5.txt")
    slow_quench_e5_time_100e5.read_cryst("hk_debug/10e5_T05_last_timestep_boxes_ev.txt")
    # slow_quench_e5_time_95e5.read_cryst("hk_debug/10e5_T05_timestep_boxes_ev_time_95e5.txt")
    # slow_quench_e5_time_90e5.read_cryst("hk_debug/10e5_T05_timestep_boxes_ev_time_90e5.txt")
    # slow_quench_e5_time_85e5.read_cryst("hk_debug/10e5_T05_timestep_boxes_ev_time_85e5.txt")
    # slow_quench_e5_time_80e5.read_cryst("hk_debug/10e5_T05_timestep_boxes_ev_time_80e5.txt")
    #print(first_timestep_e5.atom_coords.combinations)
    #last_timestep_e5.get_density_dist()
    label_matrix = slow_quench_e5_time_100e5.merge_boxes()
    # slow_quench_e5_time_95e5.merge_boxes()
    # slow_quench_e5_time_90e5.merge_boxes()
    # slow_quench_e5_time_85e5.merge_boxes()
    # slow_quench_e5_time_80e5.merge_boxes()
    bins = np.array([1,2,3,4,5,7,10,20,30,40,50,60,70])
    bin_label_matrix(label_matrix, bins)


    # plot_hk_matrix_2d(slow_quench_e5_time_100e5)
    
    
if __name__ == "__main__":
    main()