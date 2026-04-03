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
from nematic_vector import calc_nematic_tensor_2, nematic_vector_loop, nematic_vector_loop_2, compute_Q
from matplotlib.colors import ListedColormap, BoundaryNorm
from hoshenKopelmanInPython import hk_in_python

def wrap_coordinates(column, minlength, total_length):
    return (column - minlength) % total_length + minlength

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

    def __init__(self, path_to_file, cell_length = 2):
        #Read in data 
        self.path_to_file = path_to_file
        self.datapd = self.read_in_file(path_to_file)
        self.volume, self.boxlengths, self.dimensions = self.get_volume_box(path_to_file)
        self.current_timestep = self.get_timestep_from_file_name(path_to_file) 
        self.wrapped_monomers = self.wrap_coordinates_all_data()
        #Calculate box properties 
        self.n_atoms = len(self.datapd.index)
        self.no_polymers = self.datapd["mol_id"].max()
        self.polymer_length = int(self.n_atoms/self.no_polymers)
        self.bond_vectors = self.calculate_bond_vectors()
        self.bond_vectors, self.nridges = self.make_cell_grid()
        self.df_cryst = self.get_nematic_vector_5()
        self.results = results()


    def read_in_file(self, filename):
        datapd = pd.read_csv(filename, sep = " ", header = None, skiprows = 9)
        datapd.columns = ["atom_id", "mol_id", "xu", "yu", "zu"]
        datapd = datapd.sort_values("atom_id")
        datapd = datapd.set_index("atom_id")

        return datapd


    def get_volume_box(self, path_to_file):
        """Calculates volume of the total (!) simulation box."""
        datapd_first_rows = (pd.read_csv(path_to_file, sep = " ", header = None, skiprows = 5, nrows = 3))
        datapd_first_rows = datapd_first_rows.rename({0: "min", 1: "max"}, axis = 1).rename({0: "x", 1: "y", 2:"z"}, axis = 0)
        length = datapd_first_rows.iloc[:, 1] - datapd_first_rows.iloc[:, 0]
        volume = length.loc["x"]*length.loc["y"]*length.loc["z"]

        return volume, length, datapd_first_rows

    def get_timestep_from_file_name(self, path_to_file):
        """Only reads timestep from the file name."""
        pattern = r"_(\d+)\.txt$"
        timestep = re.search(pattern, path_to_file).group(1)
        return int(timestep)

    def wrap_coordinates_all_data(self):
        """Converts coordinates to wrapped coordinates. Input can be float or array (float)"""
        wrapped_coords = self.datapd
        wrapped_coords["xu"] = wrap_coordinates(self.datapd["xu"], self.dimensions.loc["x", "min"], self.boxlengths["x"])
        wrapped_coords["yu"] = wrap_coordinates(self.datapd["yu"], self.dimensions.loc["y", "min"], self.boxlengths["y"])
        wrapped_coords["zu"] = wrap_coordinates(self.datapd["zu"], self.dimensions.loc["z", "min"], self.boxlengths["z"])
        return wrapped_coords

    def calculate_bond_vectors(self):

        shifted = self.datapd.groupby('mol_id')[['xu', 'yu', 'zu']].shift(-1)
        bond_vecs = shifted - self.datapd[['xu', 'yu', 'zu']]

        bond_vecs.columns = ['bx', 'by', 'bz']
        bond_vecs = bond_vecs.dropna()
        return bond_vecs


    def make_cell_grid(self, cell_length = 2):
        nridges = (self.boxlengths/cell_length).astype(int)
        actual_cell_length = self.boxlengths/nridges
        gridx = np.linspace(self.dimensions.loc["x", "min"], self.dimensions.loc["x", "max"], nridges.loc["x"]+1)
        gridy = np.linspace(self.dimensions.loc["y", "min"], self.dimensions.loc["y", "max"], nridges.loc["y"]+1)
        gridz = np.linspace(self.dimensions.loc["z", "min"], self.dimensions.loc["z", "max"], nridges.loc["z"]+1)
        #Shift monomers to start from 0 
        #print(self.bond_vectors)
        nx = pd.cut(self.wrapped_monomers["xu"] + self.bond_vectors["bx"]/2, bins = gridx, right = False, labels = False)
        ny = pd.cut(self.wrapped_monomers["yu"] + self.bond_vectors["by"]/2, bins = gridy, right = False, labels = False)
        nz = pd.cut(self.wrapped_monomers["zu"] + self.bond_vectors["bz"]/2, bins = gridz, right = False, labels = False)
        # nx = pd.cut(self.wrapped_monomers["xu"], bins = gridx, right = False, labels = False)
        # ny = pd.cut(self.wrapped_monomers["yu"], bins = gridy, right = False, labels = False)
        # nz = pd.cut(self.wrapped_monomers["zu"], bins = gridz, right = False, labels = False)
        self.bond_vectors = self.bond_vectors.assign(
            nx=nx,
            ny=ny,
            nz=nz,
        )

        return self.bond_vectors, nridges

    def get_nematic_vector_5(self, save_string = None, cryst_cutoff = 0.8):
        #self.df_cryst = nematic_vector_loop(data, self.bond_vectors)
        self.df_cryst = (
            self.bond_vectors.groupby(['nx', 'ny', 'nz'])
            .apply(compute_Q)
            .reset_index()
        )
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

    def atom_distribution(self):
        wrapped_coords = self.atom_coords.wrap_coordinates(self.atom_coords.datapd.iloc[:, 1:4])
        #wrapped_coords = self.atom_coords.datapd.iloc[:, 1:4]
        print(wrapped_coords)

        coords = wrapped_coords[['xu', 'yu', 'zu']].to_numpy()

        # Define 33 bins between min and max for each axis
        ranges = [
            (wrapped_coords['xu'].min(), wrapped_coords['xu'].max()),
            (wrapped_coords['yu'].min(), wrapped_coords['yu'].max()),
            (wrapped_coords['zu'].min(), wrapped_coords['zu'].max()),
        ]

        # 3D histogram: 33 bins along each of x, y, z
        H, edges = np.histogramdd(
            coords,
            bins=(33, 33, 33),
            range=ranges,
        )

        print("Histogram shape:", H.shape)   # (33, 33, 33)
        print("Total count in histogram:", int(H.sum()))




        H_1d = H.ravel()
        unique_counts, freqs = np.unique(H_1d, return_counts=True)

        for c, f in zip(unique_counts, freqs):
            print(f"count={int(c)} in {int(f)} cells")



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

    def merge_boxes_2(self, ndot_cutoff = 0.97, nridges = 33, cryst_cutoff = 0.8, save = False, print_results: bool = False):
        label_matrix = hk_in_python(self.df_cryst, self.frac_cryst, ndot_cutoff = ndot_cutoff, nridges = nridges, cryst_cutoff = cryst_cutoff)

        unique_values, counts = np.unique(label_matrix, return_counts=True) #Labels and how much each label occurs


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
            print("earlier calculated frac_cryst = %f" %(self.frac_cryst))
        return label_matrix

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


def plot_hk_matrix_2d(polymer, ndot_cutoff=0.98, cryst_threshold=0.8):
    # 3D array of cluster labels
    data = polymer.merge_boxes_2(ndot_cutoff = ndot_cutoff, print_results=True)
    Nx, Ny, Nz = data.shape

    dfc = polymer.df_cryst.copy()

    # Integer indices for the grid
    ix = dfc["xid"].astype(int).to_numpy()
    iy = dfc["yid"].astype(int).to_numpy()
    iz = dfc["zid"].astype(int).to_numpy()

    # ---- Build cryst_bool, x_ev, y_ev, z_ev arrays aligned with data ----
    cryst_vals = np.full((Nx, Ny, Nz), np.nan, dtype=float)
    xev_vals   = np.full((Nx, Ny, Nz), np.nan, dtype=float)
    yev_vals   = np.full((Nx, Ny, Nz), np.nan, dtype=float)
    zev_vals   = np.full((Nx, Ny, Nz), np.nan, dtype=float)

    cryst_vals[ix, iy, iz] = dfc["cryst_bool"].to_numpy()
    xev_vals[ix, iy, iz]   = dfc["x_ev"].to_numpy()
    yev_vals[ix, iy, iz]   = dfc["y_ev"].to_numpy()
    zev_vals[ix, iy, iz]   = dfc["z_ev"].to_numpy()

    # ---- Global colour scale for cluster labels ----
    global_max = int(data.max())
    n_labels = max(global_max, 0)

    base_colors = plt.cm.viridis(np.linspace(0, 1, n_labels + 1))
    base_colors[0] = np.array([1.0, 0.0, 0.0, 1.0])  # RGBA
    cmap = ListedColormap(base_colors)
    bounds = np.arange(-0.5, n_labels + 1.5, 1)
    norm = BoundaryNorm(bounds, cmap.N)

    for k in range(Nz):
        slice_k = data[:, :, k]
        cryst_k = cryst_vals[:, :, k]

        cell_size = 0.8
        fig_width = max(4, Ny * cell_size)
        fig_height = max(4, Nx * cell_size)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        im = ax.imshow(slice_k, cmap=cmap, norm=norm, origin='lower')

        ax.set_title(
            f"Cluster labels, PVA-100, T = 0.88, Tdot = 10**(-3), z = {k}"
        )
        ax.set_xlabel("y")
        ax.set_ylabel("x")
        ax.set_xlim(-0.5, Ny - 0.5)
        ax.set_ylim(-0.5, Nx - 0.5)

        # --- Cell labels: cluster + cryst_bool, plus z-edge values above/below ---
        label_only_nonzero = False
        step_x = 1
        step_y = 1

        for x in range(0, Nx, step_x):
            for y in range(0, Ny, step_y):
                val = slice_k[x, y]
                if label_only_nonzero and val == 0:
                    continue

                cryst_val = cryst_k[x, y]

                # Central cell text: cluster + cryst_bool
                lines = [f"{int(val)}"]
                if not np.isnan(cryst_val):
                    lines.append(f"{cryst_val:.3f}")
                text_str = "\n".join(lines)

                ax.text(
                    y, x, text_str,
                    ha='center', va='center',
                    color='white' if val < global_max / 2 else 'black',
                    fontsize=5,
                )

                # If this cell is "red" (cryst_bool < threshold), skip z-val_edge
                if np.isnan(cryst_val) or cryst_val < cryst_threshold:
                    continue

                # Periodic neighbors in z
                k_above = (k + 1) % Nz   # border "above" current slice
                k_below = (k - 1) % Nz   # border "below" current slice

                # z-edge above: between (x,y,k) and (x,y,k_above)
                val_above = None
                if (not np.isnan(xev_vals[x, y, k]) and
                    not np.isnan(xev_vals[x, y, k_above]) and
                    not np.isnan(yev_vals[x, y, k]) and
                    not np.isnan(yev_vals[x, y, k_above]) and
                    not np.isnan(zev_vals[x, y, k]) and
                    not np.isnan(zev_vals[x, y, k_above])):

                    P2Q_above = (
                        (xev_vals[x, y, k] * xev_vals[x, y, k_above]) +
                        (yev_vals[x, y, k] * yev_vals[x, y, k_above]) +
                        (zev_vals[x, y, k] * zev_vals[x, y, k_above])
                    )
                    val_above = 1.5 * P2Q_above**2 - 0.5

                # z-edge below: between (x,y,k_below) and (x,y,k)
                val_below = None
                if (not np.isnan(xev_vals[x, y, k_below]) and
                    not np.isnan(xev_vals[x, y, k]) and
                    not np.isnan(yev_vals[x, y, k_below]) and
                    not np.isnan(yev_vals[x, y, k]) and
                    not np.isnan(zev_vals[x, y, k_below]) and
                    not np.isnan(zev_vals[x, y, k])):

                    P2Q_below = (
                        (xev_vals[x, y, k_below] * xev_vals[x, y, k]) +
                        (yev_vals[x, y, k_below] * yev_vals[x, y, k]) +
                        (zev_vals[x, y, k_below] * zev_vals[x, y, k])
                    )
                    val_below = 1.5 * P2Q_below**2 - 0.5

                # Z-edge values: separate texts, colored by ndot_cutoff, with z+/z- prefix
                if val_above is not None:
                    color_above = 'blue' if val_above >= ndot_cutoff else 'red'
                    ax.text(
                        y, x + 0.35, f"z+ {val_above:.2f}",
                        ha='center', va='center',
                        color=color_above,
                        fontsize=4,
                    )

                if val_below is not None:
                    color_below = 'blue' if val_below >= ndot_cutoff else 'red'
                    ax.text(
                        y, x - 0.35, f"z- {val_below:.2f}",
                        ha='center', va='center',
                        color=color_below,
                        fontsize=4,
                    )

        # --- Edge labels in x,y directions with periodic BC, colored by ndot_cutoff ---

        edge_step_x = 1
        edge_step_y = 1

        # "Horizontal" borders (between rows): (x,y)–(x_next,y), x_next periodic
        for x in range(0, Nx, edge_step_x):
            x_next = (x + 1) % Nx
            for y in range(0, Ny, edge_step_y):
                c1 = cryst_vals[x, y, k]
                # only require THIS cell not red
                if np.isnan(c1) or c1 < cryst_threshold:
                    continue

                if (np.isnan(xev_vals[x,     y, k]) or np.isnan(xev_vals[x_next, y, k]) or
                    np.isnan(yev_vals[x,     y, k]) or np.isnan(yev_vals[x_next, y, k]) or
                    np.isnan(zev_vals[x,     y, k]) or np.isnan(zev_vals[x_next, y, k])):
                    continue

                P2Q = (
                    (xev_vals[x,     y, k] * xev_vals[x_next, y, k]) +
                    (yev_vals[x,     y, k] * yev_vals[x_next, y, k]) +
                    (zev_vals[x,     y, k] * zev_vals[x_next, y, k])
                )
                val_edge = 1.5 * P2Q**2 - 0.5
                edge_color = 'blue' if val_edge >= ndot_cutoff else 'red'

                # Midpoint in x with special handling for wrap edge (Nx-1 ↔ 0)
                if x == Nx - 1 and x_next == 0:
                    mid_x = Nx - 0.5  # bottom boundary
                else:
                    mid_x = x + 0.5

                mid_y = y
                ax.text(
                    mid_y, mid_x, f"{val_edge:.2f}",
                    ha='center', va='center',
                    color=edge_color,
                    fontsize=4,
                )

        # "Vertical" borders (between columns): (x,y)–(x,y_next), y_next periodic
        for x in range(0, Nx, edge_step_x):
            for y in range(0, Ny, edge_step_y):
                y_next = (y + 1) % Ny

                c1 = cryst_vals[x, y, k]
                if np.isnan(c1) or c1 < cryst_threshold:
                    continue

                if (np.isnan(xev_vals[x, y,     k]) or np.isnan(xev_vals[x, y_next, k]) or
                    np.isnan(yev_vals[x, y,     k]) or np.isnan(yev_vals[x, y_next, k]) or
                    np.isnan(zev_vals[x, y,     k]) or np.isnan(zev_vals[x, y_next, k])):
                    continue

                P2Q = (
                    (xev_vals[x, y,     k] * xev_vals[x, y_next, k]) +
                    (yev_vals[x, y,     k] * yev_vals[x, y_next, k]) +
                    (zev_vals[x, y,     k] * zev_vals[x, y_next, k])
                )
                val_edge = 1.5 * P2Q**2 - 0.5
                edge_color = 'blue' if val_edge >= ndot_cutoff else 'red'

                # Midpoint in y with special handling for wrap edge (Ny-1 ↔ 0)
                if y == Ny - 1 and y_next == 0:
                    mid_y = Ny - 0.5  # right boundary
                else:
                    mid_y = y + 0.5

                mid_x = x
                ax.text(
                    mid_y, mid_x, f"{val_edge:.2f}",
                    ha='center', va='center',
                    color=edge_color,
                    fontsize=4,
                )

        cbar = fig.colorbar(im, ax=ax, label='cluster label')
        cbar.set_ticks(np.arange(0, n_labels + 1))
        fig.tight_layout()
        plt.savefig(f"hk_debug/test_PVA_100_polynumber_20_T088_zlayer_{k}.pdf")
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