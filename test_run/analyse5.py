import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import re 
import pandas as pd 
from tqdm import tqdm
#from find_box_id import *
import ctypes



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
        #ctypes.POINTER(ctypes.c_double), #Output cryst_array [box-id, cluster id]
        np.ctypeslib.ndpointer(dtype = np.double, flags = "C_CONTIGUOUS"),
        ctypes.POINTER(ctypes.c_double), #1d cluster array 
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int, #no. lattice points/boxes
        ctypes.c_float, #cutoff for crystallisation yes/no
        ctypes.c_float #cutoff for ndot product 
    ]
    lib.find_nearest_value.restype = ctypes.POINTER(ctypes.c_int)  # int* return type
    lib.hoshen_kopelman_crystallisation.restype = None
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
    return max_labda, max_ev, order_param

class atom_coords:

    """Used to read in files and analyse them"""

    def __init__(self, file_to_path, nridges = 33):
        #Note datapd has following columns: 
        self.datapd = self.prepare_position_data(file_to_path)
        self.n_atoms = len(self.datapd.index)
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
        data.iloc[:, 2:5] = (data.iloc[:, 2:5] - self.minlength) % self.total_volume_length + self.minlength
        # data.rename(columns={"x": "xu", "y": "yu", "z": "zu"})
        df_com.iloc[:, 0] = find_box_id(self.midpoint_ridges, data.iloc[:, 1].to_numpy())
        df_com.iloc[:, 1] = find_box_id(self.midpoint_ridges, data.iloc[:, 2].to_numpy())
        df_com.iloc[:, 2] = find_box_id(self.midpoint_ridges, data.iloc[:, 3].to_numpy())
        # for i in range(0, data.shape[0]): #Below is an all-python approach 
        #     df_com.iloc[i, 0] = find_nearest(self.midpoint_ridges, data.iloc[i, 1])
        #     df_com.iloc[i, 1] = find_nearest(self.midpoint_ridges, data.iloc[i, 2])
        #     df_com.iloc[i, 2] = find_nearest(self.midpoint_ridges, data.iloc[i, 3])

            
        data = pd.concat([data, df_com], axis=1)
        #data = data.iloc[:-1, :]
        #data.to_csv("data_test_ctypes.txt", sep = " ", mode = "w")
        return data


    def generate_box_list(self, nridges = 33):
        numbers = np.arange(0, nridges)  # 0 to 32 inclusive
        self.combinations = np.array(np.meshgrid(numbers, numbers, numbers)).T.reshape(-1, 3)
        return self.combinations

    def end_to_end_distance(self, nridges = 33, show_plot = False, save_plot_string = None):
        # Calculates end-to-end distance of each polymer 
        #print(self.datapd)
        no_polymers = self.datapd.iloc[:, 0].max()
        end_end_length = np.zeros([no_polymers, 3])
        df_end_end_length = pd.DataFrame(np.zeros([no_polymers, 4]), columns = ["xl", "yl", "zl", "end_end_length"]) #xlength etc.
        df_end_end_length.index.name = "mol_id"
        df_end_end_length.index = df_end_end_length.index + 1
        for i in range(0, no_polymers):
            #subset = data[(data['xid'] == combination[0]) & (data['yid'] == combination[1]) & (data['zid'] == combination[2])]
            subset = self.datapd[(self.datapd["mol_id"] == i+1)]
            diff = subset.diff()
            #print(i, diff.sum())
            length_3d = diff.sum()[1:4]
            df_end_end_length.iloc[i, :-1] = length_3d
            df_end_end_length.iloc[i, -1] = np.sqrt(length_3d[0]**2 + length_3d[1]**2 + length_3d[2]**2)
            #print(df_end_end_length)
            #print(i)
        #print(np.mean(df_end_end_length["end_end_length"]))
        self.mean_end_to_end_length = np.mean(df_end_end_length["end_end_length"])
        self.end_to_end_length = df_end_end_length
        values, bins, _ = plt.hist(df_end_end_length.iloc[:, -1], bins = 100)
        plt.vlines(np.mean(df_end_end_length["end_end_length"]), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "mean end-to-end length = %.2f" %self.mean_end_to_end_length)
        plt.legend()
        plt.show()
        return self.end_to_end_length;


    def get_nematic_vector_4(self, nridges = 33, save_ev = False, save_string = None):
        data = self.datapd
        # Prepare masks of all possible combinations 
        data = data[data.index % 100 != 0] # Filter out all last monomers as they do not have a bond vector per definiton
        df_cryst = pd.DataFrame(np.zeros([self.combinations.shape[0], 7]), columns = ["xid", "yid", "zid", "cryst_bool", "x_ev", "y_ev", "z_ev"])
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
                labda, ev, order_param = calc_nematic_tensor_2(subset_bond_vectors.iloc[:, 1:4])
                df_cryst.iloc[t,3] = order_param
                df_cryst.iloc[t,4:7] = ev
        self.fraction_crystallinity = fraction_crystallinity(df_cryst.iloc[:,3])
        self.df_cryst = df_cryst
        print(self.fraction_crystallinity)
        if save_ev == True:
            df_cryst.to_csv("%s" %save_string, sep = " ", mode = "w")

        return self.fraction_crystallinity


    def get_distribution_eigenvalues(self, title, nridges = 33, savestring = None):
        values, bins, _ = plt.hist(self.df_cryst.iloc[:, 3], bins = 1000)
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
        #print(self.df_cryst)
        cryst_array = self.df_cryst.iloc[:, 1:].to_numpy()
        #cryst_array_new = np.zeros([cryst_array.shape[0],cryst_array.shape[1]+1])
        #cryst_array_new[:, :-1] = cryst_array
        #cryst_array = cryst_array_new
        output_cryst = np.zeros([nridges**3, 4]) #Lattice points and cluster id 
        max_no_clusters = nridges**3
        actual_no_clusters = ctypes.c_int(0)
        cluster_id_list = np.zeros(max_no_clusters)
        input_cryst_pointer = cryst_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        #output_cryst_pointer = output_cryst.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        output_cluster_pointer = cluster_id_list.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        lib.hoshen_kopelman_crystallisation(input_cryst_pointer, cryst_array.shape[0], cryst_array.shape[1], 
            output_cryst, output_cluster_pointer, ctypes.byref(actual_no_clusters), nridges, cryst_cutoff, ndot_cutoff)
        #cluster_id_array = np.ctypeslib.as_array(out, shape = (cryst_array.shape[0], 4))
        unique_values, counts = np.unique(output_cryst[:, -1], return_counts=True)
        print(counts)
        # plt.bar(unique_values[1:], counts[1:], color='skyblue', edgecolor='black', alpha=0.7)
        # plt.xlabel("cluster-id")
        # plt.savefig("cluster_prelim_firstrun.pdf")
        # plt.show()



            

    def get_density_dist(self, nridges = 33):
        """Uses new method with combinations to calculate local density (i.e. density per box)"""
        print(self.combinations)
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
        title = r"Distribution of eigenvalues at $T = %s$, $\dot{T} = 10^7$" %temps[t]
        savestring = "plots/eigenvalue_dist_pva_100_t_%s" %temps[t]
        atom_coords.get_distribution_eigenvalues(title = title, savestring = savestring)
    return 0;

lib = c_lib_init()


# list_atom_coords_cooling = get_list_atom_coords("../../data/pva-100/cooling_tdot_e-5_time", 21, endtime= 1e7)
# plot_multiple_dists_eigenvalues(list_atom_coords_cooling)
#list_atom_coords_heating = get_list_atom_coords("../../data/pva-100/genua_heating_100_tmin_0.5_ttime_10e7",21, endtime = 1e7)
#plot_order_param(list_atom_coords_cooling, "Crystallinity vs temperature, cooling process, Tdot = 10e-5", savestring = "test_wholebox_frac_cryst_heating_100_tmin_0.5_ttime_10e7")

last_timestep_e5 = atom_coords("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
# mid_timestep_e5 = atom_coords("../../data/pva-100/cooling_tdot_e-5_time_4000000.txt")
# last_timestep_e5.calc_rdf()

#mid_timestep_e5.get_nematic_vector_4()
# mid_timestep_e5.get_distribution_eigenvalues(r"Distribution of eigenvalues at $T = 0.8$, $\dot{T} = 10^7$")
last_timestep_e5.end_to_end_distance()
#last_timestep_e5.get_nematic_vector_4()#save_ev = True, save_string = "10e5_debug_cryst.txt")
# last_timestep_e5.get_distribution_eigenvalues(r"Distribution of eigenvalues at $T = 0.5$, $\dot{T} = 10^7$")
#last_timestep_e5.read_cryst("10e5_debug_cryst.txt")
#last_timestep_e5.merge_boxes()
# #last_timestep_e5.get_density_dist()
# #plot_density_dist(last_timestep_e5, "Distribution of local densities at T = 0.5, tdot 10e-5")
# plot_volume_line(list_atom_coords_cooling, "Volume per monomer as function of temperature, PVA-100", "volume_monomer_tdot_e-5.pdf")

# list_different_tdot_t_08 = get_crystallinity_tdots("../../data/pva-100/cooling_tdot", 4000, np.array([2,3, 4, 5]))
# list_different_tdot_t_05 = get_crystallinity_tdots("../../data/pva-100/cooling_tdot", 10000, np.array([2,3, 4, 5]))

# plt.scatter(list_different_tdot_t_08[0, :], list_different_tdot_t_08[1, :], label = "T = 0.8")
# plt.scatter(list_different_tdot_t_05[0, :], list_different_tdot_t_05[1, :], label = "T = 0.5")
# plt.title("Crystallisation as function of cooling rate")
# plt.xlabel("cooling rate")
# plt.ylabel("crystallisation")
# plt.legend()
# plt.xscale("log")
# #plt.yscale("log")
# plt.savefig("cryst_tdot.pdf")
# plt.show()




