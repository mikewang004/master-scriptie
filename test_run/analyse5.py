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
    lib.find_nearest_value.restype = ctypes.POINTER(ctypes.c_int)  # int* return type
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
    #Q = 1.5 * np.mean(outer, axis = 0) - 0.5 * np.eye(3) # Sara 2015

    Q = 1.5 * outer - 0.5* np.eye(3) # Sara 2015


    #order_param = np.sqrt(1.5 * np.trace(Q**2)) #Sommer/Luo 2010
    labda, ev = np.linalg.eigh(Q)
    max_labda = np.max(labda)
    max_ev = ev[:, np.argmax(labda)]
    order_param = max_labda #Sara 2015
    return max_labda, max_ev, order_param

class atom_coords:

    """Used to read in files and analyse them"""

    def __init__(self, file_to_path):
        self.datapd = self.prepare_position_data(file_to_path)
        self.n_atoms = len(self.datapd.index)
        combinations = self.generate_box_list()
        self.bond_vectors = self.calculate_bond_vectors()
        self.get_volume_box(file_to_path)
        self.get_timestep_from_file_name(file_to_path)
        #self.box_atom_list = divide_into_box(self.datapd)
        #self.check_box_atom_list_exist()

    def read_in_temperature_scale(self, file_to_nano_slurm):
        pass
        # Should read temperature per timestep from the Slurm file 
        # Returns a 2d array with the temperature per timestep 
        
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
        #data.iloc[:, 1:] = data.iloc[:, 1:] % total_volume_length
        box_length =  (self.total_volume_length)/nridges
        df_com = pd.DataFrame(np.zeros([data.shape[0], 3]), columns = ["xid", "yid", "zid"])
        print(self.total_volume_length, self.minlength, self.maxlength)
        self.local_volume = box_length**3
        # Wrap coordinates 
        data.iloc[:, 2:5] = (data.iloc[:, 2:5] - self.minlength) % self.total_volume_length + self.minlength
        # data.rename(columns={"x": "xu", "y": "yu", "z": "zu"})

        df_com.iloc[:, 0] = find_box_id(self.midpoint_ridges, data.iloc[:, 1].to_numpy())
        df_com.iloc[:, 1] = find_box_id(self.midpoint_ridges, data.iloc[:, 2].to_numpy())
        df_com.iloc[:, 2] = find_box_id(self.midpoint_ridges, data.iloc[:, 3].to_numpy())
        #for i in range(0, data.shape[0]): #Below is an all-python approach 
            #df_com.iloc[i, 0] = find_nearest(self.midpoint_ridges, data.iloc[i, 1])
            #df_com.iloc[i, 1] = find_nearest(self.midpoint_ridges, data.iloc[i, 2])
            #df_com.iloc[i, 2] = find_nearest(self.midpoint_ridges, data.iloc[i, 3])

            
        data = pd.concat([data, df_com], axis=1)
        #data.to_csv("data_test_ctypes.txt", sep = " ", mode = "w")
        return data


    def generate_box_list(self, nridges = 33):
        numbers = np.arange(0, nridges)  # 0 to 32 inclusive
        self.combinations = np.array(np.meshgrid(numbers, numbers, numbers)).T.reshape(-1, 3)
        return self.combinations


    def get_nematic_vector_4(self, nridges = 33, save_ev = False):
        data = self.assign_center_of_mass(nridges = nridges)
        # Prepare masks of all possible combinations 
        data = data[data.index % 100 != 0] # Filter out all last monomers as they do not have a bond vector per definiton
        df_cryst = pd.DataFrame(np.zeros([self.combinations.shape[0], 7]), columns = ["xid", "yid", "zid", "cryst_bool", "x_ev", "y_ev", "z_ev"])
        df_cryst.iloc[:, :3] = self.combinations
        for t in tqdm(range(0, len(self.combinations))):
        #for t in tqdm(range(0, 10)):
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
        return self.fraction_crystallinity

    def merge_boxes(self, ndot_cutoff = 0.97, nridges = 33):
        #while clustering_done == False
        for t in tqdm(range(0, len(self.combinations))):
            combination = self.combinations[t]
            #subset = self.df_cryst[(self.df_cryst['xid'] == combination[0]) & (self.df_cryst['yid'] == combination[1]) & (self.df_cryst['zid'] == combination[2])]
            subset = filter_out_subset(self.df_cryst, combination)
            if subset.empty == False:
                # Check whether n dot n >= cutoff 
                print(subset)
                x_left = (combination + np.array([-1,0,0])) % nridges
                x_right = (combination + np.array([+1,0,0])) % nridges
                y_left = (combination + np.array([0,-1,0])) % nridges
                y_right = (combination + np.array([0,+1,0])) % nridges
                z_left = (combination + np.array([0,0,-1])) % nridges
                z_right = (combination + np.array([0,0,+1])) % nridges

                #nn_left = subset.iloc[4:7] @ filter_out_subset(self.df_cryst, x_left).iloc[4:7]
                print(filter_out_subset(self.df_cryst, x_left))



            

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
        plt.savefig(savestring)
    plt.show()

def plot_order_param(list_atom_coords, title,savestring = None, starttemp = 1.0, endtemp = 0.5):
    """Returns plot of volume as function of temperature"""
    list_order_params = []
    for t in tqdm(range(0, len(list_atom_coords))):
        print(t)
        atom_coords = list_atom_coords[t]
    #for atom_coords in list_atom_coords:
        atom_coords.get_nematic_vector_4()
        list_order_params.append((atom_coords.fraction_crystallinity))
        print(atom_coords.fraction_crystallinity)

    temps = np.linspace(starttemp, endtemp, num = len(list_atom_coords))

    plt.scatter(temps, list_order_params)
    plt.title(title)
    plt.xlabel("Temperature [unitless]")
    plt.ylabel("Fraction of crystallinity")
    if savestring == None:
        pass
    else:
        plt.savefig(savestring)
    plt.show()

lib = c_lib_init()


#list_atom_coords_cooling = get_list_atom_coords("../../data/pva-100/cooling_tdot_e-5_time", 21, endtime= 1e7)
#list_atom_coords_heating = get_list_atom_coords("../../data/pva-100/genua_heating_100_tmin_0.5_ttime_10e7",21, endtime = 1e7)
#plot_order_param(list_atom_coords_heating, "Crystallinity vs temperature, heating process", savestring = "test_wholebox_frac_cryst_heating_100_tmin_0.5_ttime_10e7.pdf")

last_timestep_e5 = atom_coords("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
last_timestep_e5.get_nematic_vector_4()
last_timestep_e5.merge_boxes()
# #last_timestep_e5.get_density_dist()
# #plot_density_dist(last_timestep_e5, "Distribution of local densities at T = 0.5, tdot 10e-5")
# plot_volume_line(list_atom_coords_cooling, "Volume per monomer as function of temperature, PVA-100", "volume_monomer_tdot_e-5.pdf")

# list_different_tdot_t_08 = get_crystallinity_tdots("../../data/pva-100/cooling_tdot", 40000, np.array([3, 4, 5]))
# list_different_tdot_t_05 = get_crystallinity_tdots("../../data/pva-100/cooling_tdot", 100000, np.array([3, 4, 5]))

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




