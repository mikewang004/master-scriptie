import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import re 
import pandas as pd 
from tqdm import tqdm


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    #return array[idx]
    return idx #Returns index instead of value


def calc_nematic_tensor_2(array):
    """Calculation for the nematic tensor of a local box. NB this is not used yet in the analysis."""
    array_length = array.shape[0]
    array = array / np.linalg.norm(array, axis = 1, keepdims = True)
    Q = np.zeros([3,3])
    # for vector in array:
    #     outer = np.outer(vector, vector)
    #     Q = Q + 1.5 * outer
    #Q = np.einsum('ni,nj->nij', array, array)

    outer = (np.einsum('ni,nj->ij', array, array)) / array_length
    #Q =  np.mean(outer  - (1/3) * np.eye(3), axis = 0) # According to Sommer/Luo Sep 2010
    #Q = 1.5 * np.mean(outer, axis = 0) - 0.5 * np.eye(3) # Sara 2015

    Q = (3/2) * outer - (1/2) * np.eye(3)


    #Q = Q / array_length - 0.5 * np.eye(3)

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
        # Wrap coordinates 
        data.iloc[:, 2:5] = (data.iloc[:, 2:5] - self.minlength) % self.total_volume_length + self.minlength
        # data.rename(columns={"x": "xu", "y": "yu", "z": "zu"})
        for i in range(0, data.shape[0]):
            df_com.iloc[i, 0] = find_nearest(self.midpoint_ridges, data.iloc[i, 1])
            df_com.iloc[i, 1] = find_nearest(self.midpoint_ridges, data.iloc[i, 2])
            df_com.iloc[i, 2] = find_nearest(self.midpoint_ridges, data.iloc[i, 3])

        data = pd.concat([data, df_com], axis=1)
        data.to_csv("data_test.txt", sep = " ", mode = "w")
        return data




    def get_nematic_vector_4(self, nridges = 33):
        data = self.assign_center_of_mass(nridges = nridges)
        # Prepare masks of all possible combinations 
        data = data[data.index % 100 != 0] # Filter out all last monomers as they do not have a bond vector per definiton
        numbers = np.arange(0, nridges)  # 0 to 32 inclusive
        combinations = np.array(np.meshgrid(numbers, numbers, numbers)).T.reshape(-1, 3)
        order_param_list = []
        for t in tqdm(range(0, len(combinations))):
            combination = combinations[t]
            #print(combination)
            subset = data[(data['xid'] == combination[0]) & (data['yid'] == combination[1]) & (data['zid'] == combination[2])]
            if subset.empty == False:
                # Get index molecules 
                indexes = subset.index
                #print(indexes)
                subset_bond_vectors = self.bond_vectors.loc[indexes]
                labda, ev, order_param = calc_nematic_tensor_2(subset_bond_vectors.iloc[:, 1:4])
                order_param_list.append(order_param)
        order_param_array = np.asarray(order_param_list)

        self.fraction_crystallinity = fraction_crystallinity(order_param_array)
        print(self.fraction_crystallinity)
        return self.fraction_crystallinity
            


    def density_dist(self):
        """Calculates density per distribution of monomers in a local box"""
        self.check_box_atom_list_exist()
        density = np.zeros(len(self.box_atom_list))
        i = 0
        for box in self.box_atom_list:
            # Get density 
            #minlength = box.min(axis = 0) #Does not work properly, volume size should be constant but returns whack histogramss
            #maxlength = box.max(axis = 0)
            #size = np.abs(maxlength - minlength)
            #volume = np.sqrt(size[0] * size[1] * size[2])

            #density[i] = len(box)/volume
            density[i] = len(box)/self.box_size #should be this
            #print(density[i])
            i = i + 1
        self.density = density
        return self.density


    def get_nematic_vector(self):
        """function not in use but should calculate nematic vector per local cube"""
        self.check_box_atom_list_exist()
        ev_list = np.zeros([len(self.box_atom_list), 3])
        labda_list = np.zeros(len(self.box_atom_list))
        order_param_list = np.zeros(len(self.box_atom_list))
        i = 0
        for array in self.box_atom_list:
            polymer_id_list = array["mol_id"].unique()
            counts = array["mol_id"].value_counts()
            polymer_vector_list = np.zeros([np.sum(counts) - len(polymer_id_list), 3])
            j = 0
            for polymer_id in polymer_id_list:
                #print(counts.iloc[l])
                #k = counts.iloc[l]
                subset = array[array["mol_id"] == polymer_id].iloc[:, 1:]
                k = subset.shape[0]
                bond_vectors = np.diff(subset, axis = 0)
                polymer_vector_list[j:j+k-1, :] = bond_vectors
                j = j + k- 1

            #print(polymer_vector_list)
            labda, ev, order_param = calc_nematic_tensor_2(polymer_vector_list)
            ev_list[i, :] = ev
            labda_list[i] = labda
            order_param_list[i] = order_param
            i = i + 1
        self.nematic_vector_list = ev_list
        self.fraction_crystallinity = fraction_crystallinity(order_param_list)
        #np.savetxt("nematic_tensor_test.txt", order_param_list)
        return self.fraction_crystallinity



def fraction_crystallinity(data, cutoff = 0.8):
    """Data 1d 1d list/array, defined as data > 0.8 -> crystallinity = 1;
    data <= 0.8 -> crystallinity = 0 as in Sommer/Luo Sep 2010"""

    mask = data > cutoff
    fraction = len(data[mask]) / len(data)
    return fraction

def plot_density_dist(data, title):
    """Returns histogram of local density per cube"""
    plt.hist(data, bins = 25)
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

def plot_order_param(list_atom_coords, title,savestring = None, n_atoms = 720000, starttemp = 1.0, endtemp = 0.5):
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



list_atom_coords_cooling = get_list_atom_coords("../../data/pva-100/cooling_tdot_e-5_time", 21)
plot_order_param(list_atom_coords_cooling, "Crystallinity vs temperature, cooling process", savestring = "test_wholebox_frac_cryst_cooling_100_tmin_0.5_ttime_10e7.pdf")

# last_timestep_e5 = atom_coords("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
# last_timestep_e5.get_nematic_vector_4()


# list_different_tdot_t_08 = get_crystallinity_tdots("../../data/pva-100/cooling_tdot", 40000, np.array([3, 4, 5]))
# list_different_tdot_t_05 = get_crystallinity_tdots("../../data/pva-100/cooling_tdot", 100000, np.array([3, 4, 5]))

plt.scatter(list_different_tdot_t_08[0, :], list_different_tdot_t_08[1, :], label = "T = 0.8")
plt.scatter(list_different_tdot_t_05[0, :], list_different_tdot_t_05[1, :], label = "T = 0.5")
plt.title("Crystallisation as function of cooling rate")
plt.xlabel("cooling rate")
plt.ylabel("crystallisation")
plt.legend()
plt.xscale("log")
#plt.yscale("log")
plt.savefig("cryst_tdot.pdf")
plt.show()
