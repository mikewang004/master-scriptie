import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import re 
import pandas as pd 

# Check periodic boundary conditions: gebruik unwrapped voor bond vectors 

def calc_nematic_tensor_2(array):
    """Calculation for the nematic tensor of a local box."""
    array_length = array.shape[0]
    array = array / np.linalg.norm(array, axis = 1, keepdims = True)
    #Q = np.zeros([3,3])

    outer = (np.einsum('ni,nj->nij', array, array)) 
    #Q =  np.mean(outer  - (1/3) * np.eye(3), axis = 0) # According to Sommer/Luo Sep 2010
    # Q = 1.5 * np.mean(outer, axis = 0) - 0.5 * np.eye(3) # Sara 2015 This works !


    Q = np.zeros([3,3])
    for i in range(array_length):
        vector = array[i, :]
        outer = np.outer(vector, vector)
        Q = Q + 1.5 * outer - 0.5 * np.eye(3)
    Q = Q / array_length

    #order_param = np.sqrt(1.5 * np.trace(Q**2)) #Sommer/Luo 2010
    labda, ev = np.linalg.eigh(Q)
    print(np.sum(labda))
    max_labda = np.max(labda)
    max_ev = ev[:, np.argmax(labda)]
    order_param = max_labda #Sara 2015
    return max_labda, max_ev, order_param

class atom_coords:

    """Used to read in files and analyse them"""

    def __init__(self, file_to_path):
        self.datapd = self.prepare_position_data(file_to_path)
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


    def divide_into_box(self, nridges=33):
        """Default nridges = 33"""
        data = self.datapd.iloc[:, 1:]
        total_volume_length = self.maxlength - self.minlength
        data = data % (total_volume_length) #For periodic boundary conditions, LAMMPS does not apply this automatically
        box_length = total_volume_length/nridges
        box_size = box_length * box_length * box_length
        boxes = []

        for i in range(nridges):
            xpos_min = self.minlength + i * box_length
            xpos_max = xpos_min + box_length
            xselect = data[(xselect := data['xu'].between(xpos_min, xpos_max))]
            
            for j in range(nridges):
                ypos_min = self.minlength + j * box_length
                ypos_max = ypos_min + box_length
                yselect = xselect[xselect['yu'].between(ypos_min, ypos_max)]

                for k in range(nridges):
                    zpos_min = self.minlength + k * box_length
                    zpos_max = zpos_min + box_length
                    zselect = yselect[yselect['zu'].between(zpos_min, zpos_max)]

                    if len(zselect) > 1:
                        boxes.append(zselect)

        self.box_atom_list = boxes
        self.box_size = box_size
        return boxes, box_size


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
                length = np.linalg.norm(bond_vectors, axis = 1)
                bond_vectors = bond_vectors / length[:, np.newaxis]
                polymer_vector_list[j:j+k-1, :] = bond_vectors
                j = j + k- 1

            labda, ev, order_param = calc_nematic_tensor_2(polymer_vector_list)
            ev_list[i, :] = ev
            labda_list[i] = labda
            order_param_list[i] = order_param
            i = i + 1
        self.nematic_vector_list = ev_list
        self.fraction_crystallinity = fraction_crystallinity(order_param_list)
        #np.savetxt("nematic_tensor_test.txt", order_param_list)
        return self.fraction_crystallinity;



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

def get_cooling_rates_time(shared_name, rates, time):
    """Wrapper to get a list of different cooling rates at the same time
    shared_name: string, rates: np.array, time: int
    rates array should be ordered from fastest to slowest rate"""
    list_atom_coords = []
    for i in range(0, len(rates)):
        list_atom_coords.append(atom_coords("%s-%i_time_%i.txt" %(shared_name, rates[i], time*10**i)))
    return list_atom_coords 


def cooling_rates_nematic_order(list_atom_coords, cooling_rates):
    frac_cryst = np.zeros([2, len(list_atom_coords)])
    frac_cryst[0, :] = 10.0**(-1 * cooling_rates)
    i = 0
    for atom_coords in list_atom_coords:
        frac = atom_coords.get_nematic_vector()
        frac_cryst[1, i] = frac
        i = i + 1
    return frac_cryst

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
    for atom_coords in list_atom_coords:
        atom_coords.get_nematic_vector()
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

# cooling_pva_100_t0 = atom_coords("pva-100/run_pos/genua_cooling_100_tmin_0.5_ttime_10e7_0.txt")
# cooling_pva_100_t10000000 = atom_coords("pva-100/run_pos/genua_cooling_100_tmin_0.5_ttime_10e7_10000000.txt")
# print("Start cooling t0")
# cooling_pva_100_t0.get_nematic_vector()
# print("End cooling t10000000")
# cooling_pva_100_t10000000.get_nematic_vector()


# list_atom_coords_cooling = get_list_atom_coords("pva-100/run_pos/genua_cooling_100_tmin_0.5_ttime_10e7", 21)
# plot_order_param(list_atom_coords_cooling, "Crystallinity vs temperature, cooling process", savestring = "test_wholebox_frac_cryst_cooling_100_tmin_0.5_ttime_10e7.pdf")
cooling_rates = np.array([3,4,5])
list_atom_coords = get_cooling_rates_time("../../data/pva-100/cooling_tdot_e", cooling_rates, 40000)
frac_cryst_08 = cooling_rates_nematic_order(list_atom_coords, cooling_rates)
frac_cryst_05 = cooling_rates_nematic_order(get_cooling_rates_time("../../data/pva-100/cooling_tdot_e", cooling_rates, 100000), cooling_rates)


plt.scatter(frac_cryst_05[0, :], frac_cryst_05[1, :], label = "T = 0.5")
plt.scatter(frac_cryst_08[0, :], frac_cryst_08[1, :], label = "T = 0.8")

plt.title("Crystallinity vs cooling rate")
plt.xlabel("Tdot")
plt.ylabel("crystallinity")
plt.xscale("log")
plt.legend()
plt.savefig("cryst_tdot.pdf")
plt.show()