import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import re 
import pandas as pd 


def calc_nematic_tensor_2(array):
    """Calculation for the nematic tensor of a local box. NB this is not used yet in the analysis."""
    array_length = array.shape[0]
    array = array / np.linalg.norm(array, axis = 1, keepdims = True)
    Q = np.zeros([3,3])

    for vector in array:
        outer = np.outer(vector, vector)
        Q = Q + 1.5 * outer

    Q = Q / array_length - 0.5 * np.eye(3)

    labda, ev = np.linalg.eigh(Q)
    max_labda = np.max(labda)
    max_ev = ev[:, np.argmax(labda)]
    #print(labda, ev)

    return max_labda, max_ev

class atom_coords:

    """Used to read in files and analyse them"""

    def __init__(self, file_to_path):
        self.datapd = self.prepare_position_data(file_to_path)
        self.get_volume_box(file_to_path)
        self.get_timestep_from_file_name(file_to_path)
        #self.box_atom_list = divide_into_box(self.datapd)
        #self.check_box_atom_list_exist()
        
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

    # def divide_into_box(self, box_length = 5.0):
    #     #Loop over entire box to calculate the limits. Returns groups of atoms located within [n] numbers of mini-boxes. 
    #     data = self.datapd.iloc[:, 2:]
    #     #Set up box ranges 
    #     box_size = (box_length * box_length *  box_length )
    #     nridges = int(np.ceil((self.maxlength - self.minlength) / box_length))
    #     print(self.minlength, self.maxlength)

    #     i = 0; j = 0; k = 0 
    #     xpos = self.minlength; ypos = self.minlength; zpos = self.minlength
    #     boxes = []
    #     for i in range(0, nridges):
    #         xselect = data[(xpos <= data['xu']) & (data['xu']<= (xpos + box_length))]
    #         xpos = xpos + box_length
    #         for j in range(0, nridges):
    #             yselect = xselect[(ypos <= data['yu']) & (data['yu']<= (ypos + box_length))]
    #             ypos = ypos + box_length
    #             for k in range(0, nridges):
    #                 zselect = yselect[(zpos <= data['zu']) & (data['zu']<= (zpos + box_length))]
    #                 zpos = zpos + box_length
    #                 if len(zselect) > 1:
    #                     #print(zselect)
    #                     boxes.append((zselect))

    #             zpos = self.minlength
    #         ypos = self.minlength
    #         #print(xpos, ypos, zpos)

    #     self.box_atom_list = boxes
    #     self.box_size = box_size
    #     return boxes, box_size

    # def divide_into_box(self, box_length=6.0):
    #     #TODO change this to be a function of nridges and calculate box length for this 
    #     data = self.datapd.iloc[:, 2:]
    #     data = data % (self.maxlength - self.minlength)
    #     box_size = box_length ** 3
    #     #nridges = int(np.ceil((self.maxlength - self.minlength) / box_length))
    #     nridges = 11
    #     boxes = []

    #     for i in range(nridges):
    #         xpos_min = self.minlength + i * box_length
    #         xpos_max = xpos_min + box_length
    #         xselect = data[(xselect := data['xu'].between(xpos_min, xpos_max))]
            
    #         for j in range(nridges):
    #             ypos_min = self.minlength + j * box_length
    #             ypos_max = ypos_min + box_length
    #             yselect = xselect[xselect['yu'].between(ypos_min, ypos_max)]

    #             for k in range(nridges):
    #                 zpos_min = self.minlength + k * box_length
    #                 zpos_max = zpos_min + box_length
    #                 zselect = yselect[yselect['zu'].between(zpos_min, zpos_max)]

    #                 if len(zselect) > 1:
    #                     boxes.append(zselect)

    #     self.box_atom_list = boxes
    #     self.box_size = box_size
    #     return boxes, box_size

    def divide_into_box(self, nridges=33):
        data = self.datapd.iloc[:, 2:]
        total_volume_length = self.maxlength - self.minlength
        data = data % (total_volume_length)
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
        i = 0
        for array in self.box_atom_list:
            labda, ev = calc_nematic_tensor_2(array.to_numpy())
            ev_list[i, :] = ev
            labda_list[i] = labda
            i = i + 1
        self.nematic_vector_list = ev_list
        print((ev_list))
        return self.nematic_vector_list


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



list_atom_coords_cooling = get_list_atom_coords("pva-100/run_pos/genua_cooling_100_tmin_0.5_ttime_10e7", 21)
#plot_volume_line(list_atom_coords_cooling, "Volume per momomer vs temperature, cooling process", savestring = "genua_cooling_100_tmin_0.5_ttime_10e7.pdf")

list_atom_coords_heating = get_list_atom_coords("pva-100/run_pos/genua_heating_100_tmin_0.5_ttime_10e7", 21)
#plot_volume_line(list_atom_coords_heating, "Volume per monomer vs temperature, heating process", savestring = "genua_heating_100_tmin_0.5_ttime_10e7.pdf", starttemp = 0.5, endtemp = 1.0)



# cooling_pva_100_0 = atom_coords("pva-100/run_pos/genua_cooling_100_tmin_0.5_ttime_10e7_0.txt")
# cooling_pva_100_50000000 = atom_coords("pva-100/run_pos/genua_cooling_100_tmin_0.5_ttime_10e7_5000000.txt")
# cooling_pva_100_100000000 = atom_coords("pva-100/run_pos/genua_cooling_100_tmin_0.5_ttime_10e7_10000000.txt")

# heating_pva_100_0 = atom_coords("pva-100/run_pos/genua_heating_100_tmin_0.5_ttime_10e7_0.txt")
# heating_pva_100_50000000 = atom_coords("pva-100/run_pos/genua_heating_100_tmin_0.5_ttime_10e7_5000000.txt")
# heating_pva_100_100000000 = atom_coords("pva-100/run_pos/genua_heating_100_tmin_0.5_ttime_10e7_10000000.txt")


# density_dist_pva_100_100000000 = cooling_pva_100_100000000.density_dist()
# density_dist_pva_100_50000000 = cooling_pva_100_50000000.density_dist()
# density_dist_pva_100_0 = cooling_pva_100_0.density_dist()


# heat_density_dist_pva_100_0  = heating_pva_100_0.density_dist()
# heat_density_dist_pva_100_50000000  = heating_pva_100_50000000.density_dist()
# heat_density_dist_pva_100_100000000  = heating_pva_100_100000000.density_dist()

# plot_density_dist(density_dist_pva_100_0, "Local density per cube, cooling cyclus, time 0, temperature 1.0")
# plot_density_dist(density_dist_pva_100_50000000, "Local density per cube, cooling cyclus, time 5e7, temperature 0.75")
# plot_density_dist(density_dist_pva_100_100000000, "Local density per cube, cooling cyclus, time 10e7, temperature 0.5")

# plot_density_dist(heat_density_dist_pva_100_0, "Local density per cube, heating cyclus, time 0, temperature 0.5")
# plot_density_dist(heat_density_dist_pva_100_50000000, "Local density per cube, heating cyclus, time 5e7, temperature 0.75")
# plot_density_dist(heat_density_dist_pva_100_100000000, "Local density per cube, heating cyclus, time 10e7, temperature 1.0")


# Testing for nematic vector

heating_pva_100_50000000 = atom_coords("pva-100/run_pos/genua_heating_100_tmin_0.5_ttime_10e7_0.txt")
heating_pva_100_50000000.get_nematic_vector()

