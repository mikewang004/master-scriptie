from analyse8 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from pathlib import Path
import pandas as pd

cryst_path_prefix = "../crystallisation"


def avrami_eq(t, a, n, b):
    return a * (t**n) + b

def avrami_fit(t, y):
    popt, pcov = sp.optimize.curve_fit(avrami_eq, t, np.log(1- y), maxfev = 100000)
    return popt, pcov



class Simulation:
    """Class for one run, holds multiple polymer objects."""

    def __init__(self, temperature: float, cooling_rate: int, path_slurm_file: str, lammps_file_name_without_time :str, run_id: int = 1,
        no_runs: int = 1, polymer_length: int = 100, home_folder: str = None, cryst_cutoff: float = 0.8, ndot_cutoff: float = 0.97, home_folder_override: bool = False):
        """lammps_file_name_without_time should not contain a time"""
        self.temp = temperature
        self.cooling_rate = cooling_rate
        self.polymer_length = polymer_length
    
        self.template_lammps_file_name = lammps_file_name_without_time
        self.time_temp_array, self.list_lammps_files = self.merge_runs(path_slurm_file, lammps_file_name_without_time, no_runs)
        self.max_temp, self.min_temp = np.max(self.time_temp_array[:, 0]), np.min(self.time_temp_array[:, 0])
        temp_str = "T" + str(self.temp).replace(".", "")
        if home_folder_override == True:
            self.path_to_home_folder = "%s/run_%i" %(home_folder, run_id)
        else:
            self.path_to_home_folder = "%s/e%i/%s/run_%i" %(home_folder, cooling_rate, temp_str, run_id)
        Path(home_folder).mkdir(parents = True, exist_ok = True) #Make home directory if not exists
        self.n_atoms = self.get_system_size() #Constant in the entire simulation
        self.cryst_cutoff, self.ndot_cutoff = cryst_cutoff, ndot_cutoff
        self.no_polymers = int(self.n_atoms/self.polymer_length)
        self.cryst = self.get_crystallisation() #2D array containing time and crystallisation at time
        #self.mean_cluster_length = self.get_mean_cluster_length()
        _ ,self.clusters_geq_2, self.no_clusters, self.mean_cryst_domain_size, self.no_cryst_grid_elements = self.get_mean_cluster_length()

    def get_polymer_by_count(self, count):
        """Gets i-th polymer"""
        if isinstance(count, int):
            #Also get corresponding crystallisation array 
            current_polymer = polymer(self.list_lammps_files[count], polymer_length= self.polymer_length)
            #print(self.time_temp_array[count, 0])
            try:
                current_polymer.read_cryst("%s/nematic_vectors/nem_vectors_time_%i.txt" %(self.path_to_home_folder, self.time_temp_array[count, 0]))
            except FileNotFoundError:
                pass
            current_polymer.timestep = self.time_temp_array[count, 0]
            #print(current_polymer.df_cryst)
            return current_polymer
        elif isinstance(count, np.ndarray):
            polymer_list = []
            for i in count:
                current_polymer = polymer(self.list_lammps_files[count], polymer_length= self.polymer_length)
                current_polymer.read_cryst("%s/nematic_vectors/nem_vectors_time_%i.txt" %(self.path_to_home_folder, self.time_temp_array[count, 1]))
                current_polymer.timestep = self.time_temp_array[count, 0]
                polymer_list.append(current_polymer)
            return polymer_list
        else:
            raise Exception("count must be an int or a numpy.ndarray")

    def get_system_size(self):
        """Returns amount of atoms in a given systen"""
        poly = self.get_polymer_by_count(0)
        return poly.atom_coords.n_atoms




    def get_time_temp_from_slurm(self, file_to_path):
        n_columns = 2
        with open(file_to_path, "r") as file:
            lines = file.readlines()

        for i, line in enumerate(lines):
            if "Per MPI rank memory allocation (min/avg/max)" in line:
                start_index = i+2
                break
        collected_rows = []
        for line in lines[start_index:]:
            # if "error: *** JOB" in line:
            #     #print(line)
            #     break
            if any(c.isalpha() for c in line):
                break
            row = line.strip().split()
            collected_rows.append(row[:n_columns])
        return np.array(collected_rows, dtype = float)


    def get_crystallisation(self, path_to_file = None):
        if path_to_file == None:
            temp_str = "T" + str(self.temp).replace(".", "")
            cooling_rate_str = "e-%i" %(np.abs(self.cooling_rate))
            #path_to_file = "../data_online/PVA-%i/%s/%s/cryst_equil_%s_%s.txt" %(self.polymer_length, self.type_simulation, temp_str, temp_str, cooling_rate_str)
            # path_to_file = "%s/cryst_equil_%s_%s.txt" %(self.path_to_home_folder, temp_str, cooling_rate_str)
            # print(path_to_file)
            path_to_file = "%s/frac_cryst.txt" %(self.path_to_home_folder)
        try:
            cryst = np.loadtxt(path_to_file)
            return cryst
        except:
            print("No crystallisation found!")
            return 0;

    def get_mean_cluster_length(self, path_to_file = None):
        if path_to_file == None:
            temp_str = "T" + str(self.temp).replace(".", "")
            cooling_rate_str = "e-%i" %(np.abs(self.cooling_rate))
            #path_to_file = "../data_online/PVA-%i/%s/%s/mean_cluster_length_%s_%s.txt" %(self.polymer_length, self.type_simulation, temp_str, temp_str, cooling_rate_str)
            #path_to_file = "%s/mean_cluster_length_%s_%s.txt" %(self.path_to_home_folder, temp_str, cooling_rate_str)
            path_to_file = "%s/domain_analysis.txt" %(self.path_to_home_folder)
        # try:
        #     cryst = pd.read_csv(path_to_file, sep = " ", index_col = 0)
        #     return cryst
        # except:
        #     print("No mean cluster length found!")
        #     return 0;

        if os.path.exists(path_to_file):
            cryst = pd.read_csv(path_to_file, sep=" ")
            
            crystallinity = cryst["crystallinity"]
            clusters_geq_2 = cryst["clusters w/ >= 2 members"]
            no_clusters = cryst["   independent clusters"]  # keep exact spaces if present
            mean_cryst_domain_size = cryst["mean size cryst domains"]
            no_cryst_grid_elements = cryst["crystalline grid elements"]

            return (crystallinity,
                    clusters_geq_2,
                    no_clusters,
                    mean_cryst_domain_size,
                    no_cryst_grid_elements)
        else:
            print("No mean cluster length found!")
            return 0, 0, 0, 0, 0;



    def merge_runs(self, path_slurm_file, lammps_file_name_without_time, part_of_run):
        """Wraps different LAMMPS runs into one file"""
        time_temp_array =  self.get_time_temp_from_slurm(path_slurm_file) 
        lammps_data_files = self.generate_list_files(lammps_file_name_without_time, time_temp_array)
        if part_of_run > 1:
            for i in range(2, part_of_run+1):
                time_temp_2 = self.get_time_temp_from_slurm(self.insert_run_number_in_string(path_slurm_file, i))
                #Add the lammps data files to the lammps data file list
                lammps_file_name_next_run = self.insert_run_number_in_string(lammps_file_name_without_time, i)
                lammps_data_files_next_run = self.generate_list_files(lammps_file_name_next_run, time_temp_2)

                lammps_data_files_next_run = lammps_data_files_next_run[1:]
                lammps_data_files = lammps_data_files + lammps_data_files_next_run
                #Add largest timestep to new time_temp array and merge them 
                time_temp_2 = time_temp_2[1:, :]
                maxtime = np.max(time_temp_array[:, 0])
                time_temp_2[:,0] = time_temp_2[:,0] + maxtime
                #Merge time/temperature arrays 
                time_temp_array = np.vstack((time_temp_array, time_temp_2))
        return time_temp_array, lammps_data_files



                


    def insert_run_number_in_string(self, string, run_number):
        """Inserts current run number in string as _run[n] with [n] the run number"""
        run_time = f"_run{run_number}"
        run_out = f"-run{run_number}"
        
        # Check for _time pattern first (exclusive with .out)
        time_index = string.find('_time')
        out_index = string.find('.out')
        
        # Insert before _time if found
        if time_index != -1:
            new_string = string[:time_index] + run_time + string[time_index:]
            return new_string
        
        # Insert before .out if found (and no _time)
        elif out_index != -1:
            new_string = string[:out_index] + run_out + string[out_index:]
            return new_string
        
        # Neither pattern found
        else:
            return string


    def generate_list_files(self, template_lammps_file_name, time_temp_array):
        list_lammps_files = []
        for i in range(0, time_temp_array.shape[0]):
            list_lammps_files.append("%s_%s.txt" %(template_lammps_file_name, int(time_temp_array[i, 0])))
        return list_lammps_files



    def check_polymer_attributes_file_exists(self, path_to_file, calc_all = False):
        if calc_all == True:
            try:
                os.remove(path_to_file)
            except OSError:
                pass
            return self.time_temp_array, self.list_lammps_files
        if os.path.exists(path_to_file) and os.path.getsize(path_to_file) > 0:
            try:
                #print(path_to_file)
                data_array = np.loadtxt(path_to_file)
                #print(data_array)
                maxtime = int(np.max(data_array[:, 0]))
            except ValueError: 
                data_array = pd.read_csv(path_to_file, sep = " ").drop(columns = "Unnamed: 0")
                #print(data_array)
                maxtime = (data_array['time'].max())
            #Only continue at given time
            shortened_time_temp_array = self.time_temp_array[self.time_temp_array[:, 0] > maxtime]
            #Also delete first [n] entries of the list_lammps_files
            current_list_lammps_files = self.list_lammps_files[(self.time_temp_array.shape[0] - shortened_time_temp_array.shape[0]):]
            return shortened_time_temp_array, current_list_lammps_files
        else:
            with open(path_to_file, "w") as file:
                file.write("")
            shortened_time_temp_array = self.time_temp_array
            current_list_lammps_files = self.list_lammps_files
        return shortened_time_temp_array, current_list_lammps_files;


    def calc_crystallisation(self, cryst_array_string: str = None, frac_cryst_save_loc: str = None, calc_all = False, cryst_cutoff: float = 0.8):
        """Calculates the crystallisation of all timesteps of a run and saves that to a file. 
        cryst_file: .txt containing the timestep and fraction of crystallisation
        polymer_cryst_files: .txt files containing crystallisation fraction and nematic vectors per small box per simulation"""

        if cryst_array_string == None: 
            cryst_array_string = "%s/%s" %(self.path_to_home_folder, "frac_cryst.txt")
        if frac_cryst_save_loc == None: 
            path_to_nem_vectors_dir = "%s/%s" %(self.path_to_home_folder, "nematic_vectors")
            Path(path_to_nem_vectors_dir).mkdir(parents = True, exist_ok = True) #Make home directory if not exists
            frac_cryst_save_loc = "%s/%s" %(path_to_nem_vectors_dir, "nem_vectors")
        local_time_temp_array, local_list_lammps_files = self.check_polymer_attributes_file_exists(cryst_array_string, calc_all = calc_all)
        for i in tqdm(range(local_time_temp_array.shape[0])):
            current_time = local_time_temp_array[i, 0]
            current_file = local_list_lammps_files[i]

            if not os.path.exists(current_file):
                print(f"File not found, skipping: {current_file}")
                continue
            current_polymer = polymer(current_file, polymer_length= self.polymer_length)
            frac_cryst = current_polymer.atom_coords.get_nematic_vector_5(
                save_string = "%s_time_%i.txt" % (frac_cryst_save_loc, current_time),
                cryst_cutoff = cryst_cutoff
            )
            with open(cryst_array_string, "a") as file:
                file.write(f"{current_time} {frac_cryst}\n")

    def calc_avg_domain_size(self, boxes_eigv_file: str = None, load_polymer_cryst_files: str = None, calc_all: bool = False, print_results: bool = False, ndot_cutoff = 0.97):
        if boxes_eigv_file == None:
            boxes_eigv_file = "%s/domain_analysis.txt" %(self.path_to_home_folder)
        local_time_temp_array, local_list_lammps_files = self.check_polymer_attributes_file_exists(boxes_eigv_file, calc_all = calc_all)
        print(local_time_temp_array)
        if local_time_temp_array.shape[0] != 0:
            cryst_domain_array = pd.DataFrame(np.zeros([local_time_temp_array.shape[0], 7]), columns = ["time", "crystallinity", "clusters w/ >= 2 members", 
            "   independent clusters", "mean size cryst domains", "crystalline grid elements", "total volume"])
            if local_time_temp_array.shape[0] < self.time_temp_array.shape[0]:
                old_cryst_array = pd.read_csv("%s" %(boxes_eigv_file), sep = " ", index_col = 0)
        else:
            self.cryst_domain_array = pd.read_csv("%s" %(boxes_eigv_file), sep = " ", index_col = 0) #Read in domain file if it exists
            return 0;
        for i in tqdm(range(0, (local_time_temp_array.shape[0]))):
            current_time = local_time_temp_array[i, 0]
            current_file = local_list_lammps_files[i]
            if not os.path.exists(current_file):
                print(f"File not found, skipping: {current_file}")
                continue
            current_polymer = polymer(current_file, polymer_length= self.polymer_length)
            if load_polymer_cryst_files == None:
                current_polymer.read_cryst("%s/nematic_vectors/nem_vectors_time_%i.txt" 
                    %(self.path_to_home_folder, current_time))
            else:
                current_polymer.read_cryst(load_polymer_cryst_files)
            label_matrix = current_polymer.merge_boxes_2(print_results = print_results, ndot_cutoff = ndot_cutoff)
            np.save("%s/nematic_vectors/label_map_time_%i.npy" %(self.path_to_home_folder, current_time), label_matrix)
            # with open(boxes_eigv_file, "a") as file:
            #     file.write(f"{current_time} {current_polymer.results.mean_cluster_size}\n")
            test = np.array([current_time, current_polymer.results.fraction_crystallinity,
                current_polymer.results.total_number_clusters, current_polymer.results.total_number_independent_clusters, 
                current_polymer.results.mean_cluster_size, current_polymer.results.total_number_crystalline_grid_elements,
                current_polymer.atom_coords.volume])
            cryst_domain_array.iloc[i, :] = test
        if local_time_temp_array.shape[0] != 0:
            if local_time_temp_array.shape[0] < self.time_temp_array.shape[0]:
                cryst_domain_array = pd.concat([old_cryst_array, cryst_domain_array], ignore_index= True)
            cryst_domain_array.to_csv("%s" %(boxes_eigv_file), sep = " ", mode = "a", header = False)




    def calc_avg_local_density(self):
        local_time_temp_array, local_list_lammps_files = self.time_temp_array, self.list_lammps_files
        avg_local_density = local_time_temp_array
        for i in tqdm(range(0, (local_time_temp_array.shape[0]))):
        #for i in range(0, 5):
            current_time = local_time_temp_array[i, 0]
            current_file = local_list_lammps_files[i]
            current_polymer = polymer(current_file, polymer_length= self.polymer_length)
            avg_local_density[i, 1] = np.mean(current_polymer.get_density_dist())
        self.avg_local_density = avg_local_density
        return avg_local_density


    def get_volume_per_monomer(self):
        """Returns a nx2 array with in the first column the time and the second volume/amount of monomers"""
        local_time_temp_array, local_list_lammps_files = self.time_temp_array, self.list_lammps_files
        volume_monomer_array = local_time_temp_array
        for i in tqdm(range(0, (local_time_temp_array.shape[0]))):
            current_time = local_time_temp_array[i, 0]
            current_polymer = polymer(local_list_lammps_files[i], polymer_length= self.polymer_length)
            volume, n_monomers = current_polymer.atom_coords.volume, current_polymer.atom_coords.n_atoms
            volume_monomer_array[i, 1] = volume/n_monomers

        return volume_monomer_array
            

        
def plot_crystallisation_different_polymer_lengths(simulation_list, save: bool = False, savestring = None, labels_list: list = None, plot_equal_length: bool = False):

    length_list =[]

    for simulation in simulation_list:
        length_list.append(simulation.cryst.shape[0])
        if plot_equal_length == True:
            length = min(length_list)
        else:
            length = max(length_list)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # First plot mean cluster length 
    for simulation in simulation_list:
        # ax1.plot(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length.iloc[:length, 3], 
        #     label = r"independent clusters, PVA-%i" %(simulation.polymer_length))
        label_cryst = r"$\phi$,  PVA-%i" %(simulation.polymer_length)
        print(simulation.cryst)
        ax1.plot(simulation.time_temp_array[:length, 0], simulation.cryst[:length, 1], 
            label = label_cryst)
        ax2.scatter(simulation.time_temp_array[:length, 0], (simulation.mean_cryst_domain_size[:length])**(1/3) *2,
            label =  r"$\langle \sigma \rangle$, PVA-%i" %(simulation.polymer_length))
        temp = simulation.temp

    ax1.set_xlabel("time")
    #ax1.set_ylabel(r"amount of independent clusters")
    ax1.set_ylabel(r"$\phi$")
    ax2.set_ylabel(r"length scale ($\sigma$)")
    #fig.tight_layout()
    fig.legend(loc = "lower right", bbox_to_anchor=(0.895, 0.115))
    #fig.suptitle(r"Independent clusters and mean domain size, $\dot{T} = 10^{%i}$" %(simulation.cooling_rate))
    fig.suptitle(r"Isothermal crystallisation, T=%.2f"%(temp))
    if savestring != None:
        fig.savefig("plots/%s" %savestring)
    #fig.savefig("plots/e%i_no_clusters_mean_domain_size.pdf"%(simulation.cooling_rate))
    plt.show()


def plot_mean_domain_size_indep_clusters(simulation_list, save: bool = False, savestring = None, labels_list: list = None, plot_equal_length: bool = False):

    length_list =[]

    for simulation in simulation_list:
        length_list.append(simulation.cryst.shape[0])
        if plot_equal_length == True:
            length = min(length_list)
        else:
            length = max(length_list)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # First plot mean cluster length 
    for simulation in simulation_list:
        print(simulation.polymer_length, simulation.n_atoms)
        # ax1.plot(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length.iloc[:length, 3], 
        #     label = r"independent clusters, PVA-%i" %(simulation.polymer_length))
        local_time_temp_array, local_list_lammps_files = simulation.time_temp_array, simulation.list_lammps_files
        volume_monomer_array = local_time_temp_array
        for i in tqdm(range(0, (local_time_temp_array.shape[0]))):
            current_polymer = polymer(local_list_lammps_files[i], polymer_length= simulation.polymer_length)
            volume_monomer_array[i, 1] = current_polymer.atom_coords.volume
        label_cryst = r"number of clusters,  PVA-%i" %(simulation.polymer_length)
        ax1.plot(simulation.time_temp_array[:length, 0], simulation.clusters_geq_2[:length]/volume_monomer_array[:length, 1], 
            label = label_cryst)
        ax2.scatter(simulation.time_temp_array[:length, 0], (simulation.mean_cryst_domain_size[:length])**(1/3) *2,
            label =  r"$\langle \sigma \rangle$, PVA-%i" %(simulation.polymer_length))
        temp = simulation.temp

    ax1.set_xlabel("time")
    #ax1.set_ylabel(r"amount of independent clusters")
    ax1.set_ylabel(r"clusters/volume")
    ax2.set_ylabel(r"length scale ($\sigma$)")
    #fig.tight_layout()
    fig.legend(loc = "lower right", bbox_to_anchor=(0.895, 0.115))
    #fig.suptitle(r"Independent clusters and mean domain size, $\dot{T} = 10^{%i}$" %(simulation.cooling_rate))
    fig.suptitle(r"Isothermal crystallisation, T=%.2f"%(temp))
    if savestring != None:
        fig.savefig("plots/%s" %savestring)
    #fig.savefig("plots/e%i_no_clusters_mean_domain_size.pdf"%(simulation.cooling_rate))
    plt.show()


def plot_no_clusters(simulation_list, save: bool = False, savestring = None, labels_list: list = None, plot_equal_length: bool = False):

    length_list =[]

    for simulation in simulation_list:
        length_list.append(simulation.cryst.shape[0])
        if plot_equal_length == True:
            length = min(length_list)
        else:
            length = max(length_list)

    for simulation in simulation_list:

        label_cryst = r"PVA-%i" %(simulation.polymer_length)
        local_time_temp_array, local_list_lammps_files = simulation.time_temp_array, simulation.list_lammps_files
        volume_monomer_array = local_time_temp_array
        for i in tqdm(range(0, (local_time_temp_array.shape[0]))):
            current_polymer = polymer(local_list_lammps_files[i], polymer_length= simulation.polymer_length)
            volume_monomer_array[i, 1] = current_polymer.atom_coords.volume
        label_cryst = r"number of clusters,  PVA-%i" %(simulation.polymer_length)
        plt.plot(simulation.time_temp_array[:length, 0], simulation.clusters_geq_2[:length]/volume_monomer_array[:length, 1], 
            label = label_cryst)
    plt.xlabel("time")
    #ax1.set_ylabel(r"amount of independent clusters")
    plt.ylabel(r"number of clusters/ total volume")
    plt.title("Cluster size, only polymer lengths")
    plt.legend()
    if savestring != None:
        plt.savefig("plots/%s" %savestring)
    #fig.savefig("plots/e%i_no_clusters_mean_domain_size.pdf"%(simulation.cooling_rate))
    plt.show()
    
def main():
    data_path = "../../data/pva-100/quick_quench/long_run"

    #cooling_e3_T07 = Simulation(0.7, -3, "%s%s" %(data_path, "/slurm-e3-T=07.out"), "%s/equil_t_07_tdot_e-3_time"%(data_path), no_runs = 3)
    #cooling_e3_T07.calc_crystallisation("../data_online/PVA-100/quick_quench/T07/cryst_equil_T07_e-3.txt", "../data_online/PVA-100/quick_quench/T07/boxes_eigenvalues/equil_t_085_tdot_e-3")
    cooling_e3_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T085.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path), no_runs = 3, 
        home_folder= "../data_online/PVA-100/quick_quench/T085")
    #cooling_e3_T085.calc_crystallisation("../data_online/PVA-100/quick_quench/T085/cryst_equil_T085_e-3.txt", "../data_online/PVA-100/quick_quench/T085/boxes_eigenvalues/equil_t_085_tdot_e-3")
    #cooling_e3_T085.calc_avg_domain_size("../data_online/PVA-100/quick_quench/T085/mean_cluster_length_T085_e-3.txt")

if __name__== "__main__":
    main()

