from analyse7 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from pathlib import Path
import pandas as pd

cryst_path_prefix = "../crystallisation"






class Simulation:
    """Class for one run, holds multiple polymer objects."""

    def __init__(self, temperature: float, cooling_rate: int, path_slurm_file: str, lammps_file_name_without_time :str, run_id: int = 1,
        no_runs: int = 1, polymer_length: int = 100, home_folder: str = None):
        """lammps_file_name_without_time should not contain a time"""
        self.temp = temperature
        self.cooling_rate = cooling_rate

        self.template_lammps_file_name = lammps_file_name_without_time
        self.time_temp_array, self.list_lammps_files = self.merge_runs(path_slurm_file, lammps_file_name_without_time, no_runs)
        self.max_temp, self.min_temp = np.max(self.time_temp_array[:, 0]), np.min(self.time_temp_array[:, 0])
        self.no_polymers = self.time_temp_array.shape[0]
        self.polymer_length = polymer_length
        self.path_to_home_folder = home_folder
        self.cryst = self.get_crystallisation() #2D array containing time and crystallisation at time
        self.mean_cluster_length = self.get_mean_cluster_length()
        Path(home_folder).mkdir(parents = True, exist_ok = True) #Make home directory if not exists

    def get_polymer_by_count(self, count):
        """Gets i-th polymer"""
        if isinstance(count, int):
            return polymer(self.list_lammps_files[count])
        elif isinstance(count, np.ndarray):
            polymer_list = []
            for i in count:
                polymer_list.append(polymer(self.list_lammps_files[i]))
            return polymer_list
        else:
            raise Exception("count must be an int or a numpy.ndarray")


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
            path_to_file = "%s/cryst_equil_%s_%s.txt" %(self.path_to_home_folder, temp_str, cooling_rate_str)
            print(path_to_file)
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
            path_to_file = "%s/mean_cluster_length_%s_%s.txt" %(self.path_to_home_folder, temp_str, cooling_rate_str)
            print(path_to_file)
        try:
            cryst = np.loadtxt(path_to_file)
            return cryst
        except:
            print("No mean cluster length found!")
            return 0;



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
            # Check what latest timestep written in file
            data_array = np.loadtxt(path_to_file)
            maxtime = int(np.max(data_array[:, 0]))
            #Only continue at given time
            shortened_time_temp_array = self.time_temp_array[self.time_temp_array[:, 0] > maxtime]
            #Also delete first [n] entries of the list_lammps_files
            current_list_lammps_files = self.list_lammps_files[(self.time_temp_array.shape[0] - shortened_time_temp_array.shape[0]):]

        else:
            with open(path_to_file, "w") as file:
                file.write("")
            shortened_time_temp_array = self.time_temp_array
            current_list_lammps_files = self.list_lammps_files
        return shortened_time_temp_array, current_list_lammps_files;


    def calc_crystallisation(self, cryst_array_string: str, frac_cryst_save_loc: str, calc_all = False):
        """Calculates the crystallisation of all timesteps of a run and saves that to a file. 
        cryst_file: .txt containing the timestep and fraction of crystallisation
        polymer_cryst_files: .txt files containing crystallisation fraction and nematic vectors per small box per simulation"""
        # Check if dirs exists and if not, create them 
        cryst_path = Path(cryst_array_string).parent
        cryst_path.mkdir(parents=True, exist_ok = True)
        boxes_eigenvalues_path = Path(frac_cryst_save_loc).parent
        boxes_eigenvalues_path.mkdir(parents=True, exist_ok = True)


        local_time_temp_array, local_list_lammps_files = self.check_polymer_attributes_file_exists(cryst_array_string, calc_all = calc_all)
        for i in tqdm(range(0, (local_time_temp_array.shape[0]))):
            current_time = local_time_temp_array[i, 0]
            current_file = local_list_lammps_files[i]
            current_polymer = polymer(current_file)
            frac_cryst = current_polymer.atom_coords.get_nematic_vector_5(save_string = "%s_time_%i.txt" %(frac_cryst_save_loc, current_time))
            with open(cryst_array_string, "a") as file:
                file.write(f"{current_time} {frac_cryst}\n")

    def calc_avg_domain_size(self, boxes_eigv_file: str, load_polymer_cryst_files: str = None, calc_all: bool = False):
        #TODO: implement more sophisticated data structure 
        local_time_temp_array, local_list_lammps_files = self.check_polymer_attributes_file_exists(boxes_eigv_file, calc_all = calc_all)
        cryst_domain_array = pd.DataFrame(np.zeros([local_time_temp_array.shape[0], 6]), columns = ["time", "crystallinity", "clusters w/ >= 2 members", 
            "independent clusters", "mean size cryst domains", "crystalline grid elements"])
        for i in tqdm(range(0, (local_time_temp_array.shape[0]))):
            current_time = local_time_temp_array[i, 0]
            current_file = local_list_lammps_files[i]
            current_polymer = polymer(current_file)
            if load_polymer_cryst_files == None:
                current_polymer.read_cryst("%s/boxes_eigenvalues_e%s/equil_t_%s_tdot_e%s_time_%i.txt" 
                    %(self.path_to_home_folder, self.cooling_rate, str(self.temp).replace(".", ""),self.cooling_rate, current_time))
            else:
                current_polymer.read_cryst(load_polymer_cryst_files)
            current_polymer.merge_boxes()

            # with open(boxes_eigv_file, "a") as file:
            #     file.write(f"{current_time} {current_polymer.results.mean_cluster_size}\n")
            test = np.array([current_time, current_polymer.results.fraction_crystallinity,
                current_polymer.results.total_number_clusters, current_polymer.results.total_number_independent_clusters, 
                current_polymer.results.mean_cluster_size, current_polymer.results.total_number_crystalline_grid_elements])
            cryst_domain_array.iloc[i, :] = test
        cryst_domain_array.to_csv("%s" %(boxes_eigv_file), sep = " ", mode = "a")




    def calc_avg_local_density(self):
        local_time_temp_array, local_list_lammps_files = self.time_temp_array, self.list_lammps_files
        avg_local_density = local_time_temp_array
        for i in tqdm(range(0, (local_time_temp_array.shape[0]))):
        #for i in range(0, 5):
            current_time = local_time_temp_array[i, 0]
            current_file = local_list_lammps_files[i]
            current_polymer = polymer(current_file)
            avg_local_density[i, 1] = np.mean(current_polymer.get_density_dist())
        self.avg_local_density = avg_local_density
        return avg_local_density



        


    
def main():
    data_path = "../../data/pva-100/quick_quench/long_run"

    #cooling_e3_T07 = Simulation(0.7, -3, "%s%s" %(data_path, "/slurm-e3-T=07.out"), "%s/equil_t_07_tdot_e-3_time"%(data_path), no_runs = 3)
    #cooling_e3_T07.calc_crystallisation("../data_online/PVA-100/quick_quench/T07/cryst_equil_T07_e-3.txt", "../data_online/PVA-100/quick_quench/T07/boxes_eigenvalues/equil_t_085_tdot_e-3")
    cooling_e3_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T085.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path), no_runs = 3, 
        home_folder= "../data_online/PVA-100/quick_quench/T085")
    #cooling_e3_T085.calc_crystallisation("../data_online/PVA-100/quick_quench/T085/cryst_equil_T085_e-3.txt", "../data_online/PVA-100/quick_quench/T085/boxes_eigenvalues/equil_t_085_tdot_e-3")
    cooling_e3_T085.calc_avg_domain_size("../data_online/PVA-100/quick_quench/T085/mean_cluster_length_T085_e-3.txt")

if __name__== "__main__":
    main()

