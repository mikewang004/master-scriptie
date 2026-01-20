from analyse7 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os


cryst_path_prefix = "../crystallisation"






class Simulation:
    """Class for one run, holds multiple polymer objects."""

    def __init__(self, temperature: float, cooling_rate: float, path_slurm_file: str, lammps_file_name_without_time :str, no_runs: int = 1):
        """lammps_file_name_without_time should not contain a time"""
        self.temp = temperature
        self.cooling_rate = cooling_rate




        self.template_lammps_file_name = lammps_file_name_without_time
        #self.time_temp_array = self.get_time_temp_from_slurm(path_slurm_file) #0th column time 1st column temperature
        #self.list_lammps_files = self.generate_list_files(self.template_lammps_file_name, self.time_temp_array)
        self.time_temp_array, self.list_lammps_files = self.merge_runs(path_slurm_file, lammps_file_name_without_time, no_runs)
        self.max_temp, self.min_temp = np.max(self.time_temp_array[:, 0]), np.min(self.time_temp_array[:, 0])


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
            if "error: *** JOB" in line:
                print(line)
                break
            row = line.strip().split()
            collected_rows.append(row[:n_columns])
        return np.array(collected_rows, dtype = float)

    def merge_runs(self, path_slurm_file, lammps_file_name_without_time, no_runs):
        """Wraps different LAMMPS runs into one file"""
        time_temp_array =  self.get_time_temp_from_slurm(path_slurm_file) 
        lammps_data_files = self.generate_list_files(lammps_file_name_without_time, time_temp_array)
        if no_runs > 1:
            for i in range(2, no_runs+1):
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



    def check_polymer_attributes_file_exists(self, path_to_file):
        if os.path.exists(path_to_file):
            # Check what latest timestep written in file
            data_array = np.loadtxt(path_to_file)
            maxtime = int(np.max(data_array[:, 0]))
            print(maxtime)
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


    def calc_crystallisation(self, cryst_file: str, frac_cryst_save_loc: str):
        """Calculates the crystallisation of all timesteps of a run and saves that to a file. 
        cryst_file: .txt containing the timestep and fraction of crystallisation
        polymer_cryst_files: .txt files containing crystallisation fraction and nematic vectors per small box per simulation"""
        cryst_array_string = cryst_file
        local_temp_time_array, local_list_lammps_files = self.check_polymer_attributes_file_exists(cryst_array_string)
        for i in range(0, (local_temp_time_array.shape[0])):
            current_time = local_temp_time_array[i, 0]
            current_file = local_list_lammps_files[i]
            current_polymer = polymer(current_file)
            frac_cryst = current_polymer.atom_coords.get_nematic_vector_5(save_string = "%s_time_%i.txt" %(frac_cryst_save_loc, current_time))
            with open(cryst_array_string, "a") as file:
                file.write(f"{current_time} {frac_cryst}\n")

    def calc_avg_domain_size(self, boxes_eigv_file: str, polymer_box_size_file: str, load_polymer_cryst_files: str = None):
        #TODO test and finish function
        with open(polymer_box_size_file, "w") as file:
            file.write("") 
        for i in range(0, (self.time_temp_array.shape[0])):
            current_time = self.time_temp_array[i, 0]
            current_file = self.list_lammps_files[i]
            current_polymer = polymer(current_file)
            if isinstance(load_polymer_cryst_files, str):
                polymer.read_cryst("%s_time_%i.txt" %(load_polymer_cryst_files, current_time))
            polymer.merge_boxes()

            # with open(cryst_array_string, "a") as file:
            #     file.write(f"{current_time} {frac_cryst}\n")

        

        


    

data_path = "../../data/pva-100/quick_quench/long_run"

# cooling_e3_T07 = Simulation(0.7, -3, "%s%s" %(data_path, "/slurm-T=0.7.out"), "%s/equil_t_07_tdot_e-3_time"%(data_path))
# cooling_e3_T07.calc_crystallisation("../crystallisation/PVA-100/quick_quench/T07/cryst_equil_T07_e-3.txt", "%s/equil_t_07_tdot_e-3"%(data_path))

cooling_e3_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T=0.85.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path), no_runs = 2)
cooling_e3_T085.calc_crystallisation("../data_online/PVA-100/quick_quench/T085/cryst_equil_T085_e-3.txt", "../data_online/PVA-100/quick_quench/T085/boxes_eigenvalues/equil_t_085_tdot_e-3")