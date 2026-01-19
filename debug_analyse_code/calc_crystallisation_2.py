from analyse7 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


cryst_path_prefix = "../crystallisation"




class Simulation:
    """Class for one run, holds multiple polymer objects."""

    def __init__(self, temperature: float, cooling_rate: float, path_slurm_file: str, lammps_file_name_without_time :str):
        """lammps_file_name_without_time should not contain a time"""
        self.temp = temperature
        self.cooling_rate = cooling_rate
        self.time_temp_array = self.get_time_temp_from_slurm(path_slurm_file) #0th column time 1st column temperature

        self.max_temp, self.min_temp = np.max(self.time_temp_array[:, 0]), np.min(self.time_temp_array[:, 0])

        self.template_lammps_file_name = lammps_file_name_without_time
        self.generate_list_files()


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


    def generate_list_files(self):
        self.list_lammps_files = []
        for i in range(0, self.time_temp_array.shape[0]):
            self.list_lammps_files.append("%s_%s.txt" %(self.template_lammps_file_name, int(self.time_temp_array[i, 0])))
        print(self.list_lammps_files)


    def calc_crystallisation(self, cryst_file: str, polymer_cryst_files: str):
        """Calculates the crystallisation of all timesteps of a run and saves that to a file. 
        cryst_file: .txt containing the timestep and fraction of crystallisation
        polymer_cryst_files: .txt files containing crystallisation fraction and nematic vectors per small box per simulation"""
        cryst_array_string = cryst_file
        with open(cryst_array_string, "w") as file:
            file.write("")
        for i in range(0, (self.time_temp_array.shape[0])):
            current_time = self.time_temp_array[i, 0]
            current_file = self.list_lammps_files[i]
            print(current_file)
            current_polymer = polymer(current_file)
            frac_cryst = current_polymer.atom_coords.get_nematic_vector_5(save_string = "%s_time_%i.txt" %(polymer_cryst_files, current_time))
            # print(current_file_name)
            with open(cryst_array_string, "a") as file:
                file.write(f"{current_time} {frac_cryst}\n")

    def calc_avg_domain_size(self, boxes_eigv_file: str, polymer_box_size_file: str):
        with open(polymer_box_size_file, "w") as file:
            file.write("") 

        


    

data_path = "../../data/pva-100/quick_quench/long_run"

# cooling_e3_T07 = Simulation(0.7, -3, "%s%s" %(data_path, "/slurm-T=0.7.out"), "%s/equil_t_07_tdot_e-3_time"%(data_path))
# cooling_e3_T07.calc_crystallisation("../crystallisation/PVA-100/quick_quench/T07/cryst_equil_T07_e-3.txt", "%s/equil_t_07_tdot_e-3"%(data_path))

cooling_e3_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-T=0.85.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path))
cooling_e3_T085.calc_crystallisation("../crystallisation/PVA-100/quick_quench/T08/cryst_equil_T08_e-3.txt", "../crystallisation/PVA-100/quick_quench/T08/equil_t_085_tdot_e-3")