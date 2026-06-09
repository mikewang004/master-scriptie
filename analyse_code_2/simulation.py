from analyse8 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import pathlib
import re
from typing import List, Tuple, Optional

#TODO: modify file so that polymer properties can be accessed easily from Simulation()

def make_folder(self, path_to_folder: str):
    """Check if folder exists, if not make folder"""
    if not os.path.exists(path_to_folder):
        os.makedirs(path_to_folder)
    return 0;


class SlurmFiles():
    """Class to handle Slurm files and merge them to a new data file containing step/temp/E_pair/E_mol/TotEng/Press/Vol
    Steps of run%i files are summed upon the last entry of run%i-1 files, so that the step forms one contineous sequence."""

    def __init__(self, path_to_data_folder, path_to_home_folder):
        self.path_to_data_folder = path_to_data_folder; self.path_to_home_folder = path_to_home_folder
        self.slurm_prefix, self.list_slurm_files = self.get_list_slurm_files()
        #self.read_slurm_file("%s/%s" %(self.path_to_home_folder, self.list_slurm_files[0]))
        self.merge_slurm_files()


    def get_list_slurm_files(self):
        """Returns a list of all slurm files present in a folder."""

        print(self.path_to_data_folder)
        folder = pathlib.Path(self.path_to_data_folder)
        pattern = re.compile(r'-run(\d+)\.out$')
        files = []
        for p in folder.glob('*-run*.out'):
            m = pattern.search(p.name)
            if m:
                files.append((int(m.group(1)), p.name))
        files.sort(key=lambda x: x[0])


        pattern = re.compile(r'^(.*)-run(\d+)\.out$')
        for q in folder.glob('*-run*.out'):
            m = pattern.search(q.name)
            if m:
                prefix = m.group(1)
                break  # first match is enough
        list_slurm_files = [name for _, name in files]
        return prefix, list_slurm_files


    def read_slurm_file(self, file_to_path):
        """Reads the columns containing step/temp/E_pair/E_mol/TotEng/Press/Vol of a singular Slurm Files"""
        n_columns = 7
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
        slurm_data = pd.DataFrame(collected_rows, columns = ["Step", "Temp", "E_pair", "E_mol", "TotEng", "Press", "Volume"])
        slurm_data = slurm_data.astype({
            "Step": "int64",
            "Temp": "float64",
            "E_pair": "float64",
            "E_mol": "float64",
            "TotEng": "float64",
            "Press": "float64",
            "Volume": "float64",
        })
        return slurm_data

    def read_slurm_summary_file(self):
        try:
            dataframe_slurm = pd.read_csv("%s/%s_sim_data.txt" %(self.path_to_home_folder, self.slurm_prefix), sep = " ")#, index_col = "num_in_run")
        except FileNotFoundError:
            dataframe_slurm is None

        #print(dataframe_slurm)
        return dataframe_slurm

    def merge_slurm_files(self):
        """Reads all slurm files, makes a new slurm file containing the step/temp/E_pair/E_mol/TotEng/Press/Vol of all slurm files run"""
        print(self.list_slurm_files)

        dataframe_slurm_existing = self.read_slurm_summary_file()

        dataframe_slurm = self.read_slurm_file("%s/%s" %(self.path_to_data_folder, self.list_slurm_files[0]))


        dataframe_slurm.insert(0, "Run", 0)
        dataframe_slurm.insert(1, "NumInRun", dataframe_slurm.index)
        step_pos = dataframe_slurm.columns.get_loc("Step")
        dataframe_slurm.insert(step_pos + 1, "StepSequence", dataframe_slurm["Step"].copy())
        #print(dataframe_slurm)
        last_time = dataframe_slurm.iloc[-1, 1]
        for i in range(1, len(self.list_slurm_files)):
            slurm_data = self.read_slurm_file("%s/%s" %(self.path_to_data_folder, self.list_slurm_files[i]))
            slurm_data.insert(0, "Run", i)
            slurm_data.insert(1, "NumInRun", slurm_data.index)
            #slurm_data.index()
            step_pos = slurm_data.columns.get_loc("Step")
            slurm_data.insert(step_pos + 1, "StepSequence", slurm_data["Step"].copy())
           # print(slurm_data)
            slurm_data = slurm_data.drop([0, 0])
            slurm_data.iloc[:, 2] += last_time
            #print(slurm_data)

            dataframe_slurm = pd.concat([dataframe_slurm, slurm_data])
            last_time = slurm_data.iloc[-1, 2]


        dataframe_slurm.reset_index(drop = True, inplace = True, names="Index")
        no_files_current_df_slurm = dataframe_slurm.shape[0]
        no_files_saved_df_slurm = dataframe_slurm_existing.shape[0]
        if no_files_current_df_slurm != no_files_saved_df_slurm:
            dataframe_slurm.to_csv("%s/%s_sim_data.txt" %(self.path_to_home_folder, self.slurm_prefix), sep = " ", index_label = "num_in_run")
        return 0;


class domain_analysis:
    """Class to calculate and save crystallisation files"""

    def __init__(self, path_to_data_folder: str, path_to_home_folder: str, list_run_files: list, slurm_file, lammps_dump_prefix: str,
        cryst_cutoff: float = 0.8, ndot_cutoff: float = 0.97):
        self.path_to_data_folder = path_to_data_folder
        self.path_to_home_folder = path_to_home_folder
        self.list_run_files = list_run_files
        self.df_slurm_sim_data = slurm_file
        self.lammps_dump_prefix = lammps_dump_prefix
        self.cryst_cutoff = cryst_cutoff; self.ndot_cutoff = ndot_cutoff




    def calc_crystallisation(self, cryst_array_string = None):
        path_to_crystallisation_folder = "%s/crystallisation" %(self.path_to_home_folder)
        make_folder(path_to_crystallisation_folder)
        if cryst_array_string == None: 
            cryst_array_string = "%s/%s" %(self.path_to_home_folder, "frac_cryst.txt")
        for i in tqdm(range(0, len(self.list_run_files))):
            #Read file 
            current_polymer = polymer("%s/%s" %(self.path_to_data_folder, self.list_run_files[i]))
            current_time = self.df_slurm_sim_data["StepSequence"].iloc[i]
            frac_cryst = current_polymer.atom_coords.get_nematic_vector_5(
                save_string = "%s/%s_cryst_time_%s.txt" %(path_to_crystallisation_folder, self.lammps_dump_prefix, current_time),
                cryst_cutoff = self.cryst_cutoff
            )

            with open(cryst_array_string, "a") as file:
                file.write(f"{current_time} {frac_cryst}\n")
        return 0;

    def calc_avg_domain_size(self):
        # Read in crystallisation files
        path_to_nematic_vectors_folder = "%s/nematic_vectors" %(self.path_to_home_folder)
        cryst_domain_array = pd.DataFrame(np.zeros([self.df_slurm_sim_data.shape[0], 7]), columns = ["time", "crystallinity", "clusters w/ >= 2 members", 
        "   independent clusters", "mean size cryst domains", "crystalline grid elements", "total volume"])
        for i in tqdm(range(0, len(self.list_run_files))):
            #Read file 
            current_polymer = polymer("%s/%s" %(self.path_to_data_folder, self.list_run_files[i]))
            current_time = self.df_slurm_sim_data["StepSequence"].iloc[i]
            try:
                current_cryst_file = pd.read_csv("%s/%s_cryst_time_%s.txt" %(path_to_crystallisation_folder, self.lammps_dump_prefix, current_time), sep = " ")
            except FileNotFoundError:
                print(f"File not found, skipping: {current_time}")
                continue
            label_matrix = current_polymer.merge_boxes_2(print_results = print_results, ndot_cutoff = ndot_cutoff)
            np.save("%s/nematic_vectors/label_map_time_%i.npy" %(self.path_to_home_folder, current_time), label_matrix)
        test = np.array([current_time, current_polymer.results.fraction_crystallinity,
            current_polymer.results.total_number_clusters, current_polymer.results.total_number_independent_clusters, 
            current_polymer.results.mean_cluster_size, current_polymer.results.total_number_crystalline_grid_elements,
            current_polymer.atom_coords.volume])
        cryst_domain_array.iloc[i, :] = test
        cryst_domain_array.to_csv("%s/domain_analysis.txt" %(self.path_to_home_folder), sep = " ", header = True)


class Simulation: 
    """Class to read in and manipulate sequences of lammps .txt dump files and corresponding slurm files """

    def __init__(self, polymer_length: int, path_to_data_folder: str, path_to_home_folder: str, cooling_rate: int = -3,
            target_temp: float = 0.88, cryst_cutoff: float = 0.8, ndot_cutoff: float = 0.97):
        self.polymer_length = polymer_length
        self.cryst_cutoff = cryst_cutoff; self.ndot_cutoff = ndot_cutoff
        self.slurm_files = SlurmFiles(path_to_data_folder, path_to_home_folder)
        self.target_temp = target_temp
        self.path_to_data_folder = path_to_data_folder
        self.path_to_home_folder = path_to_home_folder
        self.df_slurm_sim_data = pd.read_csv("%s/%s_sim_data.txt" %(self.path_to_home_folder, self.slurm_files.slurm_prefix), sep = " ")
        self.list_run_files = self.get_list_run_files()
        #print(self.list_run_files)
        self.lammps_dump_prefix = self.get_file_prefix()
        self.domain_calc = domain_analysis(path_to_data_folder, path_to_home_folder, self.list_run_files, self.df_slurm_sim_data, self.lammps_dump_prefix,
            cryst_cutoff = self.cryst_cutoff, ndot_cutoff = self.ndot_cutoff)




        

    def get_list_run_files(self):
        """Reads all slurm files, makes a new slurm file containing the step/temp/E_pair/E_mol/TotEng/Press/Vol of all slurm files run"""
        folder = pathlib.Path(self.path_to_data_folder)
        pattern = re.compile(r'_run(\d+)_time_(\d+)\.txt$')

        entries: List[Tuple[int, int, str]] = []  # (time j, run i, name)

        for p in folder.glob('*_run*_time_*.txt'):
            m = pattern.search(p.name)
            if m:
                i = int(m.group(1))   # run index
                j = int(m.group(2))   # time index

                # Skip files matching *_run%i_time_0.txt with i >= 2
                if j == 0 and i >= 2:
                    continue

                entries.append((j, i, p.name))

        # Sort by (run i, time j)
        entries.sort(key=lambda x: (x[1], x[0]))

        sorted_filenames = [name for _, _, name in entries]
        return sorted_filenames

    def get_file_prefix(self):
        folder = pathlib.Path(self.path_to_data_folder)
        pattern = re.compile(r'^(.*)_run(\d+)(.*).txt$')
        for q in folder.glob('*_run*.txt'):
            m = pattern.search(q.name)
            if m:
                prefix = m.group(1)
                break  # first match is enough

        return prefix


    def read_in_run_file(self):
        return polymer(path_to_file, polymer_length=self.polymer_length)



def main():
    PVA_200 = Simulation(100, "../../data/PVA-200/equil", "../data_online/PVA-200/icryst_T088_Tdot_e-3")
    #PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")
    PVA_200.domain_calc.calc_crystallisation()

if __name__== "__main__":
    main()