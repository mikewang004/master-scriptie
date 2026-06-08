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
            dataframe_slurm = pd.read_csv("%s/%s_sim_data.txt" %(self.path_to_home_folder, self.slurm_prefix), sep = " ", index_col = "num_in_run")
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





class Simulation: 
    """Class to read in and manipulate sequences of lammps .txt dump files and corresponding slurm files """

    def __init__(self, polymer_length: int, path_to_data_folder: str, path_to_home_folder: str, cooling_rate: int = -3,
            target_temp: float = 0.88, cryst_cutoff: float = 0.8, ndot_cutoff: float = 0.97):
        self.slurm_files = SlurmFiles(path_to_data_folder, path_to_home_folder)
        self.target_temp = target_temp
        self.path_to_data_folder = path_to_data_folder
        self.path_to_home_folder = path_to_home_folder
        self.df_slurm_sim_data = pd.read_csv("%s/%s_sim_data.txt" %(self.path_to_home_folder, self.slurm_files.slurm_prefix))
        self.list_run_files = self.get_list_run_files()





        

    def get_list_run_files(self):
        """Reads all slurm files, makes a new slurm file containing the step/temp/E_pair/E_mol/TotEng/Press/Vol of all slurm files run"""

        folder = pathlib.Path(self.path_to_data_folder)
        pattern = re.compile(r'_run(\d+)_time_(\d+)\.txt$')

        # -----------------------------------------------------------------
        # Build a list of (j, i, filename) tuples for every matching file
        # -----------------------------------------------------------------
        entries: List[tuple[int, int, str]] = []
        for p in folder.glob('*_run*_time_*.txt'):
            m = pattern.search(p.name)
            if m:
                i = int(m.group(1))   # run index
                j = int(m.group(2))   # time index
                entries.append((j, i, p.name))

        # -----------------------------------------------------------------
        # Sort by the tuple (j, i) – time first, then run
        # -----------------------------------------------------------------
        entries.sort(key=lambda x: (x[1], x[0]))

        # -----------------------------------------------------------------
        # Return only the filenames in the desired order
        # -----------------------------------------------------------------
        sorted_filenames = [name for _, _, name in entries]

        # (optional) print for debugging
        return sorted_filenames




def main():
    #PVA_200 = Simulation(100, "../../data/PVA-200/equil", "../data_online/PVA-200/icryst_T088_Tdot_e-3")
    PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")

if __name__== "__main__":
    main()