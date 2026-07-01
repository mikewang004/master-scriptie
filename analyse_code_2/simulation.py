from analyse9 import polymer, atom_coords
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
import bisect

#TODO: modify file so that polymer properties can be accessed easily from Simulation()

def make_folder(path_to_folder: str):
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
            dataframe_slurm = np.array([0,0])
            return dataframe_slurm

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
        last_time = dataframe_slurm.iloc[-1, 2]
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

        dataframe_slurm["Run"] = dataframe_slurm["Run"] + 1
        dataframe_slurm.reset_index(drop = True, inplace = True, names="Index")
        no_files_current_df_slurm = dataframe_slurm.shape[0]
        no_files_saved_df_slurm = dataframe_slurm_existing.shape[0]
        if no_files_current_df_slurm != no_files_saved_df_slurm:
            make_folder(self.path_to_home_folder)
            dataframe_slurm.to_csv("%s/%s_sim_data.txt" %(self.path_to_home_folder, self.slurm_prefix), sep = " ", index_label = "index")
        return 0;


class domain_analysis:
    """Class to calculate and save crystallisation files"""

    # def __init__(self, path_to_data_folder: str, path_to_home_folder: str, list_run_files: list, slurm_file, lammps_dump_prefix: str,
    #     cryst_cutoff: float = 0.8, ndot_cutoff: float = 0.97):
    #     self.path_to_data_folder = path_to_data_folder
    #     self.path_to_home_folder = path_to_home_folder
    #     self.list_run_files = list_run_files
    #     self.df_slurm_sim_data = slurm_file
    #     self.lammps_dump_prefix = lammps_dump_prefix
    #     self.cryst_cutoff = cryst_cutoff; self.ndot_cutoff = ndot_cutoff
    #     self.path_to_crystallisation_folder = "%s/crystallisation" %(self.path_to_home_folder)
    #     self.path_to_nematic_vectors_folder = "%s/nematic_vectors" %(self.path_to_home_folder)

        #print(self.df_slurm_sim_data, self.list_run_files)

    def __init__(self, sim):
        self.sim = sim




    def insert_cryst_line_sorted(self, file_path, current_time, frac_cryst):
        """
        Insert 'current_time frac_cryst' in the file so that lines stay sorted
        by current_time (ascending). New line is placed after the closest
        smaller or equal current_time.
        """
        new_time = float(current_time)  # or int(current_time) if always integer
        new_line = f"{current_time} {frac_cryst}\n"

        lines = []
        if os.path.exists(file_path):
            with open(file_path, "r") as f:
                lines = [ln for ln in f if ln.strip()]

        # Extract existing times from first column
        times = [float(ln.split()[0]) for ln in lines]

        # Find insertion position: after closest smaller/equal value
        idx = bisect.bisect_right(times, new_time)

        # Insert new line
        lines.insert(idx, new_line)

        # Rewrite the whole file
        with open(file_path, "w") as f:
            f.writelines(lines)


    def calc_crystallisation(self, path_to_cryst_array = None, file_override = False):
        #path_to_crystallisation_folder = "%s/crystallisation" %(self.path_to_home_folder)
        make_folder(self.sim.path_to_crystallisation_folder)
        if path_to_cryst_array == None: 
            path_to_cryst_array = "%s/%s" %(self.sim.path_to_home_folder, "frac_cryst.txt")
        for i in tqdm(range(0, len(self.sim.list_run_files))):
            current_time = self.sim.df_slurm_sim_data["Step"].iloc[i]
            #Read file 
            if not Path("%s/%s_cryst_time_%s.txt" %(self.sim.path_to_crystallisation_folder, self.sim.lammps_dump_prefix, current_time)).is_file():
                current_polymer = polymer("%s/%s" %(self.sim.path_to_data_folder, self.sim.list_run_files[i]))
                frac_cryst = current_polymer.atom_coords.get_nematic_vector_5(
                    save_string = "%s/%s_cryst_time_%s.txt" %(self.sim.path_to_crystallisation_folder, self.sim.lammps_dump_prefix, current_time),
                    cryst_cutoff = self.sim.cryst_cutoff
                )

                # with open(cryst_array_string, "a") as file:
                #     file.write(f"{current_time} {frac_cryst}\n")
                self.sim.domain_analysis.insert_cryst_line_sorted(path_to_cryst_array, current_time, frac_cryst)
        return 0;

    def calc_avg_domain_size(self):
        # Read in crystallisation files
        #path_to_nematic_vectors_folder = "%s/nematic_vectors" %(self.sim.path_to_home_folder)

        cryst_domain_array = pd.DataFrame(np.zeros([self.sim.df_slurm_sim_data.shape[0], 7]), columns = ["time", "crystallinity", "clusters w/ >= 2 members", 
        "   independent clusters", "mean size cryst domains", "crystalline grid elements", "total volume"])
        make_folder(self.sim.path_to_nematic_vectors_folder)
        for i in tqdm(range(0, len(self.sim.list_run_files))):
            #Read file 
            current_polymer = polymer("%s/%s" %(self.sim.path_to_data_folder, self.sim.list_run_files[i]))
            current_time = self.sim.df_slurm_sim_data["Step"].iloc[i]
            try:
                current_cryst_file = pd.read_csv("%s/%s_cryst_time_%s.txt" %(self.sim.path_to_crystallisation_folder, self.sim.lammps_dump_prefix, current_time), sep = " ")
            except FileNotFoundError:
                print(f"File not found, skipping: {current_time}")
                continue
            current_polymer.read_cryst("%s/%s_cryst_time_%s.txt" %(self.sim.path_to_crystallisation_folder, self.sim.lammps_dump_prefix, current_time))
            label_matrix = current_polymer.merge_boxes_2(print_results = True, ndot_cutoff = self.sim.ndot_cutoff)
            np.save("%s/label_map_time_%i.npy" %(self.sim.path_to_nematic_vectors_folder, current_time), label_matrix)
        test = np.array([current_time, current_polymer.results.fraction_crystallinity,
            current_polymer.results.total_number_clusters, current_polymer.results.total_number_independent_clusters, 
            current_polymer.results.mean_cluster_size, current_polymer.results.total_number_crystalline_grid_elements,
            current_polymer.atom_coords.volume])
        cryst_domain_array.iloc[i, :] = test
        cryst_domain_array.to_csv("%s/domain_analysis.txt" %(self.sim.path_to_home_folder), sep = " ", header = True)

    def read_crystallisation(self, path: str = None):
        if path is None:
            path = "%s/%s" %(self.sim.path_to_home_folder, "frac_cryst.txt")
        if Path(path).is_file():
            cryst_file = np.loadtxt(path)
            return cryst_file
        else:
            return None


    def read_avg_domain_size(self, path: str = None):
        if path is None:
            path = "%s/domain_analysis.txt" %(self.sim.path_to_home_folder)
        if Path(path).is_file():
            mean_domain_file = pd.read_csv(path, sep = " ").iloc[:, 1:]
            return mean_domain_file
        else:
            print("File not found")
            return 0;


    def read_domain_at_time(self, time, path: str = None):
        if path is None:
            path = "%s/label_map_time_%i.npy" %(self.sim.path_to_nematic_vectors_folder, time)
        if Path(path).is_file():
            domain = np.load(path)
            return domain
        else:
            print("File not found")
            return 0;


    def get_crossover_point(self):
        """Crossover point is defined by point with"""
        current_poly = self.sim.get_polymer_by_time(0)
        n_atoms = current_poly.atom_coords.n_atoms
        plt.scatter(self.sim.df_slurm_sim_data["Step"] * self.sim.timestep, (n_atoms/self.sim.df_slurm_sim_data["Volume"]).diff(), 
        marker = ".", label = "PVA-%i" %(self.sim.polymer_length))
        plt.xlabel(r"t/$\tau$")
        plt.show()

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

        self.path_to_crystallisation_folder = "%s/crystallisation" %(path_to_home_folder)
        self.path_to_nematic_vectors_folder = "%s/nematic_vectors" %(path_to_home_folder)
        #print(self.list_run_files)
        self.lammps_dump_prefix = self.get_file_prefix()
        #self.domain_analysis = domain_analysis(path_to_data_folder, path_to_home_folder, self.list_run_files, self.df_slurm_sim_data, self.lammps_dump_prefix,
        #    cryst_cutoff = self.cryst_cutoff, ndot_cutoff = self.ndot_cutoff)
        self.domain_analysis = domain_analysis(self)
        self.df_cryst = self.domain_analysis.read_crystallisation()
        self.timestep = 0.005




    def get_polymer_by_time(self, time, cell_length = 2.0):
        try:
            current_row = self.df_slurm_sim_data[self.df_slurm_sim_data["Step"] == time].iloc[0]
        except IndexError:
            raise ValueError("Time not in dataset, choose a different time.")
        current_polymer = polymer("%s/%s_run%i_time_%i.txt"%(self.path_to_data_folder, self.lammps_dump_prefix, current_row["Run"], current_row["StepSequence"]), polymer_length= self.polymer_length, cell_length= cell_length)
        try:
            current_cryst_file = pd.read_csv("%s/%s_cryst_time_%s.txt" %(self.path_to_crystallisation_folder, self.lammps_dump_prefix, time), sep = " ")
            current_polymer.read_cryst("%s/%s_cryst_time_%s.txt" %(self.path_to_crystallisation_folder, self.lammps_dump_prefix, time))
        except FileNotFoundError:
            print(f"Crystallisation data not found, skipping time: {time}")
        return current_polymer

    def get_mutiple_polymers_by_time(self, times: list):
        polymer_list = []
        for time in times: 
            current_polymer = self.get_polymer_by_time(time)
            polymer_list.append(current_polymer)
            try:
                current_cryst_file = pd.read_csv("%s/%s_cryst_time_%s.txt" %(self.path_to_crystallisation_folder, self.lammps_dump_prefix, time), sep = " ")
                current_polymer.read_cryst("%s/%s_cryst_time_%s.txt" %(self.path_to_crystallisation_folder, self.lammps_dump_prefix, time))
            except FileNotFoundError:
                print(f"File not found, skipping: {current_time}")
        return polymer_list

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


def calc_crystallisation_and_avg_domain_size(polymer):
    polymer.domain_analysis.calc_crystallisation()
    polymer.domain_analysis.calc_avg_domain_size()

    return 0;

def debug_merge_boxes(polymer):
    print(polymer.atom_coords.bond_vectors["nx"].max())
    print(polymer.atom_coords.bond_vectors["nz"].max())
    #print(33*33*33)
    print("For PVA-%i:" %polymer.atom_coords.polymer_length)
    polymer.merge_boxes_2(print_results= True)

def main():
    #PVA_200 = Simulation(200, "../../data/PVA-200/equil", "../data_online/PVA-200/icryst_T088_Tdot_e-3")
    #PVA_50 = Simulation(50, "../../data/PVA-50/equil", "../data_online/PVA-50/icryst_T088_Tdot_e-3")
    PVA_100 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3")
    # PVA_200 = Simulation(200, "../../data/PVA-200/equil", "../data_online/PVA-200/icryst_T088_Tdot_e-3")
    # PVA_300 = Simulation(300, "../../data/PVA-300/equil", "../data_online/PVA-300/icryst_T088_Tdot_e-3")
    # PVA_500 = Simulation(500, "../../data/PVA-500/equil", "../data_online/PVA-500/icryst_T088_Tdot_e-3")
    # PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")


    PVA_100.domain_analysis.get_crossover_point()
    # current_poly = PVA_1000.get_polymer_by_time(12000000, cell_length= 2.0)
    # debug_merge_boxes(current_poly)
    #PVA_200.get_polymer_by_time(10*120000).merge_boxes_2(print_results = True)
    # PVA_300.get_polymer_by_time(10*120000).merge_boxes_2(print_results = True)
    # PVA_500.get_polymer_by_time(10*120000).merge_boxes_2(print_results = True)
    # PVA_1000.get_polymer_by_time(10*120000).merge_boxes_2(print_results = True)
    #calc_crystallisation_and_avg_domain_size(PVA_50)
    #calc_crystallisation_and_avg_domain_size(PVA_100)
    #calc_crystallisation_and_avg_domain_size(PVA_300)
    #calc_crystallisation_and_avg_domain_size(PVA_500)
    #calc_crystallisation_and_avg_domain_size(PVA_1000)

    #PVA_200.domain_analysis.calc_crystallisation()
    #PVA_200.domain_analysis.calc_avg_domain_size()
    #PVA_200.get_polymer_by_time(1200000*100)
    #PVA_200.domain_analysis.read_avg_domain_size()


    #calc_crystallisation_and_avg_domain_size(PVA_50)
    #calc_crystallisation_and_avg_domain_size(PVA_100)
    #calc_crystallisation_and_avg_domain_size(PVA_200)
    #calc_crystallisation_and_avg_domain_size(PVA_300)
    #calc_crystallisation_and_avg_domain_size(PVA_500)
    #calc_crystallisation_and_avg_domain_size(PVA_1000)

    #current_poly = PVA_200.get_polymer_by_time(time = 21*1200000)
    #current_poly.atom_coords.get_nematic_vector_5()
    #print(current_poly.atom_coords.nridges)
    #current_poly.merge_boxes_2(print_results= True)
if __name__== "__main__":
    main()