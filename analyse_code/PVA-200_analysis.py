from analyse7 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from simulation import Simulation

ndot_cutoff = 0.97 #Threshold above which the crystalline domains can be merged 
cryst_cutoff = 0.8 #Threshold for a cell to be considered crystalline

def pva_200_analysis(data_path):
    icryst_PVA_200_T088 = Simulation(0.88, -3, "%s/slurm-PVA-200_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-200_equil_t_088_tdot_e-3_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-200/icryst_T088_Tdot_e-3")
    icryst_PVA_200_T088.calc_crystallisation()
    icryst_PVA_200_T088.calc_avg_domain_size()

def pva_50_analysis(data_path):
    icryst_PVA_50_T088 = Simulation(0.88, -3, "%s/slurm-PVA-50_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-50_equil_t_088_tdot_e-3_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-50/icryst_T088_Tdot_e-3")
    icryst_PVA_50_T088.calc_crystallisation()
    icryst_PVA_50_T088.calc_avg_domain_size()

def main():
    data_path_50 = "../../data/PVA-50"
    data_path_200 = "../../data/PVA-200"
    pva_50_analysis(data_path_50)



if __name__ == "__main__":
    main()

