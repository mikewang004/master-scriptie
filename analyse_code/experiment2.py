from analyse8 import polymer, atom_coords, plot_hk_matrix_2d
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from simulation import Simulation, plot_crystallisation_different_polymer_lengths#, avrami_fit, avrami_eq
from experiment import * 


def pva_100_analysis(data_path = "../../data/pva-100/quick_quench/long_run"):
    icryst_PVA_100_T088 = Simulation(0.88, -3, "%s%s" %(data_path, "/slurm-e3-T088.out"), "%s/equil_t_088_tdot_e-3_time"%(data_path), no_runs = 1, 
        home_folder= "../data_online/PVA-100/quick_quench")
    icryst_PVA_100_T088.calc_crystallisation()
    icryst_PVA_100_T088.calc_avg_domain_size(ndot_cutoff= 0.98)
    return icryst_PVA_100_T088



def pva_300_analysis(data_path =  "../../data/PVA-300"):
    icryst_PVA_300_T088 = Simulation(0.88, -3, "%s/slurm-PVA-300_equil_t_088_tdot_e-3_sim1.out" %(data_path), "%s/PVA-300_equil_t_088_tdot_e-3_sim1_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-300/icryst_T088_Tdot_e-3/sim1", polymer_length=300, home_folder_override= True)
    icryst_PVA_300_T088.calc_crystallisation()
    icryst_PVA_300_T088.calc_avg_domain_size()
    return icryst_PVA_300_T088

def pva_500_analysis(data_path = "../../data/PVA-500"):
    icryst_PVA_500_T088 = Simulation(0.88, -3, "%s/slurm-PVA-500_equil_t_088_tdot_e-3_sim1.out" %(data_path), "%s/PVA-500_equil_t_088_tdot_e-3_sim1_time" %(data_path), no_runs=2,
        home_folder="../data_online/PVA-500/icryst_T088_Tdot_e-3/sim1", polymer_length=500, home_folder_override= True)
    icryst_PVA_500_T088.calc_crystallisation()
    icryst_PVA_500_T088.calc_avg_domain_size()
    return icryst_PVA_500_T088

def pva_1000_analysis(data_path = "../../data/PVA-1000"):
    icryst_PVA_1000_T088 = Simulation(0.88, -3, "%s/slurm-PVA-1000_equil_t_088_tdot_e-3_sim1.out" %(data_path), "%s/PVA-1000_equil_t_088_tdot_e-3_sim1_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-1000/icryst_T088_Tdot_e-3/sim1", polymer_length=1000, home_folder_override= True)
    icryst_PVA_1000_T088.calc_crystallisation()
    icryst_PVA_1000_T088.calc_avg_domain_size()
    return icryst_PVA_1000_T088

def main():
    # icryst_PVA_300_T088 = pva_300_analysis()
    # icryst_PVA_500_T088 = pva_500_analysis()
    # icryst_PVA_1000_T088 = pva_1000_analysis()



    icryst_PVA_100_T088 = pva_100_analysis()
    poly = icryst_PVA_100_T088.get_polymer_by_count(20)
    print(poly.atom_coords.bond_vectors)
    poly.atom_coords.bond_vectors.to_csv("PVA_100_T088_poly20_bondvecs.txt", sep = " ")


   # simulation_list = [icryst_PVA_100_T088, icryst_PVA_300_T088, icryst_PVA_500_T088, icryst_PVA_1000_T088]
   # plot_crystallisation_different_polymer_lengths(simulation_list, plot_equal_length= False, save= False)
    #plot_volume_per_monomer(simulation_list, save= True, savestring = "volume_monomer_T088_Tdot_e-3.pdf")


if __name__ == "__main__":
    main()