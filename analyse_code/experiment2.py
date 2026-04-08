from analyse8 import polymer, atom_coords, plot_hk_matrix_2d
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
#from simulation import Simulation, plot_crystallisation_different_polymer_lengths#, avrami_fit, avrami_eq
from simulation import *
from experiment import * 

def pva_50_analysis(data_path =  "../../data/PVA-50"):
    icryst_PVA_50_T088 = Simulation(0.88, -3, "%s/slurm-PVA-50_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-50_equil_t_088_tdot_e-3_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-50/icryst_T088_Tdot_e-3/sim1", polymer_length=50, home_folder_override= True)
    icryst_PVA_50_T088.calc_crystallisation()
    icryst_PVA_50_T088.calc_avg_domain_size()
    return icryst_PVA_50_T088


def pva_100_analysis(data_path = "../../data/pva-100/quick_quench/long_run"):
    icryst_PVA_100_T088 = Simulation(0.88, -3, "%s%s" %(data_path, "/slurm-e3-T088.out"), "%s/equil_t_088_tdot_e-3_time"%(data_path), no_runs = 1, 
        home_folder= "../data_online/PVA-100/quick_quench")
    icryst_PVA_100_T088.calc_crystallisation()
    icryst_PVA_100_T088.calc_avg_domain_size()
    return icryst_PVA_100_T088

def pva_200_analysis(data_path =  "../../data/PVA-200"):
    icryst_PVA_200_T088 = Simulation(0.88, -3, "%s/slurm-PVA-200_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-200_equil_t_088_tdot_e-3_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-200/icryst_T088_Tdot_e-3/sim1", polymer_length=200, home_folder_override= True)
    icryst_PVA_200_T088.calc_crystallisation()
    icryst_PVA_200_T088.calc_avg_domain_size()
    return icryst_PVA_200_T088



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

def quench_PVA(data_path, slurm_name, files_name, home_folder, poly_length):
    quench = Simulation(0.88, -3, "%s/%s" %(data_path, slurm_name), "%s/%s" %(data_path, files_name), home_folder = home_folder, polymer_length= poly_length, home_folder_override= True)
    quench.calc_crystallisation()
    quench.calc_avg_domain_size()

def main():
    #icryst_PVA_300_T088 = pva_300_analysis()
    #icryst_PVA_100_T088 = pva_100_analysis()
    #icryst_PVA_500_T088 = pva_500_analysis()
    #icryst_PVA_50_T088 = pva_50_analysis()
    #icryst_PVA_200_T088 = pva_200_analysis()
    icryst_PVA_1000_T088 = pva_1000_analysis()


    #icryst_PVA_300_T088 = pva_300_analysis()
    poly = icryst_PVA_1000_T088.get_polymer_by_count(20)
    plot_hk_matrix_2d(poly, ndot_cutoff = 0.97)
    #poly.merge_boxes_2()
    #poly.bin_label_matrix()
    #print(poly.atom_coords.nridges)
    #print(mol_id, closest)
    #poly.atom_coords.make_cell_grid()
    #poly.atom_coords.datapd.to_csv("PVA_1000_T088_poly0_coords.txt", sep = " ")
    #poly.atom_coords.wrapped_monomers.to_csv("PVA_1000_T088_poly0_wrappedcoords.txt", sep = " ")
    #poly.merge_boxes_2(print_results= True)
    #poly.atom_coords.bond_vectors.to_csv("PVA_1000_T088_poly0_bondvecs.txt", sep = " ")
    # #poly.merge_boxes_2(print_results = True)

    #simulation_list = [icryst_PVA_50_T088, icryst_PVA_100_T088, icryst_PVA_200_T088, icryst_PVA_300_T088, icryst_PVA_500_T088, icryst_PVA_1000_T088]
    # for simulation in simulation_list: 
    #     poly = simulation.get_polymer_by_count(0)
    #     print(simulation.polymer_length)
    #     #print(poly.atom_coords.dimensions)
    #     #print(poly.atom_coords.nridges)
    #     print(poly.atom_coords.n_atoms)
    #plot_crystallisation_different_polymer_lengths(simulation_list, plot_equal_length= False, save= True, savestring = "T088_icryst_cryst-mean_domain_length.pdf")
    #plot_mean_domain_size_indep_clusters(simulation_list, plot_equal_length= False, savestring = "T088_icryst_no_clusters_mean_domain_length.pdf")
    #plot_no_clusters(simulation_list, plot_equal_length= False, savestring= "T088_icryst_only_clusters.pdf")
    #plot_volume_per_monomer(simulation_list, save= True, savestring = "volume_monomer_T088_Tdot_e-3.pdf")

    # quench_PVA("../../data/PVA-300","slurm-PVA-300_quench_t_088_tdot_e-3_sim1.out", "PVA-300_quench_T088_tdot_e-3_sim1_time", 
    #     "../data_online/PVA-300/quench_T088_Tdot_e-3/sim1", 300)
    # quench_PVA("../../data/PVA-200","slurm-PVA-200_quench_t_088_tdot_e-3_sim1.out", "PVA-200_quench_T088_tdot_e-3_time", 
    #     "../data_online/PVA-200/quench_T088_Tdot_e-3/sim1", 200)
    # quench_PVA("../../data/PVA-500","slurm-PVA-500_quench_t_088_tdot_e-3_sim1.out", "PVA-500_quench_T088_tdot_e-3_sim1_time", 
    #     "../data_online/PVA-500/quench_T088_Tdot_e-3/sim1", 500)
    # quench_PVA("../../data/PVA-1000","slurm-PVA-1000_quench_t_088_tdot_e-3_sim1.out", "PVA-1000_quench_T088_tdot_e-3_time", 
    #     "../data_online/PVA-1000/quench_T088_Tdot_e-3/sim1", 1000)


if __name__ == "__main__":
    main()