from analyse8 import polymer, atom_coords, plot_hk_matrix_2d
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
#from simulation import Simulation, plot_crystallisation_different_polymer_lengths#, avrami_fit, avrami_eq
from simulation import *
from experiment import * 

#TODO run PVA-1000 again

def pva_50_analysis(data_path =  "../../data/PVA-50"):
    icryst_PVA_50_T088 = Simulation(0.88, -3, "%s/slurm-PVA-50_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-50_equil_t_088_tdot_e-3_time" %(data_path), no_runs=3,
        home_folder="../data_online/PVA-50/icryst_T088_Tdot_e-3/sim1", polymer_length=50, home_folder_override= True)
    icryst_PVA_50_T088.calc_crystallisation()
    icryst_PVA_50_T088.calc_avg_domain_size()
    return icryst_PVA_50_T088


def pva_100_analysis(data_path = "../../data/pva-100/quick_quench"):
    icryst_PVA_100_T088 = Simulation(0.88, -3, "%s%s" %(data_path, "/slurm-e3-T088.out"), "%s/equil_t_088_tdot_e-3_time"%(data_path), no_runs = 2, 
        home_folder= "../data_online/PVA-100/icryst_T088_Tdot_e-3/sim1", polymer_length = 100, home_folder_override= True)
    icryst_PVA_100_T088.calc_crystallisation()
    icryst_PVA_100_T088.calc_avg_domain_size()
    return icryst_PVA_100_T088

def pva_200_analysis(data_path =  "../../data/PVA-200"):
    icryst_PVA_200_T088 = Simulation(0.88, -3, "%s/slurm-PVA-200_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-200_equil_t_088_tdot_e-3_time" %(data_path), no_runs=3,
        home_folder="../data_online/PVA-200/icryst_T088_Tdot_e-3/sim1", polymer_length=200, home_folder_override= True)
    icryst_PVA_200_T088.calc_crystallisation()
    icryst_PVA_200_T088.calc_avg_domain_size()
    return icryst_PVA_200_T088



def pva_300_analysis(data_path =  "../../data/PVA-300"):
    icryst_PVA_300_T088 = Simulation(0.88, -3, "%s/slurm-PVA-300_equil_t_088_tdot_e-3_sim1.out" %(data_path), "%s/PVA-300_equil_t_088_tdot_e-3_sim1_time" %(data_path), no_runs=3,
        home_folder="../data_online/PVA-300/icryst_T088_Tdot_e-3/sim1", polymer_length=300, home_folder_override= True)
    icryst_PVA_300_T088.calc_crystallisation()
    icryst_PVA_300_T088.calc_avg_domain_size()
    return icryst_PVA_300_T088

def pva_500_analysis(data_path = "../../data/PVA-500"):
    icryst_PVA_500_T088 = Simulation(0.88, -3, "%s/slurm-PVA-500_equil_t_088_tdot_e-3_sim1.out" %(data_path), "%s/PVA-500_equil_t_088_tdot_e-3_sim1_time" %(data_path), no_runs=4,
        home_folder="../data_online/PVA-500/icryst_T088_Tdot_e-3/sim1", polymer_length=500, home_folder_override= True)
    icryst_PVA_500_T088.calc_crystallisation()
    icryst_PVA_500_T088.calc_avg_domain_size()
    return icryst_PVA_500_T088

def pva_1000_analysis(data_path = "../../data/PVA-1000"):
    icryst_PVA_1000_T088 = Simulation(0.88, -3, "%s/slurm-PVA-1000_equil_t_088_tdot_e-3_sim1.out" %(data_path), "%s/PVA-1000_equil_t_088_tdot_e-3_sim1_time" %(data_path), no_runs=2,
        home_folder="../data_online/PVA-1000/icryst_T088_Tdot_e-3/sim1", polymer_length=1000, home_folder_override= True)
    icryst_PVA_1000_T088.calc_crystallisation()
    icryst_PVA_1000_T088.calc_avg_domain_size()
    return icryst_PVA_1000_T088

def quench_PVA(data_path, slurm_name, files_name, home_folder, poly_length, no_runs = 1):
    quench = Simulation(0.88, -3, "%s/%s" %(data_path, slurm_name), "%s/%s" %(data_path, files_name), home_folder = home_folder, polymer_length= poly_length, home_folder_override= True, 
        no_runs = no_runs)
    quench.calc_crystallisation()
    quench.calc_avg_domain_size()
    return quench

def calc_normalised_entanglement_length(ppa, poly_length):
    #ppa_mean = np.mean(ppa)
    denom = poly_length/ppa
    print(np.std(ppa))
    #ne = np.mean(poly_length / (denom-1))
    ne = np.mean(poly_length/ppa)
    err = (poly_length**2/((poly_length - np.mean(ppa))**2) * np.std(ppa))
    print(err)
    return ne, err

def calculate_ppa(icryst, poly_start_num, poly_stop_num):
    #TODO modify this to calculate all pva from all steps
    home_folder = icryst.path_to_home_folder
    ppa_folder = "%s/ppa" %home_folder
    Path(ppa_folder).mkdir(parents = True, exist_ok = True) #Make home directory if not exists
    #polymer_length = icryst.polymer_length
    for i in range(poly_start_num, poly_stop_num):
        poly = icryst.get_polymer_by_count(i)
        ppa = poly.get_entanglement_length()
        np.savetxt("%s/ppa_t_%i.txt" %(ppa_folder, poly.timestep), ppa)

def calculate_ppa_from_file(path_to_file, save_name = None):
        poly = polymer(path_to_file)
        ppa =  poly.get_entanglement_length()
        if save_name != None:
            np.savetxt("%s" %save_name, ppa)
        return ppa


def calc_order_parameter(poly):
    """Taken from 'A Kinetic View on Statistical Physics, Krapivsky et al, CUP 2010. 
        m(x,t) = l^{-d} \sum \sigma"""
    label_matrix= poly.merge_boxes_2(print_results = True)

    #print(label_matrix)
    #print(poly.df_cryst)
    unique_labels, counts = np.unique(label_matrix, return_counts=True) #Labels and how much each label occurs
    #print(unique_labels)

    cryst_lookup = (
        poly.df_cryst
        .set_index(["xid", "yid", "zid"])["cryst_bool"]
    )

    counts_labels = counts[1:]
    unique_labels_2 = unique_labels[1:]
    print(len(counts_labels[counts_labels > 1]))
    rows = []

    for label in unique_labels_2[counts_labels > 1]: 
        i, j, k = np.where(label_matrix == label)
        coords = np.stack([i, j, k], axis=1)  # shape (k, 3)
        coords_tuples = list(map(tuple, coords))
        cryst_vals = cryst_lookup.loc[coords_tuples].to_numpy()
        center = coords.mean(axis=0)  
        xc, yc, zc = center
        #print(coords)
        #print(cryst_vals)
        l = coords.shape[0]**(1/3)*2

        m = l**(-3) * np.sum(cryst_vals)
        rows.append([xc, yc, zc, l, m])

    df_centers = pd.DataFrame(rows, columns=["xc", "yc", "zc", "l", "m"])
    print(df_centers.mean())
    print(df_centers.std())
    return df_centers.mean()["m"]


def calc_order_parameter_loop(length_loop: int, icryst):
    m_array = np.zeros(length_loop)
    for i in range(1, length_loop):
        m = calc_order_parameter(icryst.get_polymer_by_count(i))
        m_array[i] = m 
    print(m_array)
    plt.plot(m_array)
    plt.show()
def main():



    icryst_PVA_50_T088 = pva_50_analysis()
    icryst_PVA_100_T088 = pva_100_analysis()
    icryst_PVA_200_T088 = pva_200_analysis()
    icryst_PVA_300_T088 = pva_300_analysis()
    icryst_PVA_500_T088 = pva_500_analysis()
    icryst_PVA_1000_T088 = pva_1000_analysis()

    # calculate_ppa(icryst_PVA_50_T088, 0, 15)
    # calculate_ppa(icryst_PVA_100_T088, 0, 15)
    # calculate_ppa(icryst_PVA_200_T088, 0, 15)
    # calculate_ppa(icryst_PVA_300_T088, 0, 15)
    # calculate_ppa(icryst_PVA_500_T088, 0, 15)
    # calculate_ppa(icryst_PVA_1000_T088, 0, 15)


    #print(np.mean(ppa_50), np.mean(ppa_200), np.mean(ppa_300), np.mean(ppa_500))
    #calc_order_parameter_loop(45, icryst_PVA_100_T088)
    #calc_order_parameter(poly)
    #plot_hk_matrix_2d(poly, ndot_cutoff = 0.97)
    #poly.atom_coords.make_cell_grid()
    #print(poly.atom_coords.bond_vectors)
    #label_boxes= poly.merge_boxes_2(print_results = True)
    #poly.bin_label_matrix()
    #print(poly.atom_coords.nridges)
    #print(mol_id, closest)
    #poly.atom_coords.make_cell_grid()
    #poly.atom_coords.datapd.to_csv("PVA_1000_T088_poly0_coords.txt", sep = " ")
    #poly.atom_coords.wrapped_monomers.to_csv("PVA_1000_T088_poly0_wrappedcoords.txt", sep = " ")
    #poly.merge_boxes_2(print_results= True)
    #poly.atom_coords.bond_vectors.to_csv("PVA_1000_T088_poly0_bondvecs.txt", sep = " ")
    # #poly.merge_boxes_2(print_results = True)

    simulation_list = [icryst_PVA_50_T088, icryst_PVA_100_T088, icryst_PVA_200_T088, icryst_PVA_300_T088, icryst_PVA_500_T088, icryst_PVA_1000_T088]
    # for simulation in simulation_list: 
    #     poly = simulation.get_polymer_by_count(0)
    #     print(simulation.polymer_length)
    #     #print(poly.atom_coords.dimensions)
    #     #print(poly.atom_coords.nridges)
    #     print(poly.atom_coords.n_atoms)
    #plot_crystallisation(simulation_list, save = True, savestring = "crystallisation_PVA-100-through-1000.pdf", fit_avrami= False)
    plot_crystallisation_different_polymer_lengths(simulation_list, plot_equal_length= False, save= True, savestring = "T088_icryst_cryst-mean_domain_length.pdf")
    #plot_mean_domain_size_indep_clusters(simulation_list, plot_equal_length= False, savestring = "T088_icryst_no_clusters_mean_domain_length.pdf")
    #plot_no_clusters(simulation_list, plot_equal_length= False, savestring= "T088_icryst_only_clusters.pdf")
    #plot_volume_per_monomer(simulation_list, save= True, savestring = "volume_monomer_T088_Tdot_e-3.pdf")

    # quench_200 = quench_PVA("../../data/PVA-200","slurm-PVA-200_quench_t_088_tdot_e-3_sim1.out", "PVA-200_quench_T088_tdot_e-3_time", 
    #     "../data_online/PVA-200/quench_T088_Tdot_e-3/sim1", 200)
    # quench_300 = quench_PVA("../../data/PVA-300","slurm-PVA-300_quench_t_088_tdot_e-3_sim1.out", "PVA-300_quench_T088_tdot_e-3_sim1_time", 
    #     "../data_online/PVA-300/quench_T088_Tdot_e-3/sim1", 300)
    # quench_500 = quench_PVA("../../data/PVA-500","slurm-PVA-500_quench_t_088_tdot_e-3_sim1.out", "PVA-500_quench_T088_tdot_e-3_sim1_time", 
    #     "../data_online/PVA-500/quench_T088_Tdot_e-3/sim1", 500)
    # quench_1000 = quench_PVA("../../data/PVA-1000","slurm-PVA-1000_quench_t_088_tdot_e-3_sim1.out", "PVA-1000_quench_T088_tdot_e-3_time", 
    #     "../data_online/PVA-1000/quench_T088_Tdot_e-3/sim1", 1000)


if __name__ == "__main__":
    main()