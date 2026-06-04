import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
from experiment2 import quench_PVA, calculate_ppa_from_file
from analyse8 import polymer

def calc_normalised_entanglement_length(ppa, poly_length):
    #ppa_mean = np.mean(ppa)
    #denom = poly_length/ppa
    #print(np.std(ppa))
    #ne = np.mean(poly_length / (denom-1))
    #ne = np.mean(poly_length/ppa)
    ne = np.mean(ppa)
    err = (poly_length**2/((poly_length - np.mean(ppa))**2) * np.std(ppa))
    #print(err)
    return ne, err


def entanglement_length_vs_crystallisation(icryst_list, start_num, stop_num):
    for icryst in icryst_list:
        cryst_array = icryst.cryst
        ne_array = np.zeros(stop_num - start_num)
        mean_domain_size = icryst.mean_cryst_domain_size[start_num:stop_num]
        #print(mean_domain_size)
        for i in range(start_num, stop_num):
            poly = icryst.get_polymer_by_count(i)
            ppa = np.loadtxt("%s/ppa/ppa_t_%i.txt" %(icryst.path_to_home_folder, poly.atom_coords.current_timestep))
            ne, err = calc_normalised_entanglement_length(ppa, poly.atom_coords.polymer_length)
            cryst = cryst_array[i, 1]
            ne_array[i] = ne

        #plt.scatter(cryst_array[start_num:stop_num, 1], ne_array, label = "PVA-%i" %icryst.polymer_length)
        plt.scatter(mean_domain_size, cryst_array[start_num:stop_num, 1], label = "PVA-%i" %icryst.polymer_length)
    plt.title(r"$\phi$ vs PPA, T = 0.88")
    plt.legend()
    #plt.xlabel(r"Crystallisation fraction $(\phi)$")
    #plt.ylabel(r"Primitive path length")
    plt.ylabel(r"Crystallisation fraction $(\phi)$")
    plt.xlabel(r"length scale /primative path length")
    plt.savefig("plots/cryst_vs_mean_cluster_size_multiple_polymers.pdf")
    #plt.savefig("plots/PPA_vs_cryst_multiple_polymers.pdf")
    plt.show()



icryst_PVA_50_T088 = quench_PVA("../../data/PVA-50", "slurm-PVA-50_equil_t_088_tdot_e-3.out", "PVA-50_equil_t_088_tdot_e-3_time", "../data_online/PVA-50/icryst_T088_Tdot_e-3/sim1", 50)
icryst_PVA_100_T088 = quench_PVA("../../data/pva-100/quick_quench/long_run", "/slurm-e3-T088.out","equil_t_088_tdot_e-3_time", "../data_online/PVA-100/icryst_T088_Tdot_e-3/sim1", 100,no_runs = 1)
icryst_PVA_200_T088 = quench_PVA("../../data/PVA-200", "slurm-PVA-200_equil_t_088_tdot_e-3.out", "PVA-200_equil_t_088_tdot_e-3_time", "../data_online/PVA-200/icryst_T088_Tdot_e-3/sim1", 200)
icryst_PVA_300_T088 = quench_PVA("../../data/PVA-300", "slurm-PVA-300_equil_t_088_tdot_e-3_sim1.out", "PVA-300_equil_t_088_tdot_e-3_sim1_time", "../data_online/PVA-300/icryst_T088_Tdot_e-3/sim1", 300)
icryst_PVA_500_T088 = quench_PVA("../../data/PVA-500", "slurm-PVA-500_equil_t_088_tdot_e-3_sim1.out", "PVA-500_equil_t_088_tdot_e-3_sim1_time", "../data_online/PVA-500/icryst_T088_Tdot_e-3/sim1", 500, no_runs =2)


#poly  = icryst_PVA_500_T088.get_polymer_by_count(0)
#poly = polymer("../../data/PVA-500/PVA-500_quench_T088_tdot_e-3_sim1_time_0.txt", polymer_length= 500)
#poly  = polymer("../../data/pva-100/quick_quench/quench_tmin_088_tdot_e-3_time_0.txt", polymer_length= 100)
#poly.end_to_end_distance()
#poly.get_kuhn_length()



#entanglement_length_vs_crystallisation([icryst_PVA_50_T088, icryst_PVA_100_T088, icryst_PVA_200_T088, icryst_PVA_300_T088, icryst_PVA_500_T088], 0, 15)


# Calculate entanglement at start 

#ne_50_start = calculate_ppa_from_file("../../data/PVA-50/PVA-50_quench_T088_tdot_e-3_time_0.txt", save_name = "pva_before_quench/PVA_50_ne.txt")
#ne_1000_start = calculate_ppa_from_file("../../data/pva-100/quick_quench/quench_tmin_088_tdot_e-3_time_0.txt", save_name = "pva_before_quench/PVA_100_ne.txt")
#ne_50_start = calculate_ppa_from_file("../../data/PVA-200/PVA-200_quench_T088_tdot_e-3_time_0.txt", save_name = "pva_before_quench/PVA_200_ne.txt")
#ne_50_start = calculate_ppa_from_file("../../data/PVA-300/PVA-300_quench_T088_tdot_e-3_sim1_time_0.txt", save_name = "pva_before_quench/PVA_300_ne.txt")
#ne_50_start = calculate_ppa_from_file("../../data/PVA-500/PVA-500_quench_T088_tdot_e-3_sim1_time_0.txt", save_name = "pva_before_quench/PVA_500_ne.txt")
#ne_50_start = calculate_ppa_from_file("../../data/PVA-1000/PVA-1000_quench_T088_tdot_e-3_time_0.txt", save_name = "pva_before_quench/PVA_1000_ne.txt")



# icryst_PVA_100_T088 = Simulation(0.88, -3, "%s%s" %(data_path, "/slurm-e3-T088.out"), "%s/equil_t_088_tdot_e-3_time"%(data_path), no_runs = 1, 
#         home_folder= "../data_online/PVA-100/icryst_T088_Tdot_e-3/sim1", polymer_length = 100, home_folder_override= True)

# ppa_50 = np.loadtxt("pva_before_quench/PVA_50_ne.txt")
# ppa_100 = np.loadtxt("pva_before_quench/PVA_100_ne.txt")
# ppa_200 = np.loadtxt("pva_before_quench/PVA_200_ne.txt")
# ppa_300 = np.loadtxt("pva_before_quench/PVA_300_ne.txt")
# ppa_500 = np.loadtxt("pva_before_quench/PVA_500_ne.txt")
# ppa_1000 = np.loadtxt("pva_before_quench/PVA_1000_ne.txt")

# ne_50, ne_err_50 = calc_normalised_entanglement_length(ppa_50, 50)
# ne_100, ne_err_100 = calc_normalised_entanglement_length(ppa_100, 100)
# ne_200, ne_err_200 = calc_normalised_entanglement_length(ppa_200, 200)
# ne_300, ne_err_300 = calc_normalised_entanglement_length(ppa_300, 300)
# ne_500, ne_err_500 = calc_normalised_entanglement_length(ppa_500, 500)
# ne_1000, ne_err_1000 = calc_normalised_entanglement_length(ppa_1000, 1000)


# polymer_lengths = np.array([50, 100, 200, 300, 500, 1000])

# #plt.scatter(polymer_lengths, np.array([np.mean(ppa_50), np.mean(ppa_100),np.mean(ppa_200), np.mean(ppa_300),np.mean(ppa_500), np.mean(ppa_1000)]))
# plt.scatter(polymer_lengths, np.array([ne_50, ne_100, ne_200, ne_300, ne_500, ne_1000])),# yerr = np.array([ne_err_50,ne_err_100,ne_err_200,ne_err_300,ne_err_500,ne_err_1000]), fmt = ".")
# plt.title("Entanglement lengths, T = 1.0")
# plt.ylabel("N_e*")
# plt.xlabel("polymer length")
# plt.savefig("PPA_before_quench_different_polymers.pdf")
# plt.show()