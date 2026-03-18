from analyse7 import polymer, atom_coords, plot_hk_matrix_2d
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from simulation import Simulation#, avrami_fit, avrami_eq

ndot_cutoff = 0.97 #Threshold above which the crystalline domains can be merged 
cryst_cutoff = 0.8 #Threshold for a cell to be considered crystalline





def pva_200_analysis(data_path):
    icryst_PVA_200_T088 = Simulation(0.88, -3, "%s/slurm-PVA-200_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-200_equil_t_088_tdot_e-3_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-200/icryst_T088_Tdot_e-3", polymer_length=200, home_folder_override= True)
    icryst_PVA_200_T088.calc_crystallisation()
    icryst_PVA_200_T088.calc_avg_domain_size()
    return icryst_PVA_200_T088

def pva_50_analysis(data_path):
    icryst_PVA_50_T088 = Simulation(0.88, -3, "%s/slurm-PVA-50_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-50_equil_t_088_tdot_e-3_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-50/icryst_T088_Tdot_e-3",polymer_length=50, home_folder_override= True)
    icryst_PVA_50_T088.calc_crystallisation()
    icryst_PVA_50_T088.calc_avg_domain_size()
    return icryst_PVA_50_T088

def pva_100_analysis(data_path):
    icryst_PVA_100_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T085.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path), no_runs = 3, 
        home_folder= "../data_online/PVA-100/quick_quench")
    icryst_PVA_100_T085.calc_crystallisation()
    icryst_PVA_100_T085.calc_avg_domain_size(ndot_cutoff= 0.98)
    return icryst_PVA_100_T085

def get_domain_distribution_polymer(polymer):
    """Returns distribution of domain sizes"""
    label_matrix = polymer.merge_boxes(print_results= True, ndot_cutoff=0.98)
    print(label_matrix)
    #Get distribution by counting label occurances
    label_matrix = label_matrix[label_matrix != 0]
    #print(label_matrix[label_matrix != 0])
    bins = [1,2,3,4,5,7,10,20,30,40,50,60,70,100]
    n, _, _ = plt.hist(label_matrix.flatten(), bins = bins)
    print(n, bins)
    plt.title(r"Crystalline domain size distribution, PVA-100, time %i, $\phi = %.3f$" %(polymer.timestep, polymer.frac_cryst))
    plt.savefig("cryst_domain_distro_pva_100_time_%i.pdf" %polymer.timestep)
    plt.show()
    return 0;


def plot_crystallisation_different_polymer_lengths(simulation_list, save: bool = False, savestring = None, labels_list: list = None, plot_equal_length: bool = False):

    length_list =[]

    for simulation in simulation_list:
        length_list.append(simulation.cryst.shape[0])
        if plot_equal_length == True:
            length = min(length_list)
        else:
            length = max(length_list)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # First plot mean cluster length 
    for simulation in simulation_list:
        # ax1.plot(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length.iloc[:length, 3], 
        #     label = r"independent clusters, PVA-%i" %(simulation.polymer_length))
        label_cryst = r"crystallinity,  PVA-%i" %(simulation.polymer_length)
        ax1.plot(simulation.cryst[:length, 0], simulation.cryst[:length, 1], 
            label = label_cryst)
        ax2.scatter(simulation.mean_cluster_length.iloc[:length, 0], (simulation.mean_cluster_length.iloc[:length, 4])**(1/3) *2,
            label =  r"mean domain size, PVA-%i" %(simulation.polymer_length))
        temp = simulation.temp

    ax1.set_xlabel("time")
    #ax1.set_ylabel(r"amount of independent clusters")
    ax1.set_ylabel(r"$\phi$")
    ax2.set_ylabel(r"length scale ($\sigma$)")
    #fig.tight_layout()
    fig.legend(loc = "lower right", bbox_to_anchor=(0.895, 0.115))
    #fig.suptitle(r"Independent clusters and mean domain size, $\dot{T} = 10^{%i}$" %(simulation.cooling_rate))
    fig.suptitle(r"Isothermal crystallisation, T=%.2f"%(temp))
    fig.savefig("plots/e%i_no_clusters_mean_domain_size.pdf"%(simulation.cooling_rate))
    plt.show()

def avrami_eq(t, n, b, a):
    #return np.log(1-t)
    return a * (t**n) + b

def avrami_fit(t, y):
    popt, pcov = sp.optimize.curve_fit(avrami_eq, t, y, maxfev = 100000)
    return popt, pcov

def plot_avrami(simulation_list: list, save: bool = False, savestring = None, show_plot = False):

    #plt.xscale("log")
    #plt.yscale("log")
    #plt.ylim(0.1, 0.6)
    #plt.xlim(10**7, 10**8)
    for simulation in simulation_list:
        print(simulation.cryst.shape)
        used_simulation_cryst = simulation.cryst[:10, :]# / 12000000
        #phi0 = np.log(1 - (used_simulation_cryst[:, 1]))
        phi0 = used_simulation_cryst[:, 1]
        #print(used_simulation_cryst)
        #print(used_simulation_cryst)
        #print(np.log(1 - (used_simulation_cryst[:, 1])))
        popt, pcov = avrami_fit(used_simulation_cryst[:, 0], phi0)
        print(popt)
        #plt.xscale("log")        # make y-axis logarithmic
        plt.plot(used_simulation_cryst[:, 0], avrami_eq(used_simulation_cryst[:, 0], *popt),
             label=f"n = {popt[0]:.4f}, b = {popt[1]:.4f}, a = {popt[2]:.10f}")
        plt.scatter(used_simulation_cryst[:, 0], phi0, label = "experimental data")

    plt.title(r"$\phi(t)$, PVA-100, T=%.2f, $\dot{T}=10^{%i}$" %(simulation.temp, simulation.cooling_rate))
    plt.xlabel("time")
    plt.ylabel("fraction of crystallinity")

    plt.legend()
    if save is True:
        if savestring is None:
            savestring = "plots/avrami_10e%i_T%s_crystallinity_early.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", ""))
        plt.savefig(savestring)
    if show_plot is True:
        plt.show()
    plt.close()


def main():
    data_path_50 = "../../data/PVA-50"
    data_path_200 = "../../data/PVA-200"
    data_path_100 = "../../data/pva-100/quick_quench/long_run"

    cooling_e4_T085 = Simulation(0.85, -4, "%s%s" %(data_path_100, "/slurm-e4-T085.out"), "%s/equil_t_085_tdot_e-4_time"%(data_path_100), no_runs = 1, 
        home_folder= "../data_online/PVA-100/quick_quench")
    cooling_e4_T085.calc_crystallisation()
    cooling_e4_T085.calc_avg_domain_size()
    #pva_50_analysis(data_path_50)
    icryst_PVA_200_T088 = pva_200_analysis(data_path_200)
    icryst_PVA_50_T088 = pva_50_analysis(data_path_50)
    icryst_PVA_100_T085 = pva_100_analysis(data_path_100)
    current_polymer = icryst_PVA_100_T085.get_polymer_by_count(100)
    print(current_polymer.atom_distribution())
    #print(current_polymer.bond_distribution())
    #current_polymer.atom_coords.get_nematic_vector_5(save_string= "test.txt")
    get_domain_distribution_polymer(current_polymer)
    #plot_hk_matrix_2d(current_polymer)
    #icryst_PVA_100_T088.get_polymer_by_count(99).atom_coords.get_nematic_vector_5()

    #polymer_list = [icryst_PVA_50_T088, icryst_PVA_100_T085, icryst_PVA_200_T088]
    #plot_crystallisation_different_polymer_lengths(polymer_list, plot_equal_length= True, save= True)
    #plot_avrami([icryst_PVA_100_T085], show_plot= True, save = False)


if __name__ == "__main__":
    main()

