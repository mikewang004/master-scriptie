from analyse7 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from simulation import Simulation



# Global variables

ndot_cutoff = 0.97 #Threshold above which the crystalline domains can be merged 
cryst_cutoff = 0.8 #Threshold for a cell to be considered crystalline

#TODO: rewrite all plotting functions to support more sophisticated data structures as implemented in 
# simulation.[py].calc_avg_domain_size 


def avrami_eq(t, a, n, b):
    return a * (t**n) + b

def avrami_fit(t, y):
    popt, pcov = sp.optimize.curve_fit(avrami_eq, t, np.log(1- y), maxfev = 100000)
    return popt, pcov


def plot_avrami(simulation_list: list, save: bool = False, savestring = None, show_plot = False):

    plt.title("Crystallisation as function of time, PVA-100")
    plt.xlabel("time")
    plt.ylabel("fraction of crystallinity")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(0.1, 0.6)
    plt.xlim(10**7, 10**8)
    for simulation in simulation_list:
        used_simulation_cryst = simulation.cryst[:, :]
        print(used_simulation_cryst)
        print(np.log(1 - (used_simulation_cryst[:, 1])))
        popt, pcov = avrami_fit(used_simulation_cryst[:, 0], np.log(1 - (used_simulation_cryst[:, 1])))
        plt.plot(used_simulation_cryst[:, 0], avrami_eq(used_simulation_cryst[:, 0], *popt))
        label = r"T=%.2f, $\dot{T}=10^{%i}$" %(simulation.temp, simulation.cooling_rate)
        plt.scatter(used_simulation_cryst[:, 0], used_simulation_cryst[:, 1], label = label)

    plt.legend()
    if save is True:
        if savestring is None:
            savestring = "plots/10e%i_T%s_crystallinity.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", ""))
        plt.savefig(savestring)
    if show_plot is True:
        plt.show()
    plt.close()


def polymer_density(simulation):
    #TODO: fix this so that two peaks can be seen.
    # Gets 4 different polymers
    nridges, bins = 33,12 
    polymer_ids = np.linspace(0, simulation.no_polymers-1, dtype = int, num = 2)
    #polymer_ids = np.array([0, 5, 8, 10, 130])
    polymer_ids = np.array([130])
    times = simulation.time_temp_array[polymer_ids, 0]
    polymers = simulation.get_polymer_by_count(polymer_ids)
    list_density_dists = []
    for polymer in polymers:
        density_dist = polymer.get_density_dist(nridges = nridges)
        list_density_dists.append(density_dist)
    
    for i in range(len(list_density_dists)):
        density = list_density_dists[i]
        # time_val = times[i]
        # time_str = rf"time = (${time_val:.1e}) \tau$"

        current_cryst = simulation.cryst[polymer_ids[i], 1]
        plt.hist(density, bins = bins)
        #density.plot(kind = "kde")
        #binwidth = density.std()/5
        #plt.hist(density, bins = np.arange(min(density), max(density) + binwidth, binwidth))
        plt.title(r"Local densities, PVA-100, $T=%.2f, \dot{T}=10^{%.1f}$, cryst = %.2f" %(simulation.temp, simulation.cooling_rate, current_cryst))
        plt.xlabel(r"local densities")
        #plt.savefig("plots/10e%i_T%s_local_density_time_%i.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", ""), times[i]))
        plt.show()
        #plt.close()


def plot_crystallisation(simulation_list: list, save: bool = False, savestring = None, show_plot = False, fit_avrami = False):

    plt.title("Crystallisation as function of time, PVA-100")
    plt.xlabel("time")
    plt.ylabel("fraction of crystallinity")
    if fit_avrami == True:
        plt.xscale("log")
        plt.yscale("log")
        plt.ylim(0.1, 0.6)
        plt.xlim(10**7, 10**8)
        for simulation in simulation_list:
            print(simulation.cryst)
            used_simulation_cryst = simulation.cryst[:, :]
            label = r"T=%.2f, $\dot{T}=10^{%i}$" %(simulation.temp, simulation.cooling_rate)
            plt.scatter(used_simulation_cryst[:, 0], used_simulation_cryst[:, 1], label = label)
            popt, pcov = avrami_fit(used_simulation_cryst[:, 0], used_simulation_cryst[:, 1])
            print(popt)
            plt.plot(used_simulation_cryst[:, 0], avrami_eq(used_simulation_cryst[:, 1], *popt))
    else:
        for simulation in simulation_list:
            label = r"T=%.2f, $\dot{T}=10^{%i}$" %(simulation.temp, simulation.cooling_rate)
            data = np.log(1-simulation.cryst[:, 1])
            plt.loglog(simulation.cryst[:, 0], data, label = label)
            plt.scatter(simulation.cryst[:, 0], simulation.cryst[:, 1], label = label)
    plt.legend()
    if save is True:
        if savestring is None:
            savestring = "plots/10e%i_T%s_crystallinity.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", ""))
        plt.savefig(savestring)
    if show_plot is True:
        plt.show()
    plt.close()


def plot_mean_cluster_length(simulation_list: list, save: bool = False, savestring = None, labels_list: list = None):
    counter = 0
    for simulation in simulation_list:
        if labels_list == None:
            label = r"T=%.2f, $\dot{T}=10^{%i}$" %(simulation.temp, simulation.cooling_rate)
        else:
            if len(labels_list) == len(simulation_list):
                label = labels_list[counter]
                counter = counter + 1
            else:
                label = None
        #data = np.log(1-simulation.cryst[:, 1])
        #plt.loglog(simulation.cryst[:, 0], data, label = label)
        plt.scatter(simulation.mean_cluster_length[:, 0], simulation.mean_cluster_length[:, 1], label = label)
    plt.title("Crystallisation as function of time, PVA-100")
    plt.xlabel("time")
    plt.ylabel("mean cluster size")
    plt.legend()
    if save is True:
        if savestring is None:
            savestring = "plots/10e%i_T%s_mean_cluster_length.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", ""))
        plt.savefig(savestring)
    plt.show()


def plot_mean_cluster_vs_crystallisation_single_temp(simulation_list: list, save: bool = False, savestring = None, labels_list: list = None, plot_equal_length: bool = False):
    if plot_equal_length == True:
        length_list =[]

        for simulation in simulation_list:
            length_list.append(simulation.cryst.shape[0])
        length = min(length_list)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # First plot mean cluster length 
    i = 0
    if plot_equal_length == True:
        for simulation in simulation_list:
            # if i == 0:
            #     label_cryst = "crystallinity, constant temperature"
            #     label_domain = "mean cluster size, constant temperature"
            # else:
            #     label_cryst = "crystallinity, instantaneous warmup from T=0.7"
            #     label_domain = "mean cluster size, instantaneous warmup from T=0.7"
            label_cryst = r"crystallinity,  $\dot{T}=10^{%i}$" %(simulation.cooling_rate)
            label_domain = r"mean domain size, $\dot{T}=10^{%i}$" %(simulation.cooling_rate)
            ax1.plot(simulation.cryst[:length, 0], simulation.cryst[:length, 1], 
                label = label_cryst)
            ax2.scatter(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length.iloc[:length, 4], 
                label =  label_domain)
            ax2.scatter(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length.iloc[:length, 3],
                label =  r"number of domains, $\dot{T}=10^{%i}$" %(simulation.cooling_rate))
    else:
        for simulation in simulation_list:
            ax1.plot(simulation.cryst[:, 0], simulation.cryst[:, 1], 
                label =  r"crystallinity $\dot{T}=10^{%i}$" %(simulation.cooling_rate))
            ax2.scatter(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length[:length, 4],
                label =  r"mean domain size $\dot{T}=10^{%i}$" %( simulation.cooling_rate))

    ax1.set_xlabel("time")
    ax1.set_ylabel(r"$\phi$")
    #ax2.set_ylabel(r"count")
    #fig.tight_layout()
    fig.legend(loc = "lower right", bbox_to_anchor=(0.895, 0.115))
    fig.suptitle(r"Crystallinity and mean domain size, T = %.2f" %(simulation.temp))
    fig.savefig("plots/T088_crystallisation_mean_domain_size.pdf")
    plt.show()


def plot_mean_cluster_vs_crystallisation_single_cooling_rate(simulation_list: list, save: bool = False, savestring = None, labels_list: list = None, plot_equal_length: bool = False):
    if plot_equal_length == True:
        length_list =[]

        for simulation in simulation_list:
            length_list.append(simulation.cryst.shape[0])
        length = min(length_list)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # First plot mean cluster length 
    i = 0
    if plot_equal_length == True:
        for simulation in simulation_list:
            if i == 0:
                iso_string = ""
            else:
                iso_string = ",heated from T=0.7"
            ax1.plot(simulation.cryst[:length, 0], simulation.cryst[:length, 1], 
                label = r"crystallinity,  $T = %.2f$%s" %(simulation.temp, iso_string))
            ax2.scatter(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length.iloc[:length, 4],
                label =  r"mean domain size, $T = %.2f$%s" %(simulation.temp, iso_string))

            #i = i + 1
    else:
        for simulation in simulation_list:
            ax1.plot(simulation.cryst[:, 0], simulation.cryst[:, 1], 
                label =  r"crystallinity $T = %.2f$" %(simulation.temp))
            ax2.scatter(simulation.mean_cluster_length.iloc[:, 0], simulation.mean_cluster_length.iloc[:, 4], 
                label =  r"mean domain size $T = %.2f$" %(simulation.temp))

    ax1.set_xlabel("time")
    ax1.set_ylabel(r"$\phi$")
    ax2.set_ylabel(r"mean domain size")
    #fig.tight_layout()
    fig.legend(loc = "lower right", bbox_to_anchor=(0.895, 0.115))
    fig.suptitle(r"Crystallinity and mean domain size, $\dot{T} = 10^{%i}$" %(simulation.cooling_rate))
    fig.savefig("plots/e%i_crystallisation_mean_domain_size.pdf"%(simulation.cooling_rate))
    plt.show()

def plot_mean_cluster_vs_no_clusters_single_cooling_rate(simulation_list: list, save: bool = False, savestring = None, labels_list: list = None, plot_equal_length: bool = False):

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
        ax1.plot(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length.iloc[:length, 3], 
            label = r"independent clusters,  $T = %.2f$" %(simulation.temp))
        ax2.scatter(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length.iloc[:length, 4],
            label =  r"mean domain size, $T = %.2f$" %(simulation.temp))

    ax1.set_xlabel("time")
    ax1.set_ylabel(r"amount of independent clusters")
    ax2.set_ylabel(r"mean domain size")
    #fig.tight_layout()
    fig.legend(loc = "lower right", bbox_to_anchor=(0.895, 0.115))
    fig.suptitle(r"Independent clusters and mean domain size, $\dot{T} = 10^{%i}$" %(simulation.cooling_rate))
    fig.savefig("plots/e%i_no_clusters_mean_domain_size.pdf"%(simulation.cooling_rate))
    plt.show()




def plot_bond_bond_correlation(simulation):
    print(simulation.no_polymers)
    polymer_ids = np.array([0, 5, 8, 10, 45])
    #polymer_ids = np.array([])
    times = simulation.time_temp_array[polymer_ids, 0]
    polymers = simulation.get_polymer_by_count(polymer_ids)
    bond_bond_list = []
    positions_list = []
    for polymer in polymers:
        polymer.bond_bond_correlation()
        bond_bond_list.append(polymer.results.bond_bond_correlation)
        positions_list.append(polymer.atom_coords.positions)
    #first_timestep_e5.bond_bond_correlation()
    # make_plot.scatter_plot(first_timestep_e5.atom_coords.positions, first_timestep_e5.results.bond_bond_correlation, xlabel = "n",
    #     ylabel = r"cos\theta(n)", title = "Distribution of bond-bond correlations, PVA-100", save_string = "plots/bond_bond_corr.pdf")

    # make_plot.scatter_plot(positions_list, bond_bond_list, labels_list, xlabel = "n",
    #     ylabel = r"$cos\theta(n)$", title = r"Distribution of bond-bond correlations, PVA-100, $\dot{T} = 10^{-5}$", save_string = "plots/bond_bond_corr_PVA_100_tdot_e5.pdf",
    #     show_plot = True, marker=".")
    for i in range(len(bond_bond_list)):
        plt.scatter(positions_list[i], bond_bond_list[i], marker = ".", label = "t=%i" %times[i])
    plt.title(r"Bond-bond correlation, PVA-%i, $T=%.2f, \dot{T}=10^{%.1f}$" %(simulation.polymer_length, simulation.temp, simulation.cooling_rate))
    plt.xlabel("monomer position (n)")
    plt.ylabel(r"$cos\theta(n)$")
    plt.ylim((-0.3, 1.1))
    plt.legend()
    plt.savefig("plots/bond_bond_e%i_T%s.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", "")))
    #plt.show()
    plt.close()



def main():
    data_path = "../../data/pva-100/quick_quench/long_run"
    home_folder_path = "../data_online"
    cooling_e3_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T085.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path), no_runs = 3, 
        home_folder= "../data_online/PVA-100/quick_quench")
    cooling_e3_T085.calc_crystallisation()
    cooling_e3_T085.calc_avg_domain_size()


    cooling_e3_T088 = Simulation(0.88, -3, "%s%s" %(data_path, "/slurm-e3-T088.out"), "%s/equil_t_088_tdot_e-3_time"%(data_path), no_runs = 1, 
        home_folder= "../data_online/PVA-100/quick_quench")
    cooling_e3_T088.calc_crystallisation()
    cooling_e3_T088.calc_avg_domain_size()

    cooling_e4_T085 = Simulation(0.85, -4, "%s%s" %(data_path, "/slurm-e4-T085.out"), "%s/equil_t_085_tdot_e-4_time"%(data_path), no_runs = 2, 
        home_folder= "../data_online/PVA-100/quick_quench")
    cooling_e4_T085.calc_crystallisation()
    cooling_e4_T085.calc_avg_domain_size()

    # cooling_e3_T070 = Simulation(0.7, -4, "%s%s" %(data_path, "/slurm-e3-T07.out"), "%s/equil_t_07_tdot_e-3_time"%(data_path), no_runs = 3, 
    #     home_folder= "../data_online/PVA-100/quick_quench")
    # cooling_e3_T070.calc_crystallisation()
    # cooling_e3_T070.calc_avg_domain_size()

    cooling_e3_T07to085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T07to085.out"), "%s/equil_t_07to085_tdot_e-3_time"%(data_path), no_runs = 1, 
        home_folder= "../data_online/PVA-100/quick_quench/e-3/T07to085/run_1", home_folder_override= True)
    cooling_e3_T07to085.calc_crystallisation()
    cooling_e3_T07to085.calc_avg_domain_size()


    simulation_list = [cooling_e3_T085, cooling_e3_T088]
    plot_mean_cluster_vs_crystallisation_single_cooling_rate(simulation_list, plot_equal_length= True)
    #plot_mean_cluster_vs_no_clusters_single_cooling_rate(simulation_list, plot_equal_length= True)

    # plot_mean_cluster_vs_crystallisation_single_temp(simulation_list, plot_equal_length= True, labels_list = [r"T=0.85, \dot{T}=10^{-3}, isothermal", 
    #     r"T=0.85, \dot{T}=10^{-4}, isothermal", r"T=0.85, \dot{T}=10^{-3}, cold crystallsation"])
    #plot_avrami([cooling_e3_T085], show_plot = True)
    # plot_bond_bond_correlation(cooling_e3_T085)
    # #plot_bond_bond_correlation(cooling_e3_T088)
    # plot_bond_bond_correlation(cooling_e4_T085)


    # simulation_list_same_temperature = [cooling_e3_T088]
    # plot_mean_cluster_vs_crystallisation_single_temp(simulation_list_same_temperature, plot_equal_length= True)
    #polymer_density(cooling_e3_T085)

    #plot_mean_cluster_vs_crystallisation_single_cooling_rate([cooling_e3_T085])


if __name__== "__main__":
    main()