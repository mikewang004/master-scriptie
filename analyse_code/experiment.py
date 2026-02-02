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

def polymer_density(simulation):
    #TODO: fix this so that two peaks can be seen.
    # Gets 4 different polymers
    nridges, bins = 5, 25
    polymer_ids = np.linspace(0, simulation.no_polymers-1, dtype = int, num = 4)
    #polymer_ids = np.array([0, 5, 8, 10, 130])
    polymer_ids = np.array([130])
    times = simulation.times[polymer_ids]
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
        #plt.hist(density, bins = bins, density = True)
        #density.plot(kind = "kde")
        #binwidth = density.std()/5
        #plt.hist(density, bins = np.arange(min(density), max(density) + binwidth, binwidth))
        plt.title(r"Local densities, PVA-100, $T=%.2f, \dot{T}=10^{%.1f}$, cryst = %.2f" %(simulation.temp, simulation.cooling_rate, current_cryst))
        plt.xlabel(r"local densities")
        #plt.savefig("plots/10e%i_T%s_local_density_time_%i.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", ""), times[i]))
        plt.show()
        #plt.close()


def plot_crystallisation(simulation_list: list, save: bool = False, savestring = None):
    for simulation in simulation_list:
        label = r"T=%.2f, $\dot{T}=10^{%i}$" %(simulation.temp, simulation.cooling_rate)
        #data = np.log(1-simulation.cryst[:, 1])
        #plt.loglog(simulation.cryst[:, 0], data, label = label)
        plt.scatter(simulation.cryst[:, 0], simulation.cryst[:, 1], label = label)
        
    plt.legend()
    #plt.set_xscale("log")
    #plt.set_yscale("log")
    #plt.ylim(0.1, 0.6)
    #plt.xlim(10**7, 10**8)
    plt.title("Crystallisation as function of time, PVA-100")
    plt.xlabel("time")
    plt.ylabel("fraction of crystallinity")
    if save is True:
        if savestring is None:
            savestring = "plots/10e%i_T%s_crystallinity.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", ""))
        plt.savefig(savestring)
    plt.close()


def plot_mean_cluster_length(simulation_list: list, save: bool = False, savestring = None, labels_list: list = None):
    counter = 0
    for simulation in simulation_list:
        if labels_list == None:
            label = r"T=%.2f, $\dot{T}=10^{%i}$" %(simulation.temp, simulation.cooling_rate)
        else:
            label = labels_list[counter]
            counter = counter + 1
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
            ax2.scatter(simulation.mean_cluster_length[:length, 0], simulation.mean_cluster_length[:length, 1], 
                label =  label_domain)
    else:
        for simulation in simulation_list:
            ax1.plot(simulation.cryst[:, 0], simulation.cryst[:, 1], 
                label =  r"crystallinity $\dot{T}=10^{%i}$" %(simulation.cooling_rate))
            ax2.scatter(simulation.mean_cluster_length[:, 0], simulation.mean_cluster_length[:, 1], 
                label =  r"mean domain size $\dot{T}=10^{%i}$" %( simulation.cooling_rate))

    ax1.set_xlabel("time")
    ax1.set_ylabel(r"$\phi$")
    ax2.set_ylabel(r"mean domain size")
    #fig.tight_layout()
    fig.legend(loc = "lower right", bbox_to_anchor=(0.895, 0.115))
    fig.suptitle(r"Crystallinity and mean domain size, T = %.2f" %(simulation.temp))
    fig.savefig("plots/T085_crystallisation_mean_domain_size.pdf")
    plt.show()


def plot_mean_cluster_vs_crystallisation_single_cooling_rate(simulation_list: list, save: bool = False, savestring = None, labels_list: list = None, plot_equal_length: bool = False):
    """This does not work with new file structure"""
    if plot_equal_length == True:
        length_list =[]

        for simulation in simulation_list:
            length_list.append(simulation.cryst.shape[0])
        length = min(length_list)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # First plot mean cluster length 
    if plot_equal_length == True:
        for simulation in simulation_list:
            ax1.plot(simulation.cryst[:length, 0], simulation.cryst[:length, 1], 
                label = r"crystallinity,  $T = .2f$" %(simulation.temp))
            ax2.scatter(simulation.mean_cluster_length[:length, 0], simulation.mean_cluster_length[:length, 1], 
                label =  r"mean domain size, $T = .2f$" %(simulation.temp))
    else:
        for simulation in simulation_list:
            ax1.plot(simulation.cryst[:, 0], simulation.cryst[:, 1], 
                label =  r"crystallinity $T = %.2f$" %(simulation.temp))
            ax2.scatter(simulation.mean_cluster_length[:, 0], simulation.mean_cluster_length[:, 1], 
                label =  r"mean domain size $T = %.2f$" %(simulation.temp))

    ax1.set_xlabel("time")
    ax1.set_ylabel(r"$\phi$")
    ax2.set_ylabel(r"mean domain size")
    #fig.tight_layout()
    fig.legend(loc = "lower right", bbox_to_anchor=(0.895, 0.115))
    fig.suptitle(r"Crystallinity and mean domain size, $\dot{T} = %i$" %(simulation.cooling_rate))
    fig.savefig("plots/%i_crystallisation_mean_domain_size.pdf"%(simulation.cooling_rate))
    plt.show()




def main():
    data_path = "../../data/pva-100/quick_quench/long_run"
    home_folder_path = "../data_online"
    cooling_e3_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T085.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path), no_runs = 3, 
        home_folder= "../data_online/PVA-100/quick_quench/T085")
    cooling_e3_T085.calc_crystallisation()
    cooling_e3_T085.calc_avg_domain_size()


    cooling_e3_T088 = Simulation(0.88, -3, "%s%s" %(data_path, "/slurm-e3-T088.out"), "%s/equil_t_088_tdot_e-3_time"%(data_path), no_runs = 1, 
        home_folder= "../data_online/PVA-100/quick_quench/T088")
    cooling_e3_T088.calc_crystallisation()
    cooling_e3_T088.calc_avg_domain_size()

    cooling_e4_T085 = Simulation(0.88, -4, "%s%s" %(data_path, "/slurm-e4-T085.out"), "%s/equil_t_085_tdot_e-4_time"%(data_path), no_runs = 2, 
        home_folder= "../data_online/PVA-100/quick_quench/T085")
    cooling_e4_T085.calc_crystallisation()
    cooling_e4_T085.calc_avg_domain_size()


    simulation_list = [cooling_e3_T085, cooling_e4_T085]
    plot_mean_cluster_vs_crystallisation_single_cooling_rate(simulation_list)


if __name__== "__main__":
    main()