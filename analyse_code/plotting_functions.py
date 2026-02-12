from analyse7 import polymer, atom_coords
from simulation import Simulation
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os

# Global variables

ndot_cutoff = 0.97 #Threshold above which the crystalline domains can be merged 
cryst_cutoff = 0.8 #Threshold for a cell to be considered crystalline

data_path = "../../data/pva-100/quick_quench/long_run"
home_folder_path = "../data_online"


class Plots():
    def __init__(simulation_list):
        for simulation in simulation_list:
            self.simulation = Simulation
            Path("%s/plots" %self.simulation.home_folder).mkdir(parents = True, exist_ok = True) #Make plots directory if not exists




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



        


def main():
    pass
    # cooling_e3_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T085.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path), no_runs = 3, 
    #     home_folder= "../data_online/PVA-100/quick_quench/T085")
    # cooling_e3_T085.calc_crystallisation()
    # cooling_e3_T085.calc_avg_domain_size()


    # cooling_e3_T088 = Simulation(0.88, -3, "%s%s" %(data_path, "/slurm-e3-T088.out"), "%s/equil_t_088_tdot_e-3_time"%(data_path), no_runs = 1, 
    #     home_folder= "../data_online/PVA-100/quick_quench/T088")
    # cooling_e3_T088.calc_crystallisation()
    # cooling_e3_T088.calc_avg_domain_size()

    # cooling_e4_T085 = Simulation(0.88, -4, "%s%s" %(data_path, "/slurm-e4-T085.out"), "%s/equil_t_085_tdot_e-4_time"%(data_path), no_runs = 2, 
    #     home_folder= "../data_online/PVA-100/quick_quench/T085")
    # cooling_e4_T085.calc_crystallisation()
    # cooling_e4_T085.calc_avg_domain_size()

if __name__== "__main__":
    main()

    
