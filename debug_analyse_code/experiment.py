from analyse7 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from simulation import Simulation



def polymer_density(simulation):
    # Gets 4 different polymers
    polymer_ids = np.linspace(0, simulation.no_polymers-1, dtype = int, num = 4)
    polymer_ids = np.array([0, 5, 8, 10, 130])
    times = simulation.times[polymer_ids]
    polymers = simulation.get_polymer_by_count(polymer_ids)
    list_density_dists = []
    for polymer in polymers:
        density_dist = polymer.get_density_dist()
        list_density_dists.append(density_dist)
    
    for i in range(len(list_density_dists)):
        density = list_density_dists[i]
        # time_val = times[i]
        # time_str = rf"time = (${time_val:.1e}) \tau$"

        current_cryst = simulation.cryst[polymer_ids[i], 1]
        plt.hist(density, bins = 15)
        plt.title(r"Local densities, PVA-100, $T=%.2f, \dot{T}=10^{%.1f}$, cryst = %.2f" %(simulation.temp, simulation.cooling_rate, current_cryst))
        plt.xlabel(r"local densities")
        plt.savefig("plots/10e%i_T%s_local_density_time_%i.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", ""), times[i]))
        #plt.show()
        plt.close()


# def plot_crystallisation(simulation):


def main():
    data_path = "../../data/pva-100/quick_quench/long_run"
    cooling_e3_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T085.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path), no_runs = 3)

    # cooling_e3_T07 = Simulation(0.7, -3, "%s%s" %(data_path, "/slurm-T=0.7.out"), "%s/equil_t_07_tdot_e-3_time"%(data_path))
    # cooling_e3_T07.calc_crystallisation("../crystallisation/PVA-100/quick_quench/T07/cryst_equil_T07_e-3.txt", "%s/equil_t_07_tdot_e-3"%(data_path))

    
    polymer_density(cooling_e3_T085)
    # cooling_e3_T085.calc_avg_local_density()
    # print(cooling_e3_T085.avg_local_density)
    # plt.scatter(cooling_e3_T085.avg_local_density[:, 0], cooling_e3_T085.avg_local_density[:, 1])
    # plt.show()
    #cooling_e3_T085.calc_crystallisation("../data_online/PVA-100/quick_quench/T085/cryst_equil_T085_e-3.txt", "../data_online/PVA-100/quick_quench/T085/boxes_eigenvalues/equil_t_085_tdot_e-3")

if __name__== "__main__":
    main()