from analyse7 import polymer, atom_coords
from polymer_plots import *

import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt

data_prefix = "../../data/pva-100"

def plot_bond_bond_correlation():
    tdot_e5_t1 = polymer(atom_coords("%s/cooling_tdot_e-5_time_0.txt" %data_prefix))
    tdot_e5_t08 = polymer(atom_coords("%s/cooling_tdot_e-5_time_4000000.txt" %data_prefix))
    tdot_e5_t07 = polymer(atom_coords("%s/cooling_tdot_e-5_time_6000000.txt" %data_prefix))
    tdot_e5_t05 = polymer(atom_coords("%s/cooling_tdot_e-5_time_10000000.txt" %data_prefix))
    tdot_list = [tdot_e5_t1, tdot_e5_t08, tdot_e5_t07, tdot_e5_t05]
    temps_list = [1,0.8,0.7,0.5]
    positions_list = []
    bond_bond_list = []
    labels_list = []
    for i in range(len(tdot_list)):
        print(tdot_list[i])
        tdot_list[i].bond_bond_correlation()
        bond_bond_list.append(tdot_list[i].results.bond_bond_correlation)
        positions_list.append(tdot_list[i].atom_coords.positions)
        labels_list = [r"T = %.1f" % temps for temps in temps_list]
    #first_timestep_e5.bond_bond_correlation()
    # make_plot.scatter_plot(first_timestep_e5.atom_coords.positions, first_timestep_e5.results.bond_bond_correlation, xlabel = "n",
    #     ylabel = r"cos\theta(n)", title = "Distribution of bond-bond correlations, PVA-100", save_string = "plots/bond_bond_corr.pdf")

    make_plot.scatter_plot(positions_list, bond_bond_list, labels_list, xlabel = "n",
        ylabel = r"cos\theta(n)", title = "Distribution of bond-bond correlations, PVA-100", save_string = "plots/bond_bond_corr_PVA_100_tdot_e5.pdf",
        show_plot = True)


def plot_end_to_end_distance(tdot_list, temps_list, tdot):
    dist_list = []
    for i in range(len(tdot_list)):
        print(tdot_list[i])
        tdot_list[i].end_to_end_distance()
        dist_list.append(tdot_list[i].results.end_to_end_distribution)
        values, bins, _ = plt.hist(tdot_list[i].results.end_to_end_distribution, bins = 200)
        plt.vlines((tdot_list[i].results.mean_squared_end_to_end), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = r"mean end-to-end length = %.2f" %tdot_list[i].results.mean_squared_end_to_end)
        plt.legend()
        plt.title(r"End-to-end length distribution, PVA-100, $\dot{T}=10^{-%s}, T = %s$" %(tdot, temps_list[i]))
        plt.xlabel("squared end-to-end distance")
        plt.savefig("plots/10e-%s_end_end_distribution_T%s.pdf" %(tdot, temps_list[i]))
        plt.close()

def plot_gyration_radius(tdot_list, temps_list, tdot):
    dist_list = []
    for i in range(len(tdot_list)):
        print(tdot_list[i])
        tdot_list[i].gyration_radius()
        dist_list.append(tdot_list[i].results.gyration_radius_distribution)
        values, bins, _ = plt.hist(tdot_list[i].results.gyration_radius_distribution, bins = 200)
        plt.vlines(np.sqrt(tdot_list[i].results.mean_gyration_radius), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "mean gyration radius = %.4f" %np.sqrt(tdot_list[i].results.mean_gyration_radius))
        plt.legend()
        plt.title(r"Gyration radius distribution, PVA-100, $\dot{T}=10^{-%s}, T = %s$" %(tdot, temps_list[i]))
        plt.xlabel("gyration radius")
        plt.savefig("plots/10e-%s_gyration_radius_T%s.pdf" %(tdot, temps_list[i]))
        plt.close()


    
def main():
    #plot_end_to_end_distance()
    tdot = 3
    tdot_e5_t1 = polymer("%s/cooling_tdot_e-5_time_0.txt" %data_prefix)
    tdot_e5_t08 = polymer("%s/cooling_tdot_e-5_time_4000000.txt" %data_prefix)
    tdot_e5_t07 = polymer("%s/cooling_tdot_e-5_time_6000000.txt" %data_prefix)
    tdot_e5_t05 = polymer("%s/cooling_tdot_e-5_time_10000000.txt" %data_prefix)

    tdot_e3_t1 = polymer("%s/cooling_tdot_e-3_time_0.txt" %data_prefix)
    tdot_e3_t08 = polymer("%s/cooling_tdot_e-3_time_40000.txt" %data_prefix)
    tdot_e3_t07 = polymer("%s/cooling_tdot_e-3_time_60000.txt" %data_prefix)
    tdot_e3_t05 = polymer("%s/cooling_tdot_e-3_time_100000.txt" %data_prefix)




    tdot_e5_list = [tdot_e5_t1, tdot_e5_t08, tdot_e5_t07, tdot_e5_t05]
    tdot_e3_list = [tdot_e3_t1, tdot_e3_t08, tdot_e3_t07, tdot_e3_t05]
    temps_list = [1,0.8,0.7,0.5]
    plot_gyration_radius(tdot_e3_list, temps_list, tdot = 3)
    plot_end_to_end_distance(tdot_e3_list, temps_list, tdot = 3)


if __name__ == "__main__":
    main()
