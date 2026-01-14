from analyse7 import polymer, atom_coords
from polymer_plots import *

import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt

data_prefix = "../../data/pva-100"

class simulation():
    """Class wrapper for multiple polymers in the same run"""
    def __init__(self, data_prefix, cooling_rate, temp, times, list_paths = None):
        list_polymers = []
        self.data_prefix = data_prefix
        self.cooling_rate = cooling_rate
        self.temp = temp
        self.times = times
        for i in range(len(times)):
            temp_str = str(temp).replace(".","")
            current_polymer = polymer("%s/equil_t_%s_tdot_e-%i_time_%i.txt" %(data_prefix, temp_str, cooling_rate, times[i]))
            list_polymers.append(current_polymer)
        self.polymers = list_polymers

def plot_bond_bond_correlation():
    tdot_e5_t1 = polymer("%s/cooling_tdot_e-5_time_0.txt" %data_prefix)
    tdot_e5_t08 = polymer("%s/cooling_tdot_e-5_time_4000000.txt" %data_prefix)
    tdot_e5_t07 = polymer("%s/cooling_tdot_e-5_time_6000000.txt" %data_prefix)
    tdot_e5_t05 = polymer("%s/cooling_tdot_e-5_time_10000000.txt" %data_prefix)
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
        ylabel = r"$cos\theta(n)$", title = r"Distribution of bond-bond correlations, PVA-100, $\dot{T} = 10^{-5}$", save_string = "plots/bond_bond_corr_PVA_100_tdot_e5.pdf",
        show_plot = True, marker=".")


def plot_end_to_end_distance(tdot_list, temps_list, tdot, time_or_temp = "temp", current_temp = None):
    dist_list = []
    for i in range(len(tdot_list)):
        print(tdot_list[i])
        tdot_list[i].end_to_end_distance()
        if time_or_temp =="time":
            save_str = "time_%s" %tdot_list[i].atom_coords.current_timestep
            time_val = tdot_list[i].atom_coords.current_timestep
            time_str = rf"time = (${time_val:.1e}) \tau$"
        elif time_or_temp == "temp":
            save_str = "T%s" %temps_list
            time_str = "T = %s" %(tdot, temps_list[i])
        else:
            raise Exception("temp_or_time should be either 'time' or 'temp'")
        dist_list.append(tdot_list[i].results.end_to_end_distribution)
        values, bins, _ = plt.hist(tdot_list[i].results.end_to_end_distribution, bins = 200)
        plt.vlines((tdot_list[i].results.mean_squared_end_to_end), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = r"mean end-to-end length = %.2f" %tdot_list[i].results.mean_squared_end_to_end)
        plt.legend()
        plt.xlabel("end-to-end distance")
        if current_temp == None:
            plt.title(r"End-to-end length distribution, PVA-100, $\dot{T}=10^{-%s} \tau$, %s" %(tdot, time_str))
        
            plt.savefig("plots/10e-%s_end_end_distribution_%s.pdf" %(tdot,save_str))
        else:
            temp_str = str(current_temp).replace(".","")
            plt.title(r"$R_e$, PVA-100, T = %.2f, $\dot{T}=10^{-%s} \tau$, %s" %(current_temp,tdot,time_str))
        
            plt.savefig("plots/10e-%s_T%s_end_end_distribution_%s.pdf" %(tdot,temp_str, save_str))
        plt.close()

def plot_gyration_radius(tdot_list, temps_list, tdot, time_or_temp = "temp", current_temp = None):
    """time_or_temp should be "time" or "temp"""
    dist_list = []
    for i in range(len(tdot_list)):
        print(tdot_list[i])
        tdot_list[i].gyration_radius()
        dist_list.append(tdot_list[i].results.gyration_radius_distribution)
        values, bins, _ = plt.hist(tdot_list[i].results.gyration_radius_distribution, bins = 200)
        plt.vlines(np.sqrt(tdot_list[i].results.mean_gyration_radius), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "mean gyration radius = %.4f" %np.sqrt(tdot_list[i].results.mean_gyration_radius))
        plt.legend()
        plt.xlabel(r"gyration radius ($R_g/\sigma$)")
        if time_or_temp =="time":
            save_str = "time_%s" %tdot_list[i].atom_coords.current_timestep
            time_val = tdot_list[i].atom_coords.current_timestep
            time_str = rf"time = (${time_val:.1e}) \tau$"
        elif time_or_temp == "temp":
            save_str = "T%s" %temps_list
            time_str = "T = %s" %(tdot, temps_list[i])
        else:
            raise Exception("temp_or_time should be either 'time' or 'temp'")
        if current_temp == None:
            plt.title(r"Gyration radius distribution, PVA-100, $\dot{T}=10^{-%s} $, %s" %(tdot,time_str))
            plt.savefig("plots/10e-%s_gyration_radius_%s.pdf" %(tdot, save_str))
        else:
            temp_str = str(current_temp).replace(".","")
            plt.title(r"$R_g$, PVA-100, $T = %.2f$, $\dot{T}=10^{-%s}$, %s" %(current_temp,tdot,time_str))
            plt.savefig("plots/10e-%s_T%s_gyration_radius_%s.pdf" %(tdot, temp_str, save_str))
        plt.close()


    
def main():
    #plot_end_to_end_distance()
    # tdot = 4
    # tdot_e5_t1 = polymer("%s/cooling_tdot_e-5_time_0.txt" %data_prefix)
    # tdot_e5_t08 = polymer("%s/cooling_tdot_e-5_time_4000000.txt" %data_prefix)
    # tdot_e5_t07 = polymer("%s/cooling_tdot_e-5_time_6000000.txt" %data_prefix)
    # tdot_e5_t05 = polymer("%s/cooling_tdot_e-5_time_10000000.txt" %data_prefix)

    # tdot_e3_t1 = polymer("%s/cooling_tdot_e-3_time_0.txt" %data_prefix)
    # tdot_e3_t08 = polymer("%s/cooling_tdot_e-3_time_40000.txt" %data_prefix)
    # tdot_e3_t07 = polymer("%s/cooling_tdot_e-3_time_60000.txt" %data_prefix)
    # tdot_e3_t05 = polymer("%s/cooling_tdot_e-3_time_100000.txt" %data_prefix)

    # tdot_e4_t1 = polymer("%s/cooling_tdot_e-4_time_0.txt" %data_prefix)
    # tdot_e4_t08 = polymer("%s/cooling_tdot_e-4_time_400000.txt" %data_prefix)
    # tdot_e4_t07 = polymer("%s/cooling_tdot_e-4_time_600000.txt" %data_prefix)
    # tdot_e4_t05 = polymer("%s/cooling_tdot_e-4_time_1000000.txt" %data_prefix)
    # tdot_e5_list = [tdot_e5_t1, tdot_e5_t08, tdot_e5_t07, tdot_e5_t05]
    # tdot_e3_list = [tdot_e3_t1, tdot_e3_t08, tdot_e3_t07, tdot_e3_t05]
    # tdot_e4_list = [tdot_e4_t1, tdot_e4_t08, tdot_e4_t07, tdot_e4_t05]
    # temps_list = [1,0.8,0.7,0.5]


    times = [0, 1.2e7, 3.6e7, 6.0e7]
    current_temp = 0.85
    tdot = 3
    quick_quench_e3_t085_long_run = simulation("../../data/pva-100/quick_quench/long_run",tdot,current_temp, times).polymers
    
    plot_end_to_end_distance(quick_quench_e3_t085_long_run, times, tdot, time_or_temp = "time", current_temp=current_temp)
    plot_gyration_radius(quick_quench_e3_t085_long_run, times, tdot, time_or_temp = "time",current_temp=current_temp)


    #plot_gyration_radius(tdot_e3_list, temps_list, tdot = 3, time_or_temp = "time")
    #plot_end_to_end_distance(tdot_e3_list, temps_list, tdot = 3, time_or_temp = "time")
    #plot_bond_bond_correlation()


if __name__ == "__main__":
    main()
