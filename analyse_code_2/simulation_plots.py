from analyse9 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from simulation import Simulation
import scienceplots

plt.style.use('science')



def plot_crystallisation_vs_volume():
    pass

def crystallisation_vs_time(simulations: list):
    for simulation in simulations:
        #plt.scatter(simulation.)
        print(simulation.df_slurm_sim_data)


def plot_gyration_radius(current_poly):

    values, bins, _ = plt.hist((current_poly.results.gyration_radius_distribution), bins = 100, density = True)
    plt.vlines(np.sqrt(current_poly.results.mean_gyration_radius), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "mean gyration radius = %.4f" %np.sqrt(current_poly.results.mean_gyration_radius))
    plt.legend()
    plt.show()


def plot_2x2_gyration_radius(simulation, times, save_string = None):

    fig, axes = plt.subplots(2, 2, figsize=(12, 12), sharex= True, sharey= True)
    axes = axes.ravel()

    for i in range(0,4):
        current_poly = simulation.get_polymer_by_time(times[i])
        current_poly.gyration_radius()
        values, bins, __ = axes[i].hist(current_poly.results.gyration_radius_distribution, bins=100, color="steelblue", density = True)
        axes[i].vlines(np.sqrt(current_poly.results.mean_gyration_radius), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "mean gyration radius = %.4f" %np.sqrt(current_poly.results.mean_gyration_radius))
        axes[i].legend()
        if times[i] < 1.1e6:
            axes[i].set_title(rf"$t = {times[i]}\ \mathrm{{\tau}}$")
        else:
            axes[i].set_title(rf"$t = {times[i]:.1e}\ \mathrm{{\tau}}$")
        axes[i].tick_params(labelbottom=True, labelleft=True)
        axes[i].set_xlabel(r"$R_g$")
        


    fig.suptitle(r"$R_g$, PVA-%i" %(simulation.polymer_length))
    plt.tight_layout()
    if save_string is not None:
        plt.savefig(save_string)
    plt.show()
        #axes[i].set_title(col, fontsize=10)


def plot_2x2_end_end_radius(simulation, times, save_string = None):

    fig, axes = plt.subplots(2, 2, figsize=(12, 12), sharex= True, sharey= True)
    axes = axes.ravel()

    for i in range(0,4):
        current_poly = simulation.get_polymer_by_time(times[i])
        current_poly.end_to_end_distance()
        values, bins, __ = axes[i].hist(current_poly.results.end_to_end_distribution, bins=100, color="steelblue", density = True)
        axes[i].vlines((current_poly.results.mean_squared_end_to_end), ymin = 0, ymax = np.max(values), linestyles ="dashed", color = "red", label = "mean end-to-end radius = %.4f" %np.sqrt(current_poly.results.mean_squared_end_to_end))
        axes[i].legend()
        if times[i] < 1.1e6:
            axes[i].set_title(rf"$t = {times[i]}\ \mathrm{{\tau}}$")
        else:
            axes[i].set_title(rf"$t = {times[i]:.1e}\ \mathrm{{\tau}}$")
        axes[i].tick_params(labelbottom=True, labelleft=True)
        axes[i].set_xlabel(r"$R_e$")
        


    fig.suptitle(r"$R_e$, PVA-%i" %(simulation.polymer_length))
    plt.tight_layout()
    if save_string is not None:
        plt.savefig(save_string)
    else:
        plt.show()
        #axes[i].set_title(col, fontsize=10)

def main():
    PVA_100 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3")
    print(PVA_100.df_slurm_sim_data)
    PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")
    simulations = [PVA_100, PVA_1000]
    times = np.array([0, 5, 10, 30])*1200000
    # current_poly = PVA_100.get_polymer_by_time(times[2])
    # current_poly.gyration_radius()
    # plot_gyration_radius(current_poly)


    #plot_2x2_gyration_radius(PVA_100, times, save_string="plots/PVA_100_Rg.pdf")
    plot_2x2_end_end_radius(PVA_100, times, save_string="plots/PVA_100_Re.pdf")
    plot_2x2_end_end_radius(PVA_1000, times, save_string="plots/PVA_1000_Re.pdf")
        

        

if __name__== "__main__":
    main()