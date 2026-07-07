from analyse9 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from simulation import Simulation
import scienceplots
from matplotlib.lines import Line2D
import pandas as pd


plt.style.use('science')
plt_cm_to_in = 1/2.54
plt_max_x = 15.5
plt_max_y = 22
plt_caption_font = 9 #pt

plt_colours_chain_length = "viridis"
plt_colours_time = "magma"



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


def plot_2x2_monomer_density(simulation, times, save_string = None, bins = 18):

    fig, axes = plt.subplots(2, 2, figsize=(12, 12), sharex= True, sharey= True)
    axes = axes.ravel()

    for i in range(0,4):
        current_poly = simulation.get_polymer_by_time(times[i])
        monomer_density = current_poly.atom_coords.assign_monomers_to_box()

        triplet_counts = (
            monomer_density.groupby(['nx', 'ny', 'nz'])
            .size()
            .reset_index(name='count')
        )
        print(triplet_counts)
        values, bins, __ = axes[i].hist(triplet_counts["count"], bins=bins, color="steelblue", density = True)
        if times[i] < 1.1e6:
            axes[i].set_title(rf"$t = {times[i]}\ \mathrm{{\tau}}$")
        else:
            axes[i].set_title(rf"$t = {times[i]:.1e}\ \mathrm{{\tau}}$")
        axes[i].tick_params(labelbottom=True, labelleft=True)
        axes[i].set_xlabel(r"$N_\text{monomer}/(2\sigma)^3$")
        


    fig.suptitle(r"Local density, PVA-%i" %(simulation.polymer_length))
    plt.tight_layout()
    if save_string is not None:
        plt.savefig(save_string)
    else:
        plt.show()
        #axes[i].set_title(col, fontsize=10)

def plot_monomer_density(current_polymer):

    monomer_density = current_polymer.atom_coords.assign_monomers_to_box()

    triplet_counts = (
        monomer_density.groupby(['nx', 'ny', 'nz'])
        .size()
        .reset_index(name='count')
    )

    local_density = triplet_counts["count"]/current_polymer.atom_coords.local_volume
    #print(triplet_counts.min(), triplet_counts.max(), triplet_counts.mean())
    mu, sigma = sp.stats.norm.fit(local_density)
    values, bins, __ = plt.hist(local_density, bins = int(triplet_counts["count"].max()-triplet_counts["count"].min()+1), color="steelblue", density = True)
    x_kde = np.linspace(local_density.min(), local_density.max(), 100)
    pdf = sp.stats.norm.pdf(x_kde, mu, sigma)
    plt.plot(x_kde, pdf)
    plt.title("Local density, PVA-%i" %(current_polymer.atom_coords.polymer_length))
    plt.show()

def plot_multiple_monomer_densities(simulation_list, times: list):
    plt.figure(figsize=(8, 5))
    marker_styles = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
    colors = plt.cm.tab10(np.linspace(0, 1, len(times)))  # or any other colormap
    for i in range(0, len(simulation_list)):
        simulation = simulation_list[i]
        for j in range(0, len(times)):
            current_time = times[j]
            current_poly = simulation.get_polymer_by_time(current_time)
            monomer_count = current_poly.atom_coords.assign_monomers_to_box()

            triplet_counts = (
                monomer_count.groupby(['nx', 'ny', 'nz'])
                .size()
                .reset_index(name='count')
            )
            local_density = triplet_counts["count"]/current_poly.atom_coords.local_volume
            mu, sigma = sp.stats.norm.fit(local_density)
            x_kde = np.linspace(local_density.min(),local_density.max(), 100)
            pdf = sp.stats.norm.pdf(x_kde, mu, sigma)
            plt.plot(x_kde, pdf, 
                color = colors[j], linestyle = "-", marker=marker_styles[i % len(marker_styles)], markersize = 3,
                label = r"PVA-%i, time %i $\tau$" %(current_poly.atom_coords.polymer_length, current_time*simulation.timestep))
    #plt.legend(prop={'size': 6})
    plt.legend()
    plt.xlabel(r"monomer count/$V_\text{local}$")
    plt.title("Local density, normalised probability")
    plt.savefig("plots/local_density_pva-100_1000.pdf")
    plt.show()

        



def plot_volume_vs_density(simulations: list):
    for simulation in simulations:
        current_poly = simulation.get_polymer_by_time(0)
        n_atoms = current_poly.atom_coords.n_atoms
        plt.scatter(simulation.df_slurm_sim_data["Step"] * simulation.timestep, n_atoms/simulation.df_slurm_sim_data["Volume"], 
        marker = ".", label = "PVA-%i" %(simulation.polymer_length))
    plt.legend()
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"$n_\text{atoms}/\sigma^3$")
    plt.title("Density of different chains vs time")
    plt.savefig("plots/volume_monomer_time_different_chains.pdf")




def plot_bond_bond_correlation_different_times(simulation, times: list, savestring: str = None):
    plt.figure(figsize=(12, 4))
    for time in times:
        current_poly = simulation.get_polymer_by_time(time)
        bond_bond_corr = current_poly.bond_bond_correlation_2()
        n = np.arange(1, len(bond_bond_corr)+1)
        plt.scatter(n, bond_bond_corr, label = "t = %i" %(time * simulation.timestep), marker = ".")
        print("time %i done !" %time)
    plt.xlabel("n")
    plt.ylabel(r"$cos\theta(n)$")
    plt.title("Bond-bond correlation, PVA-%i" %(simulation.polymer_length))
    plt.legend()
    if savestring is not None:
        plt.savefig(savestring)
    plt.show()
        #axes[i].set

def plot_bond_bond_correlation_different_times_multiple_simulations(simulations: list, times: list, savestring: str = None):
    fig, axes = plt.subplots(2, 1, figsize=(12, 12), sharex= True, sharey= True)
    axes = axes.ravel()

    for i in range(0,len(simulations)):
        simulation = simulations[i]
        for time in times:
            current_poly = simulation.get_polymer_by_time(time)
            positions, diag_means = current_poly.bond_bond_correlation()
            axes[i].scatter(positions, diag_means, label = "t = %i" %(time * simulation.timestep))
        axes[i].set_xlabel("n")
        axes[i].set_ylabel(r"$cos\theta(n)$")
        axes[i].set_title("PVA-%i" %(simulation.polymer_length))
    fig.suptitle("Bond-bond correlation")

                             


    fig.suptitle(r"$R_g$, PVA-%i" %(simulation.polymer_length))
    plt.tight_layout()
    if savestring is not None:
        plt.savefig(savestring)
    plt.show()
        #axes[i].set_title(col, fontsize=10)


def plot_avg_domain_size(simulations: list, savestring = None):
    plt.figure(figsize=(8, 5))
    for i in range(0, len(simulations)):
        simulation = simulations[i]
        current_domain_file = simulation.domain_analysis.read_avg_domain_size()
        print(current_domain_file)

        plt.scatter(current_domain_file["time"]*simulation.timestep, current_domain_file["mean size cryst domains"]**(1/3), label = "PVA-%i" %(simulation.polymer_length))
    plt.legend()
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"$\sigma$")
    plt.title("Mean domain size, various chains")
    if savestring is not None:
        plt.savefig(savestring)
    plt.show()


def plot_distribution_crystalline_domains(simulations: list, times: list, savestring = None):
    plt.figure(figsize=(8, 5))
    marker_styles = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
    colors = plt.cm.tab10(np.linspace(0, 1, len(times)))  # or any other colormap
    for i in range(0, len(simulations)):
        simulation = simulations[i]
        for j in range(0, len(times)):
            time = times[j]
            current_domain = simulation.domain_analysis.read_domain_at_time(time)
            current_poly = simulation.get_polymer_by_time(time)
            #print(current_domain)

            values, counts = np.unique(current_domain, return_counts=True)

            # for v, c in zip(values, counts):
            #     print(f"value {v} occurs {c} times")

            counts = counts[counts > 1] 
            counts = (counts * current_poly.atom_coords.local_volume)**(1/3)
            mu, sigma = sp.stats.norm.fit(counts)

            min_c = counts.min()
            max_c = counts.max()

            nbins = 50  # choose number of bins
            bins = np.logspace(np.log10(min_c), np.log10(max_c), nbins + 1)
            #values, bins, __ = plt.hist(local_density, bins = int(triplet_counts["count"].max()-triplet_counts["count"].min()+1), color="steelblue", density = True)
            # x_kde = np.linspace(counts.min(), counts.max(), 1000)
            # pdf = sp.stats.norm.pdf(x_kde, mu, sigma)
            # plt.plot(x_kde, pdf, label = r"PVA-%i, %i $\tau$" %(simulation.polymer_length, time * simulation.timestep))
            # plt.legend()

            loc, scale = sp.stats.expon.fit(counts, floc=0)
            lambda_hat = 1.0 / scale
            x_exp = np.linspace((min_c), (max_c), 1000)
            pdf_exp = sp.stats.expon.pdf(x_exp, loc=loc, scale=scale)
            plt.plot(x_exp, pdf_exp, color = colors[j], linestyle = "-", marker=marker_styles[i % len(marker_styles)], markersize = 3,
               label = r"PVA-%i, %i $\tau$" %(simulation.polymer_length, time * simulation.timestep))
            #plt.xlim(-1000, 2000)
    #plt.legend()
            #values,_, _ = plt.hist(counts, bins=bins, label = r"PVA-%i, %i $\tau$" %(simulation.polymer_length, time * simulation.timestep))   # or choose e.g. bins=20
    #print(values, bins)
    plt.xscale("log")
    plt.xlabel(r"Domain size $\sigma$")
    plt.ylabel("Frequency")
    #plt.ylim(0, 60)
    plt.legend()
    plt.title("Domain size distribution")
    if savestring is not None:
        plt.savefig(savestring)
    plt.show()


def plot_mean_rg(simulations: list, savestring = None):
    plt.figure(figsize=(8, 5))
    for i in range(0, len(simulations)):
        simulation = simulations[i]
        gyration_radii = []
        #for j in tqdm(range(0, 2)):
        for j in tqdm(range(0, len(simulation.list_run_files))):
            current_time = simulation.df_slurm_sim_data["Step"].iloc[j]
            current_poly = simulation.get_polymer_by_time(current_time)
            current_poly.gyration_radius()
            print(np.sqrt(current_poly.results.mean_gyration_radius))
            gyration_radii.append(np.sqrt(current_poly.results.mean_gyration_radius))
        plt.scatter(simulation.df_slurm_sim_data["Step"]*simulation.timestep, gyration_radii, label = r"PVA-%i" %(simulation.polymer_length))
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"$R_g$")
    plt.title("Mean gyration radius")
    plt.legend()
    if savestring is not None:
        plt.savefig(savestring)
    plt.show()


def plot_distribution_nematic_eigenvalues(simulations: list, times: list, savestring: str = None):
    plt.figure(figsize=(8, 5))
    for i in range(0, len(simulations)):
        simulation = simulations[i]
        for j in range(0, len(times)):
            time = times[j]
            current_poly = simulation.get_polymer_by_time(time)
            cryst_array = pd.read_csv("%s/%s_cryst_time_%s.txt" %(simulation.path_to_crystallisation_folder, simulation.lammps_dump_prefix, time), sep = " ", 
                index_col = "Unnamed: 0")
            print(cryst_array)
            values, bins, _ = plt.hist(cryst_array["cryst_bool"], label = r"PVA-%i, time = %i $\tau$" %(simulation.polymer_length, time*simulation.timestep))
    plt.legend()
    plt.show()

def plot_crossover_values(simulations: list):
    #plt.figure(figsize=(8, 5))
    plt.figure(figsize = (plt_max_x/1.5 * plt_cm_to_in,plt_max_y/3 * plt_cm_to_in))
    times = []
    crossovers = []
    polymer_lengths = []
    idx_list = []
    for i in range(0, len(simulations)):

        simulation = simulations[i]
        poly = simulation.get_polymer_by_time(0)
        time = simulation.df_slurm_sim_data["Step"] * simulation.timestep
        monomer_density = (poly.atom_coords.n_atoms/simulation.df_slurm_sim_data["Volume"])
        idx, xkn, ykn = simulation.domain_analysis.get_crossover_point_cutoff()
        # times.append(time[crossover_index])
        # crossovers.append(monomer_density[crossover_index])
        idx_list.append(idx)
        times.append(xkn); crossovers.append(ykn);
        polymer_lengths.append(simulation.polymer_length)
        monomer_density = monomer_density/monomer_density.iloc[-1]
        # plt.scatter(time, monomer_density, label = "PVA-%i" %simulation.polymer_length)
        # plt.scatter(xkn, ykn, marker = "x", color = "black")

    tc = pd.DataFrame(
        {
            "index": idx_list,
            "time": times,
            "monomer density": crossovers,
        },
        index=polymer_lengths,
    )

    tc.index.name = "polymer_lengths"

    print(tc.head())
    tc.to_csv("../data_online/crossover_times.txt", sep = " ")
    plt.scatter(polymer_lengths, times)
    #plt.xscale("log")
    #plt.yscale("log")
    #plt.scatter(polymer_lengths, times)
    plt.xlabel("Chain lengths")
    #plt.xlabel(r"$t/\tau$")
    #plt.ylabel(r"$n_\text{monomers}/\sigma^3$")
    plt.ylabel(r"$t_\text{crossover } [t/\tau]$")
    plt.legend(fontsize=plt_caption_font)
    plt.savefig("plots/crossover_density_vs_time_different_chains.pdf")
    plt.show()


def plot_monomer_density_and_crossover_values(simulations: list):
    width = plt_max_x/1.5 * plt_cm_to_in
    height = plt_max_y/3 * plt_cm_to_in

    # Two vertical subplots; total figure height = 2 * height
    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(width, 2 * height),
        sharex=True
    )

    times = []
    crossovers = []
    polymer_lengths = []
    idx_list = []

    for simulation in simulations:
        poly = simulation.get_polymer_by_time(0)
        time = simulation.df_slurm_sim_data["Step"] * simulation.timestep
        monomer_density = (poly.atom_coords.n_atoms/simulation.df_slurm_sim_data["Volume"])
        idx, xkn, ykn = simulation.domain_analysis.get_crossover_point_cutoff()

        idx_list.append(idx)
        times.append(xkn)
        crossovers.append(ykn)
        polymer_lengths.append(simulation.polymer_length)



        # Top subplot
        ax1.scatter(time, monomer_density, label=f"PVA-{simulation.polymer_length}")

        monomer_density = monomer_density/monomer_density.iloc[-1]

        # Bottom subplot (same plot for now)
        ax2.scatter(time, monomer_density, label=f"PVA-{simulation.polymer_length}")
        ax2.scatter(xkn, ykn, marker="x", color="black")

    tc = pd.DataFrame(
        {
            "index": idx_list,
            "time": times,
            "monomer density": crossovers,
        },
        index=polymer_lengths,
    )
    tc.index.name = "polymer_lengths"

    print(tc.head())
    tc.to_csv("../data_online/crossover_times.txt", sep=" ")

    # Horizontal lines in both subplots
    ax2.hlines(0.99, time.iloc[0], time.iloc[-1], color="black")

    # Labels & legend
    ax2.set_xlabel(r"$t/\tau$")
    ax1.set_ylabel(r"$n_\text{monomers}/\sigma^3$")
    ax2.set_ylabel(r"$\frac{n_\text{monomers}}{\sigma^3}/(\frac{n_\text{monomers}}{\sigma^3})_\text{max}$")

    ax1.legend(fontsize=plt_caption_font)
    ax2.legend(fontsize=plt_caption_font)

    ax1.text(
        0.02, 0.95, "(a)",
        transform=ax1.transAxes,
        fontsize=plt_caption_font,
        va="top", ha="left"
    )
    ax2.text(
        0.02, 0.95, "(b)",
        transform=ax2.transAxes,
        fontsize=plt_caption_font,
        va="top", ha="left"
    )


    fig.tight_layout()
    fig.savefig("crossover_density_vs_time_different_chains_subplots.pdf")
    plt.show()



def main():

    PVA_50 = Simulation(50, "../../data/PVA-50/equil", "../data_online/PVA-50/icryst_T088_Tdot_e-3")
    PVA_100 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3")
    PVA_200 = Simulation(200, "../../data/PVA-200/equil", "../data_online/PVA-200/icryst_T088_Tdot_e-3")
    PVA_300 = Simulation(300, "../../data/PVA-300/equil", "../data_online/PVA-300/icryst_T088_Tdot_e-3")
    PVA_500 = Simulation(500, "../../data/PVA-500/equil", "../data_online/PVA-500/icryst_T088_Tdot_e-3")
    PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")

    simulations = [PVA_50, PVA_100, PVA_200, PVA_300, PVA_500, PVA_1000]

    print(PVA_100.tc_time)

    plot_monomer_density_and_crossover_values(simulations)
    plot_crossover_values(simulations)
    #plot_volume_vs_density(simulations)

    #PVA_100 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3")
    #PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")

    # #print(PVA_100.df_slurm_sim_data)
    #simulations = [PVA_100, PVA_1000]
    #times = np.array([0, 8, 30])*1200000
    #times = np.array([0, 5, 10, 30])*1200000
    #times = np.array([0, 10, 91])*1200000
    #times = np.array([0])
    #plot_avg_domain_size([PVA_100, PVA_200, PVA_300, PVA_500, PVA_1000], savestring="plots/mean_domain_size.pdf")
    #plot_bond_bond_correlation_different_times(PVA_1000, times, savestring= "plots/bond_bond_correlation_2_PVA_1000.pdf")

    #plot_2x2_gyration_radius(PVA_100, times, save_string="plots/PVA_100_Rg.pdf")
    # plot_2x2_end_end_radius(PVA_100, times, save_string="plots/PVA_100_Re.pdf")
    # plot_2x2_end_end_radius(PVA_1000, times, save_string="plots/PVA_1000_Re.pdf")

    #current_poly = PVA_100.get_polymer_by_time(times[3])
    #plot_2x2_monomer_density(PVA_100, times, save_string="plots/PVA_100_local_density.pdf")
    #plot_2x2_monomer_density(PVA_1000, times, save_string="plots/PVA_1000_local_density.pdf", bins = 18)
    #current_poly.gyration_radius()
    #plot_gyration_radius(current_poly)
        
    #plot_distribution_crystalline_domains([PVA_500], times)#,savestring="plots/cryst_size_dist_pva_500.pdf")
    #plot_mean_rg([PVA_100, PVA_200, PVA_300], savestring = "plots/Rg_vs_time_PVA_100_200_300.pdf")
    #plot_multiple_monomer_densities(simulations, times)

        

    #plot_distribution_nematic_eigenvalues([PVA_100], np.array([0, 5])*1200000)

if __name__== "__main__":
    main()
