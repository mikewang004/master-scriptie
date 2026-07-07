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


class simulation_plots():

    def __init__(self, simulations):
        self.cm_to_in = 1/2.54
        self.max_x = 15.5
        self.max_y = 22
        self.caption_font = 9
        self.save_folder = "plots"

        self.simulations = simulations
        #self.times = times 



        self.simulation_colours = self.fix_colour_schemes(simulations, "plasma")
        #self.times_colours = self.fix_colour_schemes(times, "magma")
        self.std_width = self.max_x/1.5 * plt_cm_to_in
        self.std_height = self.max_y/3 * plt_cm_to_in

        self.times_colour_list = ["0", "1", "tc", "3", "end"]
        self.times_colours = self.fix_colour_schemes(["0", "1", "tc", "3", "end"], "magma")




    def fix_colour_schemes(self, values_list, cmap_name):
        cmap = plt.get_cmap(cmap_name)
        colors = cmap(np.linspace(0, 1, len(values_list)))

        values_colors_dict = {item: col for item, col in zip(values_list, colors)}

        return values_colors_dict


    def plot_crossover_values_vs_chain_length(self, show_plot = True, save_string = "plots/dc_vs_chain_length.pdf"):
        plt.figure(figsize = (self.std_width, self.std_height))

        times = []
        crossovers = []
        polymer_lengths = []
        idx_list = []
        for i in range(0, len(self.simulations)):

            simulation = self.simulations[i]
            poly = simulation.get_polymer_by_time(0)
            time = simulation.df_slurm_sim_data["Step"] * simulation.timestep
            monomer_density = (poly.atom_coords.n_atoms/simulation.df_slurm_sim_data["Volume"])
            idx, xkn, ykn = simulation.domain_analysis.get_crossover_point_cutoff()
            # times.append(time[crossover_index])
            # crossovers.append(monomer_density[crossover_index])
            idx_list.append(idx)
            times.append(xkn); crossovers.append(ykn);
            monomer_density = monomer_density/monomer_density.iloc[-1]
            plt.scatter(simulation.polymer_length, ykn, c = self.simulation_colours[simulation], label=f"PVA-{simulation.polymer_length}")
        plt.xlabel("Chain lengths")
        #plt.ylabel(r"$t_\text{crossover } [t/\tau]$")
        plt.ylabel(r"$n_\text{monomers}/\sigma^3$")
        plt.legend(fontsize=self.caption_font)
        plt.savefig(save_string)
        if show_plot == True:
            plt.show()

    def plot_monomer_density_and_crossover_values(self, show_plot= True):
        width = self.max_x/1.5 * plt_cm_to_in
        height = self.max_y/3 * plt_cm_to_in

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

        for simulation in self.simulations:
            poly = simulation.get_polymer_by_time(0)
            time = simulation.df_slurm_sim_data["Step"] * simulation.timestep
            monomer_density = (poly.atom_coords.n_atoms/simulation.df_slurm_sim_data["Volume"])




            # Top subplot
            ax1.scatter(time, monomer_density, label=f"PVA-{simulation.polymer_length}", color = self.simulation_colours[simulation])

            monomer_density = monomer_density/monomer_density.iloc[-1]

            # Bottom subplot (same plot for now)
            ax2.scatter(time, monomer_density, label=f"PVA-{simulation.polymer_length}", color = self.simulation_colours[simulation])

        for simulation in self.simulations:
            idx, xkn, ykn = simulation.domain_analysis.get_crossover_point_cutoff()
            idx_list.append(idx)
            times.append(xkn)
            crossovers.append(ykn)
            polymer_lengths.append(simulation.polymer_length)

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
        ax2.set_xlabel(r"$t/\tau$", fontsize = self.caption_font)
        ax1.set_ylabel(r"$n_\text{monomers}/\sigma^3$", fontsize = self.caption_font)
        ax2.set_ylabel(r"$n_\text{monomers}/\sigma^3$ (normalised)", fontsize = self.caption_font)

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

        ax1.set_xscale("log")
        ax2.set_xscale("log")
        fig.tight_layout()
        fig.savefig("%s/crossover_density_vs_time_different_chains_subplots_log.pdf" %self.save_folder)
        if show_plot == True:
            plt.show()


    def plot_avg_domain_size(self, savestring = "mean_domain_size.pdf", show_plot = True):
        plt.figure(figsize = (self.std_width*1.25, self.std_height*1.5))
        for i in range(0, len(self.simulations)):
            simulation = self.simulations[i]
            current_domain_file = simulation.domain_analysis.read_avg_domain_size()
            time = simulation.get_simulation_time()
            plt.scatter(time, current_domain_file["mean size cryst domains"]**(1/3), 
                label = "PVA-%i" %(simulation.polymer_length), c= self.simulation_colours[simulation], marker = ".")
        plt.legend(fontsize=self.caption_font)
        plt.xlabel(r"$t/\tau_c$", fontsize = self.caption_font)
        plt.ylabel(r"$l/\sigma$", fontsize = self.caption_font)
        #plt.title("Mean domain size, various chains")
        if savestring is not None:
            plt.savefig("%s/%s" %(self.save_folder, savestring))
        if show_plot == True:
            plt.show()


    def plot_rg_two_polymers_three_times(self, mode = "rg", savestring_default = True, show_plot = True):
        """Mode can either be 'rg' or "re". """
        fig, axes = plt.subplots(
            3, 1,
            figsize=(self.std_width, 3 * self.std_height),
            sharex=True
        )

        PVA_100 = self.simulations[1]; PVA_1000 = self.simulations[-1]

        index_t_start = 0; index_t_end = 117
        times_PVA_100 = [index_t_start, PVA_100.tc_idex, index_t_end]
        times_PVA_1000 = [index_t_start, PVA_1000.tc_idex, index_t_end]
        letter_subplot_list = ["(a)", "(b)", "(c)"]

        for i in range(0, len(axes)):
            print(times_PVA_100[i])
            current_poly_100 = PVA_100.get_polymer_by_time(int(times_PVA_100[i]*1200000))
            current_poly_1000 = PVA_1000.get_polymer_by_time(int(times_PVA_1000[i]*1200000))

            if mode == "rg":
                current_poly_100.gyration_radius()
                current_poly_1000.gyration_radius()
                values, bins, __ = axes[i].hist(current_poly_100.results.gyration_radius_distribution/
                    np.sqrt(current_poly_100.results.mean_gyration_radius), bins=100, color=self.simulation_colours[PVA_100], density = True, histtype = "step", label = "PVA-%i" %PVA_100.polymer_length)
                values, bins, __ = axes[i].hist(current_poly_1000.results.gyration_radius_distribution/
                    np.sqrt(current_poly_1000.results.mean_gyration_radius), bins=100, color=self.simulation_colours[PVA_1000], density = True, histtype = "step", label = "PVA-%i" %PVA_1000.polymer_length)
                axes[i].set_xlabel(r"$R_g/ \langle R_g \rangle$", fontsize = self.caption_font)
                savestring = "plots/rg_pva_100_1000.pdf"

            elif mode == "re":
                current_poly_100.end_to_end_distance()
                current_poly_1000.end_to_end_distance()
                values, bins, __ = axes[i].hist(current_poly_100.results.end_to_end_distribution/
                    (current_poly_100.results.mean_squared_end_to_end), bins=100, color=self.simulation_colours[PVA_100], density = True, histtype = "step", label = "PVA-%i" %PVA_100.polymer_length)
                values, bins, __ = axes[i].hist(current_poly_1000.results.end_to_end_distribution/
                    (current_poly_1000.results.mean_squared_end_to_end), bins=100, color=self.simulation_colours[PVA_1000], density = True, histtype = "step", label = "PVA-%i" %PVA_1000.polymer_length)
                axes[i].set_xlabel(r"$R_e/ \langle R_e \rangle$", fontsize = self.caption_font)
                savestring = "plots/re_pva_100_1000.pdf"
            axes[i].tick_params(labelbottom=True, labelleft=True)
            axes[i].text(0.02, 0.95, letter_subplot_list[i],
                transform=axes[i].transAxes,
                fontsize=plt_caption_font,
                va="top", ha="left")
            axes[i].legend(fontsize = self.caption_font)

        axes[0].set_title(r"$t = 0$")
        axes[1].set_title(r"$t = t_c$")
        axes[2].set_title(r"$t = t_{end}$")

        #plt.legend(fontsize=self.caption_font)
        plt.tight_layout()
        if savestring_default == True:
            plt.savefig(savestring)
        if show_plot == True:
            plt.show()


def main():

    PVA_50 = Simulation(50, "../../data/PVA-50/equil", "../data_online/PVA-50/icryst_T088_Tdot_e-3")
    PVA_100 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3")
    PVA_200 = Simulation(200, "../../data/PVA-200/equil", "../data_online/PVA-200/icryst_T088_Tdot_e-3")
    PVA_300 = Simulation(300, "../../data/PVA-300/equil", "../data_online/PVA-300/icryst_T088_Tdot_e-3")
    PVA_500 = Simulation(500, "../../data/PVA-500/equil", "../data_online/PVA-500/icryst_T088_Tdot_e-3")
    PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")


    simulations = [PVA_50, PVA_100, PVA_200, PVA_300, PVA_500, PVA_1000]

    sp = simulation_plots(simulations)

    #sp.plot_monomer_density_and_crossover_values(show_plot=False)
    sp.plot_rg_two_polymers_three_times(mode = "re")
    #sp.plot_avg_domain_size()
    #sp.plot_crossover_values_vs_chain_length()


if __name__== "__main__":
    main()
