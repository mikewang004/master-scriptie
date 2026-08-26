from analyse9 import polymer, atom_coords, get_time_temp_from_slurm
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from simulation import Simulation, fit_functions
import scienceplots
from matplotlib.lines import Line2D
import pandas as pd


plt.style.use('science')
plt_cm_to_in = 1/2.54
plt_max_x = 15.5
plt_max_y = 22
plt_caption_font = 9 #pt

plt_colours_chain_length = "viridis"
plt_colours_time = "cividis"


def avrami_fit(t, a, n, b):
    #return a * t**n + b
    return np.log(a) + n*np.log(t) 

class avrami():
    def __init__(self, time, crystallisation):

        self.time = time
        self.phi = crystallisation

        self.log_time = np.log(time)
        self.log_log_phi = np.log(np.log(1/(1-crystallisation))) #returns log(log(1/(1-phi)))
        #self.log_log_phi = np.log(1/(1-crystallisation))




class simulation_plots():

    def __init__(self, simulations, caption_font = 12, legend_font = 9, title_font = 14):
        self.cm_to_in = 1/2.54
        self.max_x = 15.5
        self.max_y = 22
        self.caption_font = caption_font
        self.legend_font = legend_font
        self.title_font = title_font
        self.save_folder = "plots"

        self.simulations = simulations
        #self.times = times 

        plt.rcParams["font.size"] = self.caption_font
        plt.rcParams["legend.fontsize"] = self.legend_font
        plt.rcParams["axes.titlesize"] = self.title_font



        self.simulation_colours = self.fix_colour_schemes(simulations, "plasma")
        #self.times_colours = self.fix_colour_schemes(times, "magma")
        self.std_width = self.max_x/1.5 * plt_cm_to_in
        self.std_height = self.max_y/3 * plt_cm_to_in
        plt.rcParams["figure.figsize"] = [self.std_width, self.std_height]
        self.times_colour_list = ["0", "1", "2", "3", "4"]
        self.times_colours = self.fix_colour_schemes(["0", "1", "2", "3", "4"], "inferno", cmap_max= 0.8)
        self.path_to_latex_plots_folder = "../../master-thesis-latex/content/plots"

        self.end_time_index =  {
            0: 140,  # PVA-50
            1: 140,  # PVA-100
            2: 140,  # PVA-200
            3: 133,  # PVA-300
            4: 99,  # PVA-500
            5: 119,  # PVA-1000
        }




    def fix_colour_schemes(self, values_list, cmap_name, cmap_max = 1):
        cmap = plt.get_cmap(cmap_name)
        colors = cmap(np.linspace(0, cmap_max, len(values_list)))

        values_colors_dict = {item: col for item, col in zip(values_list, colors)}

        return values_colors_dict


    def plot_crossover_values_vs_chain_length(self, show_plot = True, save_string = None):
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
            #idx, xkn, ykn = simulation.domain_analysis.get_crossover_point_cutoff(cutoff=0.985)
            idx, xkn, ykn, popt = simulation.domain_analysis.find_knee(monomer_density)
            # times.append(time[crossover_index])
            # crossovers.append(monomer_density[crossover_index])
            # idx_list.append(idx)
            # times.append(xkn); crossovers.append(ykn);
            mean_cryst = np.mean(simulation.df_cryst[:, 1])
            #monomer_density = monomer_density/monomer_density.iloc[-1]
            #plt.scatter(time[idx], mean_cryst, c = self.simulation_colours[simulation], label=f"PVA-{simulation.polymer_length}")
            plt.scatter(simulation.polymer_length, xkn, c = self.simulation_colours[simulation], label=f"PVA-{simulation.polymer_length}")
        plt.xlabel("Chain lengths")
        plt.ylabel(r"$t_\text{crossover } [t/\tau]$")
        #plt.ylabel(r"$n_\text{monomers}/\sigma^3$")
        plt.legend()
        if save_string == None:
            save_string = "%s/crossover_point/tc_vs_chain_length.pdf" %(self.path_to_latex_plots_folder)
        plt.savefig(save_string)
        if show_plot == True:
            plt.show()


    def plot_monomer_density_and_crossover_values(self, show_plot=True, mode="a", marker_size=None):
        """
        Parameters
        ----------
        show_plot : bool
            Whether to call plt.show() at the end.
        mode : str
            "a" – both subplots (original behaviour)
            "b" – only the monomer-density panel (no fit) as a standalone figure
            "c" – only the density + fit + crossover panel as a standalone figure
        marker_size : float or None
            Marker size for all scatter plots. If None, uses matplotlib's default.
        """

        # ------------------------------------------------------------------ #
        #  Collect data (always needed)                                        #
        # ------------------------------------------------------------------ #
        times           = []
        crossovers      = []
        polymer_lengths = []
        idx_list        = []

        pva_50      = self.simulations[0]
        max_density = (pva_50.get_polymer_by_time(0).atom_coords.n_atoms
                    / pva_50.df_slurm_sim_data["Volume"].iloc[-1])

        # Pre-compute everything so we can reuse it for any mode
        sim_data = []
        for simulation in self.simulations:
            poly             = simulation.get_polymer_by_time(0)
            time             = simulation.df_slurm_sim_data["Step"] * simulation.timestep
            monomer_density  = poly.atom_coords.n_atoms / simulation.df_slurm_sim_data["Volume"]
            idx, xkn, ykn, popt = simulation.domain_analysis.find_knee(monomer_density)

            polymer_lengths.append(simulation.polymer_length)
            times.append(xkn)
            crossovers.append(ykn)
            idx_list.append(idx)

            sim_data.append(dict(
                simulation      = simulation,
                time            = time,
                monomer_density = monomer_density,
                idx             = idx,
                xkn             = xkn,
                ykn             = ykn,
                popt            = popt,
            ))

        tc = pd.DataFrame(
            {"index": idx_list, "time": times, "monomer density": crossovers},
            index=polymer_lengths,
        )
        tc.index.name = "polymer_lengths"
        print(tc.head())
        tc.to_csv("../data_online/crossover_times.txt", sep=" ")

        # ------------------------------------------------------------------ #
        #  Helper: draw panel (a) — raw density scatter                       #
        # ------------------------------------------------------------------ #
        def _draw_panel_a(ax):
            for d in sim_data:
                ax.scatter(
                    d["time"], d["monomer_density"],
                    label=f"PVA-{d['simulation'].polymer_length}",
                    color=self.simulation_colours[d["simulation"]],
                    s=marker_size,
                )
            ax.set_ylabel(r"$n_\text{monomers}/\sigma^3$", fontsize=self.caption_font)
            ax.legend(fontsize=plt_caption_font)

        # ------------------------------------------------------------------ #
        #  Helper: draw panel (c) — density + fit + crossover markers         #
        # ------------------------------------------------------------------ #
        def _draw_panel_c(ax):
            for d in sim_data:
                time_con = np.linspace(0, d["time"].max(), 50000)
                ax.scatter(
                    d["time"], d["monomer_density"],
                    label=f"PVA-{d['simulation'].polymer_length}",
                    color=self.simulation_colours[d["simulation"]],
                    s=marker_size,
                )
                ax.plot(
                    time_con,
                    fit_functions.double_exp(time_con, *d["popt"]),
                    color="r", linestyle="dashed",
                )
            for d in sim_data:
                ax.scatter(d["xkn"], d["ykn"], marker="x", color="black", zorder=3, s=marker_size)
            ax.set_xlabel(r"$t/\tau$",                      fontsize=self.caption_font)
            ax.set_ylabel(r"$n_\text{monomers}/\sigma^3$",  fontsize=self.caption_font)
            ax.legend(fontsize=plt_caption_font)

        # ------------------------------------------------------------------ #
        #  Build figure according to mode                                      #
        # ------------------------------------------------------------------ #
        width  = self.max_x / 1.5 * plt_cm_to_in
        height = self.max_y / 3   * plt_cm_to_in

        if mode == "a":
            # ── original: two stacked subplots ──────────────────────────────
            fig, (ax1, ax2) = plt.subplots(
                2, 1,
                figsize=(width, 2 * height),
                sharex=True,
            )
            _draw_panel_a(ax1)
            _draw_panel_c(ax2)

            ax1.text(0.02, 0.95, "(a)", transform=ax1.transAxes,
                    fontsize=plt_caption_font, va="top", ha="left")
            ax2.text(0.02, 0.95, "(b)", transform=ax2.transAxes,
                    fontsize=plt_caption_font, va="top", ha="left")

            fig.tight_layout()
            fig.savefig(
                "%s/crossover_point/crossover_density_vs_time_different_chains_subplots.pdf"
                % self.path_to_latex_plots_folder
            )

        elif mode == "b":
            # ── standalone panel (a): raw density only ──────────────────────
            fig, ax = plt.subplots(figsize=(width, height))
            _draw_panel_a(ax)
            ax.set_xlabel(r"$t/\tau$", fontsize=self.caption_font)
            fig.tight_layout()
            fig.savefig(
                "%s/crossover_point/crossover_density_vs_time_panel_a.pdf"
                % self.path_to_latex_plots_folder
            )

        elif mode == "c":
            # ── standalone panel (c): density + fit + crossover markers ─────
            fig, ax = plt.subplots(figsize=(width, height))
            _draw_panel_c(ax)
            fig.tight_layout()
            fig.savefig(
                "%s/crossover_point/crossover_density_vs_time_panel_c.pdf"
                % self.path_to_latex_plots_folder
            )

        else:
            raise ValueError(f"Unknown mode '{mode}'. Choose 'a', 'b', or 'c'.")

        if show_plot:
            plt.show()




    def plot_crystallinity(self, savestring = None, show_plot = True):
        plt.figure(figsize = (self.std_width*1.25, self.std_height*1.5))
        for i in range(0, len(self.simulations)):
            simulation = self.simulations[i]
            current_domain_file = simulation.domain_analysis.read_avg_domain_size()
            time = simulation.get_simulation_time()
            

            
            plt.scatter(time[:len(current_domain_file["crystallinity"])], current_domain_file["crystallinity"].iloc[:len(time)], 
                label = "PVA-%i" %(simulation.polymer_length), c= self.simulation_colours[simulation], marker = ".")


        plt.legend(fontsize=self.caption_font)
        plt.xlabel(r"$t_c$", fontsize = self.caption_font)
        plt.ylabel(r"$\phi$", fontsize = self.caption_font)
        #plt.xscale("log")
        #plt.title("Mean domain size, various chains")
        if savestring is None:
            savestring = "%s/domain_analysis/crystallinity_vs_time.pdf" %self.path_to_latex_plots_folder
        plt.savefig("%s" %( savestring))
        if show_plot == True:
            plt.show()

    def plot_crystallinity_avrami(self, savestring = "plots/avrami_crystallinity_vs_time.pdf", show_plot = True):
        plt.figure(figsize = (self.std_width*1.25, self.std_height*1.5))
        for i in range(0, len(self.simulations)):
            simulation = self.simulations[i]
            current_domain_file = simulation.domain_analysis.read_avg_domain_size()
            time = simulation.get_simulation_time()
            av = avrami(time[:len(current_domain_file["crystallinity"])], current_domain_file["crystallinity"].iloc[:len(time)])
            avrami_late_regime_start_index = 20
            av_start_regmine_stop_index = 10
            plt.scatter(av.log_time, av.log_log_phi, 
                label = "PVA-%i" %(simulation.polymer_length), c= self.simulation_colours[simulation], marker = ".")
            popt, pcov = sp.optimize.curve_fit(avrami_fit, av.log_time[avrami_late_regime_start_index:], av.log_log_phi[avrami_late_regime_start_index:], maxfev = 50000)
            #popt, pcov = sp.optimize.curve_fit(avrami_fit, av.log_time[:av_start_regmine_stop_index], av.log_log_phi[:av_start_regmine_stop_index], maxfev = 5000000)
            plt.plot(av.log_time[avrami_late_regime_start_index:], avrami_fit(av.log_time[avrami_late_regime_start_index:], *popt), color = self.simulation_colours[simulation], 
                label = "PVA-%i, a = %.2f, n = %.2f" %(simulation.polymer_length, popt[0], popt[1]))
            # plt.plot(av.log_time[:av_start_regmine_stop_index], avrami_fit(av.log_time[:av_start_regmine_stop_index], *popt), color = self.simulation_colours[simulation], 
            #     label = "PVA-%i, a = %.2f, n = %.2f" %(simulation.polymer_length, popt[0], popt[1]))
            print(simulation.polymer_length)
            print(popt)
        plt.legend(fontsize=self.caption_font)
        plt.xlabel(r"$\log(t/\tau_c)$", fontsize = self.caption_font)
        plt.ylabel(r"$\log(\log(1/(1-\phi$)))", fontsize = self.caption_font)
        #plt.xscale("log")
        #plt.title("Mean domain size, various chains")
        if savestring is not None:
            plt.savefig("%s" %( savestring))
        if show_plot == True:
            plt.show()


    def plot_avg_domain_size(self, savestring = "mean_domain_size.pdf", show_plot = True):
        plt.figure(figsize = (self.std_width*1.25, self.std_height*1.5))
        for i in range(0, len(self.simulations)):
            simulation = self.simulations[i]
            current_domain_file = simulation.domain_analysis.read_avg_domain_size()
            time = simulation.get_simulation_time()
            plt.scatter(time[:len(current_domain_file["mean size cryst domains"])], current_domain_file["mean size cryst domains"].iloc[:len(time)]**(1/3), 
                label = "PVA-%i" %(simulation.polymer_length), c= self.simulation_colours[simulation], marker = ".")


        plt.legend(fontsize=self.caption_font)
        plt.xlabel(r"$t/\tau_c$", fontsize = self.caption_font)
        plt.ylabel(r"$l/\sigma$", fontsize = self.caption_font)
        #plt.xscale("log")
        #plt.title("Mean domain size, various chains")
        if savestring is not None:
            plt.savefig("%s/%s" %(self.save_folder, savestring))
        if show_plot == True:
            plt.show()


    
    def plot_avg_domain_size_and_crossover_values(self, show_plot= True):
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
        # Calc max density (is at PVA-50)
        pva_50 = self.simulations[0]
        max_density = pva_50.get_polymer_by_time(0).atom_coords.n_atoms/pva_50.df_slurm_sim_data["Volume"].iloc[-1]

        for simulation in self.simulations:
            current_domain_file = simulation.domain_analysis.read_avg_domain_size()
            time = simulation.get_simulation_time(normalised_by_tc= False)
            mean_domain_size = current_domain_file["mean size cryst domains"].iloc[:len(time)]**(1/3)
            ax1.scatter(time, mean_domain_size, 
                label = "PVA-%i" %(simulation.polymer_length), c= self.simulation_colours[simulation], marker = ".")

            idx, xkn, ykn, popt = simulation.domain_analysis.find_knee(mean_domain_size)
            polymer_lengths.append(simulation.polymer_length)
            times.append(xkn); crossovers.append(ykn); idx_list.append(idx)
            time_con = np.linspace(0, time.max(), 50000)
            ax2.scatter(time, mean_domain_size, label=f"PVA-{simulation.polymer_length}", color = self.simulation_colours[simulation], marker = ".", s = 4.0)
            ax2.plot(time_con, fit_functions.double_exp(time_con, *popt), color = self.simulation_colours[simulation], linewidth = 3.0)

        tc = pd.DataFrame(
            {
                "index": idx_list,
                "time": times,
                "mean domain size": crossovers,
            },
            index=polymer_lengths,
        )
        tc.index.name = "polymer_lengths"

        for i in range(0,len(times)):
            xkn = times[i]; ykn = crossovers[i]
            ax2.scatter(xkn, ykn, marker="x", color="black", zorder = 3)

        print(tc.head())
        tc.to_csv("../data_online/mean_domain_size.txt", sep=" ")


        # Horizontal lines in both subplots

        # Labels & legend
        ax2.set_xlabel(r"$t/\tau$", fontsize = self.caption_font)
        ax1.set_ylabel(r"$l/\sigma$", fontsize = self.caption_font)
        ax2.set_ylabel(r"$l/\sigma$", fontsize = self.caption_font)

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

        #ax1.set_xscale("log")
        #ax2.set_xscale("log")
        fig.tight_layout()
        fig.savefig("%s/crossover_point/crossover_mean_domain_size_vs_time_different_chains_subplots.pdf" %self.path_to_latex_plots_folder)
        if show_plot == True:
            plt.show()
    

    def choose_two_polymers(self, N1 = 1, N2 = 5):
        """PVA-50: 0, -100: 1, -200: 2, 300: 3, 500: 4, 1000: 5"""

        sim1 = self.simulations[N1]; sim2 = self.simulations[N2]
        index_t_start = 0
        time_by_index = {
            0: 140,  # PVA-50
            1: 140,  # PVA-100
            2: 140,  # PVA-200
            3: 133,  # PVA-300
            4: 119,  # PVA-500
            5: 119,  # PVA-1000
        }

        times_poly_1 = [index_t_start, sim1.tc_idx, time_by_index[N1]]
        times_poly_2 = [index_t_start, sim2.tc_idx, time_by_index[N2]]
        polymer_list = [sim1, sim2]
        times_different_PVA = [times_poly_1, times_poly_2]
        return polymer_list, times_different_PVA



    def plot_rg_two_polymers_three_times(self, index_poly_1 = None, index_poly_2 = None, mode = "rg", savestring_default = True, show_plot = True):
        """Mode can either be 'rg' or "re". 
        Note: PVA-100 has 174 items, PVA-1000 119."""
        #TODO review smoothening method
        fig, axes = plt.subplots(
            2, 1,
            figsize=(1.5 * self.std_width, 2 * self.std_height),
            sharex=True
        )


        if index_poly_1 == None:
            PVA_100 = self.simulations[1]; PVA_1000 = self.simulations[-1]

            index_t_start = 0; index_t_end_1000 = 119; index_t_end_100 = 140 #Also 140 for PVA-50, -100, 200
            times_PVA_100 = [index_t_start, PVA_100.tc_idx, index_t_end_100]
            times_PVA_1000 = [index_t_start, PVA_1000.tc_idx, index_t_end_1000]
            times_different_PVA = [times_PVA_100, times_PVA_1000]
            polymer_list = [PVA_100, PVA_1000]
        else:
            polymer_list, times_different_PVA = self.choose_two_polymers(index_poly_1, index_poly_2)
        letter_subplot_list = ["(a)", "(b)", "(c)"]
        #print(times_PVA_100[i])
        #current_poly_100 = PVA_100.get_polymer_by_time(int(times_PVA_100[i]*1200000))
        #current_poly_1000 = PVA_1000.get_polymer_by_time(int(times_PVA_1000[i]*1200000))
        start_int = 0 
        if mode == "cryst_domain_dist":
            start_int = 1
        for i in range(0, len(times_different_PVA)):
            polymer_times = times_different_PVA[i]
            first_poly = polymer_list[i].get_polymer_by_time(0)
            for j in range(start_int, len(polymer_times)):
            #for j in range(0, 1): 
                current_time = int(polymer_times[j])*1200000
                current_poly = polymer_list[i].get_polymer_by_time(current_time)
                if mode == "rg":
                    first_poly.gyration_radius()
                    current_poly.gyration_radius()
                    # values, bins, __ = axes[i].hist(current_poly.results.gyration_radius_distribution/
                    #     np.sqrt(current_poly.results.mean_gyration_radius), bins=50, 
                    #     color=self.times_colours["%i" %(2*j)], density = True, histtype = "step", 
                    #     label = r"$%i t_c$" %(int(current_poly.atom_coords.current_timestep/polymer_list[i].tc_time)))
                    counts, bin_edges = np.histogram(current_poly.results.gyration_radius_distribution/
                         np.sqrt(first_poly.results.mean_gyration_radius), bins = 80, density = True)

                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) /2

                    smooth_counts = sp.ndimage.gaussian_filter1d(counts, sigma = 1.7)
                    axes[i].plot(bin_centers, smooth_counts, color=self.times_colours["%i" %(2*j)], linestyle = "-", marker = ".",
                    label = r"$%i t_c$" %(int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time)))
                    axes[i].set_xlabel(r"$R_g/ \sqrt{\langle R_{g, t = 0 tc}^2 \rangle}$", fontsize = self.caption_font)
                    axes[i].set_ylabel(r"$P(R_g/ \sqrt{\langle R_{g, t = 0 tc}^2 \rangle})$")
                    savestring = "%s/polymer_conformation/rg_pva_%i_%i.pdf" %(self.path_to_latex_plots_folder, polymer_list[0].polymer_length, polymer_list[1].polymer_length)

                elif mode == "re":
                    first_poly.end_to_end_distance()
                    current_poly.end_to_end_distance()
                    counts, bin_edges = np.histogram(current_poly.results.end_to_end_distribution/first_poly.results.mean_squared_end_to_end, bins = 100, density= True)
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) /2
                    smooth_counts = sp.ndimage.gaussian_filter1d(counts, sigma = 2.0)
                    axes[i].plot(bin_centers, smooth_counts, color=self.times_colours["%i" %(2*j)], linestyle = "-", marker = ".",
                    label = r"$%i t_c$" %(int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time)))
                    axes[i].set_xlabel(r"$R_e/ \sqrt{\langle R_{e, t = 0 tc}^2 \rangle}$", fontsize = self.caption_font)
                    axes[i].set_ylabel(r"$P(R_e/ \sqrt{\langle R_{e, t = 0 tc}^2 \rangle})$")
                    savestring = "%s/polymer_conformation/re_pva_%i_%i.pdf" %(self.path_to_latex_plots_folder, polymer_list[0].polymer_length, polymer_list[1].polymer_length)
                #values, bins, __ = axes[i].hist(current_poly_1000.results.end_to_end_distribution/
                #    (current_poly_1000.results.mean_squared_end_to_end), bins=50, color=self.simulation_colours[PVA_1000], density = True, histtype = "step", label = "PVA-%i" %PVA_1000.polymer_length)
                elif mode == "nematic":

                    current_poly.nematic_distributuion()
                    print("nem dist calc done!")
                    counts, bin_edges = np.histogram(current_poly.results.nematic_value_dist, bins = 100, density = True)
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) /2
                    smooth_counts = sp.ndimage.gaussian_filter1d(counts, sigma = 1.0)
                    axes[i].plot(bin_centers, smooth_counts, color=self.times_colours["%i" %(2*j)], linestyle = "-", marker = ".",
                    label = r"$%i t_c$" %(int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time)))
                # values, bins, __ = axes[i].hist(current_poly_1000.results.nematic_value_dist, 
                #     bins=100, color=self.simulation_colours[PVA_1000], density = True, histtype = "step", label = "PVA-%i" %PVA_1000.polymer_length)
                    axes[i].set_xlabel(r"$\lambda$")
                    axes[i].set_ylabel(r"$P(\lambda)$")
                    
                    savestring = "%s/nematic_dist/nem_value_pva_%i_%i.pdf" %(self.path_to_latex_plots_folder, polymer_list[0].polymer_length, polymer_list[1].polymer_length)


                elif mode == "cryst_domain_dist":
                    cryst_domain_dist = pd.read_csv("%s/domain_distribution/%s_domain_dist_%s.txt" 
                        %(polymer_list[i].path_to_home_folder, polymer_list[i].lammps_dump_prefix, current_time), sep = " ")
                    #print(cryst_domain_dist["volume_pdf"])

                    axes[i].plot(cryst_domain_dist["volume"], cryst_domain_dist["volume_pdf"], color=self.times_colours["%i" %(2*j)], 
                        label = r"$%i t_c$" %(int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time)))
                    savestring = "%s/domain_analysis/nem_value_pva_%i_%i.pdf" %(self.path_to_latex_plots_folder, polymer_list[0].polymer_length, polymer_list[1].polymer_length)
                    axes[i].set_xscale("log")
                    axes[i].set_yscale("log")
                    axes[i].set_xlabel(r"$V_\text{domain}$")
                    axes[i].set_ylabel(r"$d\phi/dV$")


                elif mode == "local_monomer_density_dist":
                    monomer_count = current_poly.atom_coords.assign_monomers_to_box()
                    print(monomer_count.iloc[:, :7])
                    triplet_counts = (
                        monomer_count.groupby(['nx', 'ny', 'nz'])
                        .size()
                        .reset_index(name='count')
                    )

                    n_bins = 26
                    Lx, Ly, Lz = current_poly.atom_coords.boxlengths  # box dimensions
                    coords = monomer_count[['xu', 'yu', 'zu']].values
                    x_min, y_min, z_min = coords.min(axis=0)  # or use known box origin

                    bins_x = np.linspace(x_min, x_min + Lx, n_bins + 1)
                    bins_y = np.linspace(y_min, y_min + Ly, n_bins + 1)
                    bins_z = np.linspace(z_min, z_min + Lz, n_bins + 1)
                    counts, edges = np.histogramdd(coords, bins=(bins_x, bins_y, bins_z))

                    # --- voxel volume ---
                    voxel_vol = (Lx / n_bins) * (Ly / n_bins) * (Lz / n_bins)

                    # --- local density: monomers per unit volume ---
                    local_density = counts / voxel_vol  # shape: (n_bins, n_bins, n_bins)
                    global_density = current_poly.atom_coords.n_atoms/current_poly.atom_coords.volume
                    # print(f"Voxel volume:       {voxel_vol:.4f}")
                    # print(f"Max local density:  {local_density.max():.4f}")
                    # print(f"Mean local density: {local_density.mean():.4f}")
                    # print(f"Global density check: {len(coords) / (Lx * Ly * Lz):.4f}")
                    density_flat = local_density.flatten()
                    density_nonzero = density_flat[density_flat > 0]
                    counts, bin_edges = np.histogram(density_nonzero, bins = n_bins, density = True)
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                    counts_smooth = sp.ndimage.gaussian_filter1d(counts.astype(float), sigma=1)

                    axes[i].plot(bin_centers, counts_smooth, color=self.times_colours["%i" %(2*j)],
                        label=r"$%i t_c$" %(int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time)))
                    if j == 2:
                        axes[i].vlines(global_density, 0, 2.5, color = "red", label = r"$\rho_\text{global}$ = %.2f at %i $t_c$" %(global_density,
                        int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time)), linestyle = "dashed")
                    # axes[i].hist(density_nonzero, bins = n_bins, color = self.times_colours["%i" %(2*j)],
                    #      label = r"$%i t_c$" %(int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time)))

                    print(polymer_list[i].polymer_length, int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time))
                    #local_density = triplet_counts["count"]/current_poly.atom_coords.local_volume
                    #print(np.histogram(local_density, bins = 20))
                    mu, sigma = sp.stats.norm.fit(density_nonzero)
                    x_kde = np.linspace(local_density.min(),local_density.max(), 100)
                    pdf = sp.stats.norm.pdf(x_kde, mu, sigma)
                    # axes[i].plot(x_kde, pdf, 
                    #     color = self.times_colours["%i" %(2*j)], linestyle = "-", markersize = 3,
                    #     label = r"$%i t_c$" %(int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time)))
                    savestring = "%s/polymer_conformation/local_density_pva-100_1000.pdf" %self.path_to_latex_plots_folder
                    axes[i].set_xlabel(r"$N_\text{local monomers}/V_\text{local}$")
                    axes[i].set_ylabel(r"$P(N_\text{local monomers}/V_\text{local}$)")

                elif mode == "bond_bond_corr":
                    bond_bond_corr = current_poly.bond_bond_correlation_2()
                    n = np.arange(1, len(bond_bond_corr)+1)
                    axes[i].scatter(n, bond_bond_corr, marker = ".",
                        label = r"$%i t_c$" %(int(current_poly.atom_coords.current_timestep*polymer_list[i].timestep/polymer_list[i].tc_time)),
                        color=self.times_colours["%i" %(2*j)])
                    savestring = "%s/polymer_conformation/bond_bond_correlation_pva_100_1000.pdf" %self.path_to_latex_plots_folder
                    axes[i].set_xlabel("n")
                    axes[i].set_ylabel(r"$cos\theta(n)$")
                    axes[i].set_xlim((0, 100))
            axes[i].tick_params(labelbottom=True, labelleft=True)
            axes[i].text(0.05, 0.95, letter_subplot_list[i],
                transform=axes[i].transAxes,
                fontsize=plt_caption_font,
                va="top", ha="right")
            ha = "left"
            if mode == "nematic":
                ha = "right"
                axes[i].vlines(self.simulations[1].cryst_cutoff, 0, np.max(counts), color = "red",linestyles = "dotted", label = r"$\lambda_\text{cutoff} = 0.8$")
            axes[i].legend(fontsize = self.caption_font)



        axes[0].set_title("PVA-%i" %(polymer_list[0].polymer_length))
        axes[1].set_title("PVA-%i" %(polymer_list[1].polymer_length))

        axes[0].legend(fontsize=self.caption_font)
        axes[1].legend(fontsize=self.caption_font)
        plt.tight_layout()
        if savestring_default == True:
            plt.savefig(savestring)
        if show_plot == True:
            plt.show()

    def plot_rg_vs_re(self, savestring = None):
        list_times = [0,1,7]
        for j in range(0, 3):
            gyration_radius_list = np.zeros(len(self.simulations))
            end_end_list = np.zeros(len(self.simulations))
            polymer_length = np.zeros(len(self.simulations))
            for i in tqdm(range(0, len(self.simulations))):
                simulation = self.simulations[i]
                #current_poly = simulation.get_polymer_by_time(simulation.df_slurm_sim_data["Step"].iloc[simulation.tc_idx])
                if j == 0:
                    time = 0
                    label = r"$0 t_c$"
                elif j ==1:
                    time = simulation.df_slurm_sim_data["Step"].iloc[simulation.tc_idx]
                    label = r"$1 t_c$"
                elif j == 2:
                    time = simulation.df_slurm_sim_data["Step"].iloc[self.end_time_index[i]]
                    label = r"$7 t_c$"
                current_poly = simulation.get_polymer_by_time(time)
                current_poly.end_to_end_distance()
                current_poly.gyration_radius()
                gyration_radius_list[i] = ((current_poly.results.mean_gyration_radius))
                end_end_list[i] = ((current_poly.results.mean_squared_end_to_end)**2)
                polymer_length[i] = simulation.polymer_length
                if i == 4 and j == 2:
                    gyration_radius_list[i] = np.nan
                    end_end_list[i] = np.nan
                    polymer_length[i] = np.nan
                print(end_end_list[i], gyration_radius_list[i])

            plt.scatter(polymer_length, gyration_radius_list / 6, label = r"$\langle R_g^2 \rangle / 6, %i t_c$ "%(list_times[j]))
            plt.scatter(polymer_length, end_end_list, label = r"$\langle R_e^2 \rangle, %i t_c$" %(list_times[j]))
            #plt.scatter(gyration_radius_list, end_end_list/6, label = label, marker = ".")
        #plt.xlabel(r"$\langle R_g^2 \rangle / 6$")
        #plt.ylabel(r"$\langle R_e^2 \rangle$")
        plt.xlabel(r"N")
        plt.xscale("log")
        plt.yscale("log")
        plt.legend()
        if savestring == None:
            savestring = "%s/polymer_conformation/rg_vs_re_at_tc.pdf" %self.path_to_latex_plots_folder
        plt.savefig(savestring)
            




    def plot_stem_length(self, savestring = None, show_plot = True):
        plt.figure(figsize = (self.std_width*1.25, self.std_height*1.5))

        # simulation = self.simulations
        # print(simulation.df_slurm_sim_data)
        # positions = []
        # min_stem_lengths = []
        #max_time_index = 50
        # for i in range(0, max_time_index):
        #     current_poly = simulation.get_polymer_by_time(simulation.df_slurm_sim_data["Step"].iloc[i])
        #     position, min_stem_length = current_poly.bond_bond_corr_min_value()

        #     positions.append(position)
        #     min_stem_lengths.append(min_stem_length)
        for i in tqdm(range(0, len(self.simulations))):
            simulation = self.simulations[i]
            max_time_index = simulation.df_slurm_sim_data["index"].shape[0]
            positions = []
            min_stem_lengths = []
            for j in tqdm(range(0, max_time_index)):

                current_poly = simulation.get_polymer_by_time(simulation.df_slurm_sim_data["Step"].iloc[j])
                position, min_stem_length = current_poly.bond_bond_corr_min_value("%s/bond_bond_correlation/%s_bond_bond_corr_%s.txt" 
                    %(simulation.path_to_home_folder, simulation.lammps_dump_prefix, current_poly.atom_coords.current_timestep))


                if j > 1:
                    if position < positions[j-1]:
                        position = positions[j-1]
                min_stem_lengths.append(min_stem_length)
                positions.append(position)
            
            time = simulation.get_simulation_time()
            print(time.shape, len(positions))
            plt.plot(time[1:], positions[1:], color=self.simulation_colours[simulation], label = "PVA-%i" %simulation.polymer_length, alpha = 0.5)
        #plt.ylim((0, 50))
        plt.legend()
        plt.ylim((15, 40))
        plt.ylabel("Stem length")
        plt.xlabel("$t_c$")
        if savestring == None:
            savestring = "%s/stem_lengths/stem_lengths_all_polymer_chains.pdf" %(self.path_to_latex_plots_folder)
        plt.savefig(savestring)
        plt.show()


    def plot_crystallinity_different_quench_temps(self):
        data_path = "../../data/pva-100/quick_quench/long_run"
        cryst_path = "../test_run/"
        time_temp_T08 = get_time_temp_from_slurm(data_path + "/slurm-e3-T08.out")
        time_temp_T07 = get_time_temp_from_slurm(data_path + "/slurm-e3-T07-run1.out")
        time_temp_T085 = get_time_temp_from_slurm(data_path + "/slurm-e3-T085-run1.out")
        time_temp_T075 = get_time_temp_from_slurm(data_path + "/slurm-e3-T075.out")
        time_temp_T088 = get_time_temp_from_slurm(data_path + "/slurm-e3-T088.out")

        cryst_T07 = "%s/all_times_cryst_equil_t_07_tdot_e-3.txt" %cryst_path
        cryst_T08 = "%s/all_times_cryst_equil_t_08_tdot_e-3.txt" %cryst_path
        cryst_T075 = "%s/all_times_cryst_equil_t_075_tdot_e-3.txt" %cryst_path
        cryst_T085 = "%s/all_times_cryst_equil_t_085_tdot_e-3.txt" %cryst_path
        cryst_T088 = "%s/all_times_cryst_equil_t_088_tdot_e-3.txt" %cryst_path
        cryst_list = [cryst_T07, cryst_T075, cryst_T08, cryst_T085, cryst_T088]
        colors = {0: 'tab:blue',1: 'tab:red', 2: 'tab:purple', 3: 'tab:green', 4: 'tab:orange'}
        temps = [0.7, 0.75, 0.8, 0.85, 0.88]
        i = 0
        for item in cryst_list:
            array = np.loadtxt(item, delimiter="," )
            print(array)
            time = array[:, 0]; cryst = array[:, 1]
            color = colors.get(i, 'tab:gray')  # fallback for unspecified iterations
            plt.scatter(time*0.005, cryst, label = "T = %.2f" %temps[i], marker = ".", color = color)
            i = i + 1
        
        plt.title(r"Crystallisation after quench at $\dot{T} = 10^{-3}$")
        plt.xlabel(r"$t/\tau$")
        plt.ylabel(r"$\phi(t)$")
        plt.legend()
        plt.savefig("%s/quench_temps/crystallisation_after_quick_quench.pdf" %(self.path_to_latex_plots_folder))
        plt.show()


    def plot_length_tie_chains(self, savestring = None, mode = "N_tie"):
        for i in range(0, len(self.simulations)):
            simulation = self.simulations[i]
            path_to_tie_chain_file = "%s/tie_chains.txt" %simulation.path_to_home_folder
            tie_chains = pd.read_csv(path_to_tie_chain_file, sep = " ")
            time = simulation.get_simulation_time()
            if mode == "N_tie":
                y = tie_chains["mean_tie_length"].iloc[1:]
                ylabel = r"$N_\text{tie}$"
                savestring = "%s/tie_chains/mean_tie_chain_lengths.pdf" %self.path_to_latex_plots_folder
            elif mode == "f_tie":
                y = tie_chains["f_tie"].iloc[1:]
                ylabel = r"f_\text{tie}"
                savestring = "%s/tie_chains/f_tie.pdf" %self.path_to_latex_plots_folder
            plt.plot(time[1:],y,
                color=self.simulation_colours[simulation], label = "PVA-%i" %simulation.polymer_length)
        plt.xlabel(r"$tc$")
        plt.ylabel(ylabel)
        plt.legend()
        plt.savefig(savestring)
        return 0;


    def plot_n_vs_tie(self, savestring = None):
        for i in range(0, len(self.simulations)):
            simulation = self.simulations[i]
            path_to_tie_chain_file = "%s/tie_chains.txt" %simulation.path_to_home_folder
            tie_chains = pd.read_csv(path_to_tie_chain_file, sep = " ")
            time = simulation.get_simulation_time()
            y = tie_chains["mean_tie_length"].iloc[1:]
            ylabel = r"$N_\text{tie}$"
            savestring = "%s/tie_chains/n_vs_mean_tie_length.pdf" %self.path_to_latex_plots_folder
            plt.scatter(simulation.polymer_length,y[100],
                color=self.simulation_colours[simulation], label = "PVA-%i" %simulation.polymer_length)
        plt.xlabel(r"$N$")
        plt.ylabel(ylabel)
        plt.legend()
        plt.savefig(savestring)
        return 0;



def main():

    PVA_50 = Simulation(50, "../../data/PVA-50/equil", "../data_online/PVA-50/icryst_T088_Tdot_e-3")
    PVA_100 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3")
    PVA_200 = Simulation(200, "../../data/PVA-200/equil", "../data_online/PVA-200/icryst_T088_Tdot_e-3")
    PVA_300 = Simulation(300, "../../data/PVA-300/equil", "../data_online/PVA-300/icryst_T088_Tdot_e-3")
    PVA_500 = Simulation(500, "../../data/PVA-500/equil", "../data_online/PVA-500/icryst_T088_Tdot_e-3")
    PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")


    simulations = [PVA_50, PVA_100, PVA_200, PVA_300, PVA_500, PVA_1000]

    simp = simulation_plots(simulations)

    #simp.plot_monomer_density_and_crossover_values(show_plot=False, mode = "b", marker_size = 10.0)
    #simp.plot_rg_two_polymers_three_times(mode = "re", index_poly_1= 1, index_poly_2= 5)
    #simp.plot_crystallinity()
    #simp.plot_avg_domain_size()
    #simp.plot_crossover_values_vs_chain_length()

    #simp.plot_stem_length()

    #simp.plot_crystallinity_different_quench_temps()
    #simp.plot_length_tie_chains(mode = "N_tie")
    #simp.plot_length_tie_chains(mode = "f_tie")

    #simp.plot_avg_domain_size_and_crossover_values()

    #simp.plot_rg_vs_re()
    simp.plot_n_vs_tie()


if __name__== "__main__":
    main()
