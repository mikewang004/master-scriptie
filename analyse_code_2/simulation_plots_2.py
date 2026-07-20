from analyse9 import polymer, atom_coords
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

    def __init__(self, simulations, caption_font = 11, legend_font = 9, title_font = 18):
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
            idx, xkn, ykn, popt = simulation.domain_analysis.find_knee()
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
        # Calc max density (is at PVA-50)
        pva_50 = self.simulations[0]
        max_density = pva_50.get_polymer_by_time(0).atom_coords.n_atoms/pva_50.df_slurm_sim_data["Volume"].iloc[-1]

        for simulation in self.simulations:
            poly = simulation.get_polymer_by_time(0)
            time = simulation.df_slurm_sim_data["Step"] * simulation.timestep
            monomer_density = (poly.atom_coords.n_atoms/simulation.df_slurm_sim_data["Volume"])
            # Top subplot
            ax1.scatter(time, monomer_density, label=f"PVA-{simulation.polymer_length}", color = self.simulation_colours[simulation])

            idx, xkn, ykn, popt = simulation.domain_analysis.find_knee()
            polymer_lengths.append(simulation.polymer_length)
            times.append(xkn); crossovers.append(ykn); idx_list.append(idx)
            time_con = np.linspace(0, time.max(), 50000)
            ax2.scatter(time, monomer_density, label=f"PVA-{simulation.polymer_length}", color = self.simulation_colours[simulation], marker = ".")
            ax2.plot(time_con, fit_functions.double_exp(time_con, *popt), color = self.simulation_colours[simulation], linewidth = 3.0)

        tc = pd.DataFrame(
            {
                "index": idx_list,
                "time": times,
                "monomer density": crossovers,
            },
            index=polymer_lengths,
        )
        tc.index.name = "polymer_lengths"

        for i in range(0,len(times)):
            xkn = times[i]; ykn = crossovers[i]
            ax2.scatter(xkn, ykn, marker="x", color="black", zorder = 3)

        print(tc.head())
        tc.to_csv("../data_online/crossover_times.txt", sep=" ")


        # Horizontal lines in both subplots

        # Labels & legend
        ax2.set_xlabel(r"$t/\tau$", fontsize = self.caption_font)
        ax1.set_ylabel(r"$n_\text{monomers}/\sigma^3$", fontsize = self.caption_font)
        ax2.set_ylabel(r"$n_\text{monomers}/\sigma^3$", fontsize = self.caption_font)

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
        fig.savefig("%s/crossover_point/crossover_density_vs_time_different_chains_subplots.pdf" %self.path_to_latex_plots_folder)
        if show_plot == True:
            plt.show()


    def plot_crystallinity(self, savestring = "plots/crystallinity_vs_time.pdf", show_plot = True):
        plt.figure(figsize = (self.std_width*1.25, self.std_height*1.5))
        for i in range(0, len(self.simulations)):
            simulation = self.simulations[i]
            current_domain_file = simulation.domain_analysis.read_avg_domain_size()
            time = simulation.get_simulation_time()
            

            
            plt.scatter(time[:len(current_domain_file["crystallinity"])], current_domain_file["crystallinity"].iloc[:len(time)], 
                label = "PVA-%i" %(simulation.polymer_length), c= self.simulation_colours[simulation], marker = ".")


        plt.legend(fontsize=self.caption_font)
        plt.xlabel(r"$t/\tau_c$", fontsize = self.caption_font)
        plt.ylabel(r"$\phi$", fontsize = self.caption_font)
        #plt.xscale("log")
        #plt.title("Mean domain size, various chains")
        if savestring is not None:
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
        plt.xscale("log")
        #plt.title("Mean domain size, various chains")
        if savestring is not None:
            plt.savefig("%s/%s" %(self.save_folder, savestring))
        if show_plot == True:
            plt.show()

    


    def plot_rg_two_polymers_three_times(self, mode = "rg", savestring_default = True, show_plot = True):
        """Mode can either be 'rg' or "re". 
        Note: PVA-100 has 174 items, PVA-1000 119."""
        fig, axes = plt.subplots(
            2, 1,
            figsize=(1.5 * self.std_width, 2 * self.std_height),
            sharex=True
        )

        PVA_100 = self.simulations[1]; PVA_1000 = self.simulations[-1]

        index_t_start = 0; index_t_end_1000 = 119; index_t_end_100 = 140
        times_PVA_100 = [index_t_start, PVA_100.tc_idx, index_t_end_100]
        times_PVA_1000 = [index_t_start, PVA_1000.tc_idx, index_t_end_1000]
        times_different_PVA = [times_PVA_100, times_PVA_1000]
        polymer_list = [PVA_100, PVA_1000]
        letter_subplot_list = ["(a)", "(b)", "(c)"]
        #print(times_PVA_100[i])
        #current_poly_100 = PVA_100.get_polymer_by_time(int(times_PVA_100[i]*1200000))
        #current_poly_1000 = PVA_1000.get_polymer_by_time(int(times_PVA_1000[i]*1200000))
        for i in range(0, len(times_different_PVA)):
            polymer_times = times_different_PVA[i]
            for j in range(0, len(polymer_times)): 
                current_time = int(polymer_times[j])*1200000
                current_poly = polymer_list[i].get_polymer_by_time(current_time)
                if mode == "rg":
                    current_poly.gyration_radius()
                    # values, bins, __ = axes[i].hist(current_poly.results.gyration_radius_distribution/
                    #     np.sqrt(current_poly.results.mean_gyration_radius), bins=50, 
                    #     color=self.times_colours["%i" %(2*j)], density = True, histtype = "step", 
                    #     label = r"$%i t/t_c$" %(int(current_poly.atom_coords.current_timestep/polymer_list[i].tc_time)))
                    counts, bin_edges = np.histogram(current_poly.results.gyration_radius_distribution/
                         np.sqrt(current_poly.results.mean_gyration_radius), bins = 80, density = True)

                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) /2

                    smooth_counts = sp.ndimage.gaussian_filter1d(counts, sigma = 2.0)
                    axes[i].plot(bin_centers, smooth_counts, color=self.times_colours["%i" %(2*j)], linestyle = "-", marker = ".",
                    label = r"$%i t/t_c$" %(int(current_poly.atom_coords.current_timestep/polymer_list[i].tc_time)))
                    axes[i].set_xlabel(r"$R_g/ \sqrt{\langle R_g^2 \rangle}$", fontsize = self.caption_font)
                    savestring = "%s/polymer_conformation/rg_pva_100_1000.pdf" %self.path_to_latex_plots_folder

                elif mode == "re":
                    current_poly.end_to_end_distance()
                    counts, bin_edges = np.histogram(current_poly.results.end_to_end_distribution/current_poly.results.mean_squared_end_to_end, bins = 100, density= True)
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) /2
                    smooth_counts = sp.ndimage.gaussian_filter1d(counts, sigma = 2.0)
                    axes[i].plot(bin_centers, smooth_counts, color=self.times_colours["%i" %(2*j)], linestyle = "-", marker = ".",
                    label = r"$%i t/t_c$" %(int(current_poly.atom_coords.current_timestep/polymer_list[i].tc_time)))
                #values, bins, __ = axes[i].hist(current_poly_1000.results.end_to_end_distribution/
                #    (current_poly_1000.results.mean_squared_end_to_end), bins=50, color=self.simulation_colours[PVA_1000], density = True, histtype = "step", label = "PVA-%i" %PVA_1000.polymer_length)
                    axes[i].set_xlabel(r"$R_e/ \sqrt{\langle R_e^2 \rangle}$", fontsize = self.caption_font)
                    savestring = "%s/polymer_conformation/re_pva_100_1000.pdf" %self.path_to_latex_plots_folder
                elif mode == "nematic":
                    current_poly.nematic_distributuion()
                    counts, bin_edges = np.histogram(current_poly.results.nematic_value_dist, bins = 100, density = True)
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) /2
                    smooth_counts = sp.ndimage.gaussian_filter1d(counts, sigma = 1.0)
                    axes[i].plot(bin_centers, smooth_counts, color=self.times_colours["%i" %(2*j)], linestyle = "-", marker = "o",
                    label = r"$%i t/t_c$" %(int(current_poly.atom_coords.current_timestep/polymer_list[i].tc_time)))
                # values, bins, __ = axes[i].hist(current_poly_1000.results.nematic_value_dist, 
                #     bins=100, color=self.simulation_colours[PVA_1000], density = True, histtype = "step", label = "PVA-%i" %PVA_1000.polymer_length)
                    axes[i].set_xlabel(r"$\lambda$", fontsize = self.caption_font)
                    
                    savestring = "%s/polymer_conformation/nem_value_pva_100_1000.pdf" %self.path_to_latex_plots_folder


                elif mode == "cryst_domain_dist":
                    cryst_domain_dist = pd.read_csv("%s/domain_distribution/%s_domain_dist_%s.txt" 
                        %(polymer_list[i].path_to_home_folder, polymer_list[i].lammps_dump_prefix, current_time), sep = " ")
                    #print(cryst_domain_dist["volume_pdf"])

                    axes[i].plot(cryst_domain_dist["volume"], cryst_domain_dist["volume_pdf"], 
                        label = r"$%i t/t_c$" %(int(current_poly.atom_coords.current_timestep/polymer_list[i].tc_time)))
                    savestring = "%s/domain_analysis/nem_value_pva_100_1000.pdf" %self.path_to_latex_plots_folder
                    axes[i].set_xscale("log")
                    axes[i].set_yscale("log")
                    axes[i].set_xlabel(r"$V_\text{domain}$")
                    axes[i].set_ylabel(r"$d\phi/dV$")

            axes[i].tick_params(labelbottom=True, labelleft=True)
            axes[i].text(0.02, 0.95, letter_subplot_list[i],
                transform=axes[i].transAxes,
                fontsize=plt_caption_font,
                va="top", ha="left")
            if mode == "nematic":
                axes[i].vlines(PVA_100.cryst_cutoff, 0, np.max(counts), color = "red",linestyles = "dotted", label = r"crystallisation cutoff")
            axes[i].legend(fontsize = self.caption_font)



        axes[0].set_title("PVA-100")
        axes[1].set_title("PVA-1000")

        axes[0].legend(fontsize=self.caption_font)
        axes[1].legend(fontsize=self.caption_font)
        plt.tight_layout()
        if savestring_default == True:
            plt.savefig(savestring)
        if show_plot == True:
            plt.show()

    def plot_stem_length(self, savestring = "plots/stem_length_pva_1000.pdf", show_plot = True):
        plt.figure(figsize = (self.std_width*1.25, self.std_height*1.5))

        simulation = self.simulations[-1]
        print(simulation.df_slurm_sim_data)
        positions = []
        min_stem_lengths = []
        max_time_index = 50
        for i in range(0, max_time_index):
            current_poly = simulation.get_polymer_by_time(simulation.df_slurm_sim_data["Step"].iloc[i])
            position, min_stem_length = current_poly.bond_bond_corr_min_value()

            positions.append(position)
            min_stem_lengths.append(min_stem_length)

        time = simulation.get_simulation_time()
        plt.scatter(time[1:max_time_index], positions[1:], color=self.simulation_colours[simulation], label = "PVA-%i" %simulation.polymer_length)
        plt.ylim((0, 50))
        plt.legend()
        plt.ylabel("Stem length")
        plt.xlabel("$t/t_c$")
        plt.savefig(savestring)
        plt.show()




def main():

    PVA_50 = Simulation(50, "../../data/PVA-50/equil", "../data_online/PVA-50/icryst_T088_Tdot_e-3")
    PVA_100 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3")
    PVA_200 = Simulation(200, "../../data/PVA-200/equil", "../data_online/PVA-200/icryst_T088_Tdot_e-3")
    PVA_300 = Simulation(300, "../../data/PVA-300/equil", "../data_online/PVA-300/icryst_T088_Tdot_e-3")
    PVA_500 = Simulation(500, "../../data/PVA-500/equil", "../data_online/PVA-500/icryst_T088_Tdot_e-3")
    PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")


    simulations = [PVA_50, PVA_100, PVA_200, PVA_300, PVA_500, PVA_1000]

    simp = simulation_plots(simulations)

    #simp.plot_monomer_density_and_crossover_values(show_plot=True)
    simp.plot_rg_two_polymers_three_times(mode = "cryst_domain_dist")
    #simp.plot_crystallinity_avrami(savestring= None)
    #simp.plot_avg_domain_size()
    #simp.plot_crossover_values_vs_chain_length()

    #simp.plot_stem_length()


if __name__== "__main__":
    main()
