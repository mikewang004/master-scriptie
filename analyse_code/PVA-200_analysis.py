from analyse7 import polymer, atom_coords, plot_hk_matrix_2d
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from simulation import Simulation#, avrami_fit, avrami_eq

ndot_cutoff = 0.97 #Threshold above which the crystalline domains can be merged 
cryst_cutoff = 0.8 #Threshold for a cell to be considered crystalline





def pva_200_analysis(data_path):
    icryst_PVA_200_T088 = Simulation(0.88, -3, "%s/slurm-PVA-200_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-200_equil_t_088_tdot_e-3_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-200/icryst_T088_Tdot_e-3", polymer_length=200, home_folder_override= True)
    icryst_PVA_200_T088.calc_crystallisation()
    icryst_PVA_200_T088.calc_avg_domain_size()
    return icryst_PVA_200_T088

def pva_50_analysis(data_path):
    icryst_PVA_50_T088 = Simulation(0.88, -3, "%s/slurm-PVA-50_equil_t_088_tdot_e-3.out" %(data_path), "%s/PVA-50_equil_t_088_tdot_e-3_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-50/icryst_T088_Tdot_e-3",polymer_length=50, home_folder_override= True)
    icryst_PVA_50_T088.calc_crystallisation()
    icryst_PVA_50_T088.calc_avg_domain_size()
    return icryst_PVA_50_T088

def pva_100_analysis(data_path):
    icryst_PVA_100_T088 = Simulation(0.88, -3, "%s%s" %(data_path, "/slurm-e3-T088.out"), "%s/equil_t_088_tdot_e-3_time"%(data_path), no_runs = 1, 
        home_folder= "../data_online/PVA-100/quick_quench")
    icryst_PVA_100_T088.calc_crystallisation()
    icryst_PVA_100_T088.calc_avg_domain_size(ndot_cutoff= 0.98)
    return icryst_PVA_100_T088



def pva_300_analysis(data_path):
    icryst_PVA_300_T088 = Simulation(0.88, -3, "%s/slurm-PVA-300_equil_t_088_tdot_e-3_sim1.out" %(data_path), "%s/PVA-300_equil_t_088_tdot_e-3_sim1_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-300/icryst_T088_Tdot_e-3/sim1", polymer_length=300, home_folder_override= True)
    icryst_PVA_300_T088.calc_crystallisation()
    icryst_PVA_300_T088.calc_avg_domain_size()
    return icryst_PVA_300_T088

def pva_500_analysis(data_path):
    icryst_PVA_500_T088 = Simulation(0.88, -3, "%s/slurm-PVA-500_equil_t_088_tdot_e-3_sim1.out" %(data_path), "%s/PVA-500_equil_t_088_tdot_e-3_sim1_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-500/icryst_T088_Tdot_e-3/sim1", polymer_length=500, home_folder_override= True)
    icryst_PVA_500_T088.calc_crystallisation()
    icryst_PVA_500_T088.calc_avg_domain_size()
    return icryst_PVA_500_T088

def pva_1000_analysis(data_path):
    icryst_PVA_1000_T088 = Simulation(0.88, -3, "%s/slurm-PVA-1000_equil_t_088_tdot_e-3_sim1.out" %(data_path), "%s/PVA-1000_equil_t_088_tdot_e-3_sim1_time" %(data_path), no_runs=1,
        home_folder="../data_online/PVA-1000/icryst_T088_Tdot_e-3/sim1", polymer_length=1000, home_folder_override= True)
    icryst_PVA_1000_T088.calc_crystallisation()
    icryst_PVA_1000_T088.calc_avg_domain_size()
    return icryst_PVA_1000_T088

def get_domain_distribution_polymer(polymer):
    """Returns distribution of domain sizes"""
    label_matrix = polymer.merge_boxes(print_results= True, ndot_cutoff=0.98)
    print(label_matrix)
    #Get distribution by counting label occurances
    label_matrix = label_matrix[label_matrix != 0]
    #print(label_matrix[label_matrix != 0])
    bins = [1,2,3,4,5,7,10,20,30,40,50,60,70,100]
    n, _, _ = plt.hist(label_matrix.flatten(), bins = bins)
    print(n, bins)
    plt.title(r"Crystalline domain size distribution, PVA-100, time %i, $\phi = %.3f$" %(polymer.timestep, polymer.frac_cryst))
    plt.savefig("cryst_domain_distro_pva_100_time_%i.pdf" %polymer.timestep)
    plt.show()
    return 0;


def plot_crystallisation_different_polymer_lengths(simulation_list, save: bool = False, savestring = None, labels_list: list = None, plot_equal_length: bool = False):

    length_list =[]

    for simulation in simulation_list:
        length_list.append(simulation.cryst.shape[0])
        if plot_equal_length == True:
            length = min(length_list)
        else:
            length = max(length_list)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # First plot mean cluster length 
    for simulation in simulation_list:
        # ax1.plot(simulation.mean_cluster_length.iloc[:length, 0], simulation.mean_cluster_length.iloc[:length, 3], 
        #     label = r"independent clusters, PVA-%i" %(simulation.polymer_length))
        label_cryst = r"crystallinity,  PVA-%i" %(simulation.polymer_length)
        ax1.plot(simulation.cryst[:length, 0], simulation.cryst[:length, 1], 
            label = label_cryst)
        ax2.scatter(simulation.mean_cluster_length.iloc[:length, 0], (simulation.mean_cluster_length.iloc[:length, 4])**(1/3) *2,
            label =  r"mean domain size, PVA-%i" %(simulation.polymer_length))
        temp = simulation.temp

    ax1.set_xlabel("time")
    #ax1.set_ylabel(r"amount of independent clusters")
    ax1.set_ylabel(r"$\phi$")
    ax2.set_ylabel(r"length scale ($\sigma$)")
    #fig.tight_layout()
    fig.legend(loc = "lower right", bbox_to_anchor=(0.895, 0.115))
    #fig.suptitle(r"Independent clusters and mean domain size, $\dot{T} = 10^{%i}$" %(simulation.cooling_rate))
    fig.suptitle(r"Isothermal crystallisation, T=%.2f"%(temp))
    fig.savefig("plots/e%i_no_clusters_mean_domain_size.pdf"%(simulation.cooling_rate))
    plt.show()

def avrami_eq(t, n, b, a):
    #return np.log(1-t)
    return a * (t**n) 
    #return np.exp(-a*t**n) + b


def avrami_fit(t, y):
    popt, pcov = sp.optimize.curve_fit(avrami_eq, t, y, maxfev = 1000000)
    return popt, pcov

def plot_avrami(simulation_list: list, save: bool = False, savestring = None, show_plot = False):

    #plt.xscale("log")
    #plt.yscale("log")
    #plt.ylim(0.1, 0.6)
    #plt.xlim(10**7, 10**8)
    for simulation in simulation_list:
        print(simulation.cryst.shape)
        used_simulation_cryst = simulation.cryst[:50, :]# / 12000000
        used_simulation_cryst[:, 0] /= 10000000
        #phi0 = np.log(1 - (used_simulation_cryst[:, 1]))
        phi0 = used_simulation_cryst[:, 1]
        #print(used_simulation_cryst)
        #print(used_simulation_cryst)
        #print(np.log(1 - (used_simulation_cryst[:, 1])))
        popt, pcov = avrami_fit(used_simulation_cryst[:, 0], phi0)
        print(popt)
        sim_times = np.linspace(np.min(used_simulation_cryst[:, 0]), np.max(used_simulation_cryst[:, 0]), 10000)
        #plt.xscale("log")        # make y-axis logarithmic
        plt.plot(sim_times, avrami_eq(sim_times, *popt),
             label=f"n = {popt[0]:.4f}, b = {popt[1]:.4f}, a = {popt[2]:.10f}")
        plt.scatter(used_simulation_cryst[:, 0], phi0, label = "experimental data")
        #plt.xscale("log")
        #plt.yscale("log")

    plt.title(r"$\phi(t)$, PVA-100, T=%.2f, $\dot{T}=10^{%i}$" %(simulation.temp, simulation.cooling_rate))
    plt.xlabel("time")
    plt.ylabel("fraction of crystallinity")

    plt.legend()
    if save is True:
        if savestring is None:
            savestring = "plots/avrami_10e%i_T%s_crystallinity_early.pdf" %(simulation.cooling_rate, str(simulation.temp).replace(".", ""))
        plt.savefig(savestring)
    if show_plot is True:
        plt.show()
    plt.close()



def log_binning(imax):
    """
    Create bin edges corresponding to the C++ binning scheme
    for integer cluster sizes 1..imax.

    Bins:
        1,2,3,4,5
        6-7
        8-10
        11-20, 21-30, ..., 61-70
        71-100
        101-200, 201-300, ..., 601-700
        701-1000
        1001-2000, 2001-3000, ... up to imax

    Uses [edge[k], edge[k+1]) with 0.5 offsets so integers fall correctly.
    """
    edges = [0.5]  # so size 1 is in [0.5, 1.5)

    # 1–5 : one bin per integer
    for i in range(1, 6):
        if i > imax:
            edges.append(imax + 0.5)
            return np.array(edges)
        edges.append(i + 0.5)

    # 6–7
    if imax >= 6:
        edges.append(min(7, imax) + 0.5)
        if imax <= 7:
            return np.array(edges)
    else:
        edges.append(imax + 0.5)
        return np.array(edges)

    # 8–10
    if imax >= 8:
        edges.append(min(10, imax) + 0.5)
        if imax <= 10:
            return np.array(edges)
    else:
        edges.append(imax + 0.5)
        return np.array(edges)

    # 11–20, 21–30, ..., 61–70 (width 10)
    for upper in range(20, 71, 10):
        lower = upper - 9
        if imax > lower:
            edges.append(min(upper, imax) + 0.5)
            if imax <= upper:
                return np.array(edges)
        else:
            edges.append(imax + 0.5)
            return np.array(edges)

    # 71–100
    if imax > 70:
        edges.append(min(100, imax) + 0.5)
        if imax <= 100:
            return np.array(edges)
    else:
        edges.append(imax + 0.5)
        return np.array(edges)

    # 101–200, 201–300, ..., 601–700 (width 100)
    for upper in range(200, 701, 100):
        lower = upper - 99
        if imax > lower:
            edges.append(min(upper, imax) + 0.5)
            if imax <= upper:
                return np.array(edges)
        else:
            edges.append(imax + 0.5)
            return np.array(edges)

    # 701–1000
    if imax > 700:
        edges.append(min(1000, imax) + 0.5)
        if imax <= 1000:
            return np.array(edges)
    else:
        edges.append(imax + 0.5)
        return np.array(edges)

    # 1001–2000, 2001–3000, ... up to imax (step 1000)
    upper = 2000
    while upper < imax:
        edges.append(min(upper, imax) + 0.5)
        if imax <= upper:
            return np.array(edges)
        upper += 1000

    if edges[-1] < imax + 0.5:
        edges.append(imax + 0.5)

    return np.array(edges)

def cpp_style_cluster_hist(polymer, plot=True):
    """
    Parameters
    ----------
    label_array : np.ndarray
        3D array (e.g. 33x33x33) with integer labels in [0..n],
        where 0 is background and 1..n are cluster labels.
    plot : bool
        If True, plot a histogram using Matplotlib.

    Returns
    -------
    bin_centers : np.ndarray
        Center of each bin (approximate cluster size).
    bin_counts : np.ndarray
        Number of clusters in each bin.
    bin_edges : np.ndarray
        Edges used for the histogram.
    """

    # ---- 1. Compute cluster sizes ----
    label_array = polymer.merge_boxes(print_results= True, ndot_cutoff=0.98)
    flat = label_array.ravel()
    max_label = flat.max()

    # counts_per_label[i] = number of voxels with label i
    counts_per_label = np.bincount(flat, minlength=max_label + 1)

    # Assume label 0 is background
    cluster_sizes = counts_per_label[1:]  # sizes for labels 1..max_label
    # Remove any empty labels (should be rare but safe)
    cluster_sizes = cluster_sizes[cluster_sizes > 0]

    if cluster_sizes.size == 0:
        raise ValueError("No non-zero cluster labels found in the array.")

    imax = int(cluster_sizes.max())

    # ---- 2. Build C++-style bin edges ----


    bin_edges = log_binning(imax)

    # ---- 3. Histogram of cluster sizes ----
    bin_counts, edges_used = np.histogram(cluster_sizes, bins=bin_edges)

    bin_centers = 0.5 * (edges_used[:-1] + edges_used[1:])
    bin_widths = np.diff(edges_used)



    # ---- 4. Optional plotting ----
    if plot:
        #plt.figure(figsize=(6, 4))
        plt.bar(bin_centers, bin_counts,
                width=bin_widths, align='center', edgecolor='k')
        plt.xlabel("Cluster size")
        plt.ylabel("Number of clusters in bin")
        plt.xscale("log")  # optional but typically useful here
        #plt.tight_layout()
        print(polymer.timestep)
        plt.title(r"Crystalline domain size distribution, PVA-%i, $\phi = %.3f$" %(polymer.atom_coords.polymer_length, polymer.frac_cryst))
        plt.savefig("cryst_domain_distro_pva_200_quench_time_%i.pdf" %polymer.timestep)
        plt.show()

    return bin_centers, bin_counts, edges_used


def main():
    data_path_50 = "../../data/PVA-50"
    data_path_200 = "../../data/PVA-200"
    data_path_100 = "../../data/pva-100/quick_quench/long_run"
    data_path_300 = "../../data/PVA-300"
    data_path_500 = "../../data/PVA-500"
    data_path_1000 = "../../data/PVA-1000"

    # cooling_e4_T085 = Simulation(0.85, -4, "%s%s" %(data_path_100, "/slurm-e4-T085.out"), "%s/equil_t_085_tdot_e-4_time"%(data_path_100), no_runs = 1, 
    #     home_folder= "../data_online/PVA-100/quick_quench")
    # cooling_e4_T085.calc_crystallisation()
    # cooling_e4_T085.calc_avg_domain_size()
    #pva_50_analysis(data_path_50)
    icryst_PVA_200_T088 = pva_200_analysis(data_path_200)
    icryst_PVA_50_T088 = pva_50_analysis(data_path_50)
    icryst_PVA_100_T085 = pva_100_analysis(data_path_100)
    icryst_PVA_300_T088 = pva_300_analysis(data_path_300)
    icryst_PVA_500_T088 = pva_500_analysis(data_path_500)
    icryst_PVA_1000_T088 = pva_1000_analysis(data_path_1000)

    # quench_e5 = Simulation(0.5, -5, "../../data/pva-100/genua_cooling_100_ttime_10e7.out", "../../data/pva-100/genua_cooling_100_tmin_0.5_ttime_10e7",
    #     no_runs=1, home_folder = "../data_online/PVA-100/quench/cryst_quench_e-5", home_folder_override= True)
    # quench_e5.calc_crystallisation()
    # quench_e5.calc_avg_domain_size()
    current_polymer = icryst_PVA_100_T085.get_polymer_by_count(20)
    #current_polymer.merge_boxes()
    #current_polymer.merge_boxes_2(ndot_cutoff= 0.97, print_results= True)
    #print(current_polymer.atom_distribution())
    #cpp_style_cluster_hist(current_polymer)
    #print(current_polymer.bond_distribution())
    #current_polymer.atom_coords.get_nematic_vector_5(save_string= "test.txt")
    #get_domain_distribution_polymer(current_polymer)
    plot_hk_matrix_2d(current_polymer, ndot_cutoff= 0.97)
    #icryst_PVA_100_T088.get_polymer_by_count(99).atom_coords.get_nematic_vector_5()

    #polymer_list = [icryst_PVA_50_T088, icryst_PVA_200_T088, icryst_PVA_300_T088, icryst_PVA_500_T088, icryst_PVA_1000_T088]
    #polymer_list = [icryst_PVA_100_T085]
    #plot_crystallisation_different_polymer_lengths(polymer_list, plot_equal_length= False, save= False)
    #plot_avrami([icryst_PVA_100_T085], show_plot= True, save = False)



if __name__ == "__main__":
    main()

