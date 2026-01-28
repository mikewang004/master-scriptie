from analyse7 import * 
import pandas as pd
from simulation import Simulation
import matplotlib.pyplot as plt

rng = np.random.default_rng(seed = 42)

def convert_input_lattice_to_cryst_array(input_lattice_file, ndot_cutoff = 0.97):
    """Converts a .txt with the lattice sites either 1 or 0 to a df array with columns
    xid yid zid cryst xev yev zev"""

    lattice = np.loadtxt(input_lattice_file)
    lattice_2 = np.loadtxt("debug_cryst_lattice_2.txt")
    print(lattice)

    #Add another layer of zid to lattice 
    zid_zero_layer = np.zeros_like(lattice)
    lattice = np.stack([lattice, lattice_2, zid_zero_layer, zid_zero_layer, zid_zero_layer, zid_zero_layer], axis = 2)

    yid, xid, zid = np.indices(lattice.shape)
    print(lattice.shape)

    xid_flat = xid.flatten()
    yid_flat = yid.flatten()
    zid_flat = zid.flatten()
    cryst_flat = lattice.flatten()
    zeros_help_array = np.zeros(xid_flat.shape)
    zeros_help_array_int = np.zeros(xid_flat.shape, dtype = int)

    cryst_array = pd.DataFrame({
    'xid': xid_flat,
    'yid': yid_flat,
    "zid": zid_flat,
    'cryst': cryst_flat,
    "xev": zeros_help_array,
    "yev": zeros_help_array,
    "zev": zeros_help_array
    })

    #Add RNG for eigenvalues simulation: if cryst = 1, then xev, yev,zev should be at least 0.8 


    for l in range(0, cryst_array.shape[0]):
        if cryst_array.iloc[l, 3] > ndot_cutoff:
            cryst_array.iloc[l, 4:] = rng.uniform(0.9, 1.0, size = 3)



    cryst_array = cryst_array.sort_values(['zid', 'xid', 'yid']).reset_index(drop=True)

    #print(cryst_array)
    cryst_array.to_csv("debug_cryst_analyse_2D.txt", sep = " ", mode = "w")




def plot_label_matrix(label_matrix):
    print(label_matrix)



nridges = 33
#nridges = 6

# convert_input_lattice_to_cryst_array("debug_cryst_lattice.txt")


#last_timestep_long_quench = atom_coords("../../data/pva-100/quick_quench/long_run/equil_t_08_tdot_e-3_time_58800000.txt")
#last_timestep_e5 = atom_coords("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
#last_timestep_e5.get_distribution_eigenvalues(r"Distribution of largest eigenvalues, $T = 0.5, \dot{T} = 10^{-5}$", readfile = "10e5_debug_labdas.txt", savestring = "10e5_T05_labda_dist")
#last_timestep_e5.get_nematic_vector_4(save_ev = True,save_string = "10e5_debug_cryst.txt")
#last_timestep_e5.read_cryst("10e5_debug_cryst.txt")
#last_timestep_e5.read_cryst("debug_cryst_analyse_2D.txt") #Should be independent from actual last_timestep used
#last_timestep_e5.merge_boxes(nridges = nridges)

#last_timestep_e5.get_distribution_eigenvalues(r"Distribution of eigenvalues at $T = 0.5$, $\dot{T} = 10^{-7}$")

#last_timestep_e5.gyration_radius(show_plot= True)
#last_timestep_e5.gyration_radius_debug()



#first_timestep_e5 = atom_coords("../../data/pva-100/cooling_tdot_e-5_time_0.txt")
#first_timestep_e5.bond_bond_correlation()
# #first_timestep_e5.get_nematic_vector_4(save_ev = True,save_string = "10e5_T1_debug_cryst.txt")
# #first_timestep_e5.get_distribution_eigenvalues(r"Distribution of largest eigenvalues, $T = 1.0, \dot{T} = 10^{-5}$", readfile = "10e5_T1_debug_labdas.txt", savestring = "10e5_T1_labda_dist")

# first_timestep_e5.gyration_radius(show_plot = True)
# first_timestep_e5.end_to_end_distance(show_plot = True)



# label_matrix = np.load(hk_label_matrix.npy)






# slow_quench_e5_time_90e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_9000000.txt")
# slow_quench_e5_time_90e5.read_cryst("10e5_T05_timestep_boxes_ev_time_90e5.txt")
# slow_quench_e5_time_90e5.merge_boxes(save = True)

# slow_quench_e5_time_90e5.label_matrix = pd.read_csv("10e5_T05_hk_labels_time_90e5.txt", sep = " ").drop(columns=["Unnamed: 0"])





# hk_label_matrix = slow_quench_e5_time_90e5.label_matrix
# cryst_matrix = slow_quench_e5_time_90e5.df_cryst


# # print(hk_label_matrix[hk_label_matrix.iloc[:, 3]> 0].index)
# # print(cryst_matrix[cryst_matrix.iloc[:, 3] > 0.8].index)


# hk_label_cryst = hk_label_matrix[hk_label_matrix["label"]> 0]
# cryst_matrix_cryst = cryst_matrix[cryst_matrix.iloc[:, 3] > 0.8]

#print(hk_label_cryst.loc[19602])
# print(cryst_matrix_cryst.loc[19602])

# list_with_both = list(set(hk_label_cryst.index) ^ set(cryst_matrix_cryst.index))
# print(list_with_both)


data_path = "../../data/pva-100/quick_quench/long_run"
home_folder_path = "../data_online"

cooling_e3_T088 = Simulation(0.88, -3, "%s%s" %(data_path, "/slurm-e3-T088.out"), "%s/equil_t_088_tdot_e-3_time"%(data_path), 
    home_folder = "../data_online/PVA-100/quick_quench/T088")
cooling_e3_T088.cryst_boxes = pd.read_csv("../data_online/PVA-100/quick_quench/T088/mean_cluster_length_T088_e-3.txt", sep = " ").drop(columns=["Unnamed: 0"])
print(cooling_e3_T088.cryst_boxes)

#plt.scatter(cooling_e3_T088.cryst_boxes.iloc[:, 0], cooling_e3_T088.cryst_boxes.iloc[:, 2], label = "clusters w/ >= 2 members")
plt.scatter(cooling_e3_T088.cryst_boxes.iloc[:, 0], cooling_e3_T088.cryst_boxes.iloc[:, 3], label = "independent clusters")
plt.scatter(cooling_e3_T088.cryst_boxes.iloc[:, 0], cooling_e3_T088.cryst_boxes.iloc[:, 4], label = "mean size cryst domains")
#plt.scatter(cooling_e3_T088.cryst_boxes.iloc[:, 0], cooling_e3_T088.cryst_boxes.iloc[:, 5], label = "crystalline grid elements")

plt.title(r"Crystallisation properties, PVA-100, T = 0.88, $\dot{T} = 10^{-3}$")
plt.xlabel("time")
plt.legend()
plt.savefig("T088_e-3_independent_clusters_mean_size_cryst_domains.pdf")
plt.show()