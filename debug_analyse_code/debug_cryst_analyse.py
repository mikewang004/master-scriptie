from analyse7 import * 


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




def plot_label_matrix(lattice):
    occupied = lattice > 0
    x, y, z = np.where(occupied)
    values = lattice[occupied]  # could be used for coloring

    # --- 3. Plot as 3D scatter on cubic lattice ---
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Color by value; remove "c=values" if you just want uniform color
    sc = ax.scatter(x, y, z,
                    c=values,
                    cmap="viridis",
                    marker="s",  # square markers look more like voxels
                    s=20)

    # Make axes look cubic
    ax.set_box_aspect([1, 1, 1])  # equal aspect ratio

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Site value")

    plt.tight_layout()
    plt.show()

def plot_only_largest_clusters(lattice, n_keep):
    labels = lattice

    # ---------- 3. Count size of each label ----------
    sizes = np.bincount(labels.ravel())          # sizes[i] = voxels with label i
    if len(sizes) > 0:
        sizes[0] = 0                             # label 0 = background, EXCLUDE

    # List of non-zero labels that actually appear
    nonzero_labels = np.nonzero(sizes)[0]        # labels with size > 0
    n_clusters = len(nonzero_labels)

    # Adjust n_keep if user asks for more than exist
    n_keep = min(n_keep, n_clusters)

    # Indices of the n_keep largest *non-zero* labels
    largest_labels = np.argsort(sizes)[-n_keep:]  # these are != 0 because sizes[0]=0
    print("Largest labels (excluding 0):", largest_labels)
    print("Their sizes:", sizes[largest_labels])

    # ---------- 4. Mask: only voxels in those clusters (still excludes label 0) ----------
    mask = np.isin(labels, largest_labels)   # True only for selected cluster labels
    x, y, z = np.where(mask)
    cluster_ids = labels[mask]               # used for coloring

    # ---------- 5. Plot ONLY the selected clusters (no label 0) ----------
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(
        x, y, z,
        c=cluster_ids,
        cmap="tab20",     # each cluster a different color
        marker="s",
        s=30
    )



    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Cluster label")

    plt.tight_layout()
    plt.show()



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



label_matrix = np.load("hk_label_matrix.npy")
plot_label_matrix(label_matrix)
plot_only_largest_clusters(label_matrix, 10)

