from analyse6 import * 


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








nridges = 33

# convert_input_lattice_to_cryst_array("debug_cryst_lattice.txt")


#last_timestep_long_quench = atom_coords("../../data/pva-100/quick_quench/long_run/equil_t_08_tdot_e-3_time_58800000.txt")
last_timestep_e5 = atom_coords("../../data/pva-100/cooling_tdot_e-5_time_0.txt")
last_timestep_e5.read_cryst("10e5_debug_cryst.txt")
#last_timestep_e5.read_cryst("debug_cryst_analyse_2D.txt") #Should be independent from actual last_timestep used
last_timestep_e5.merge_boxes(nridges = nridges)

#last_timestep_e5.gyration_tensor(show_plot= True)
#last_timestep_e5.gyration_radius_debug()


