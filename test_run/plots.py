from analyse5 import *

def plot_volume_monomer_multiple_temps():
    list_atom_coords_cooling_tdot_e3 = get_list_atom_coords("../../data/pva-100/cooling_tdot_e-3_time", 21, endtime= 1e5)
    list_atom_coords_cooling_tdot_e4 = get_list_atom_coords("../../data/pva-100/cooling_tdot_e-4_time", 21, endtime= 1e6)
    list_atom_coords_cooling_tdot_e5 = get_list_atom_coords("../../data/pva-100/cooling_tdot_e-5_time", 21, endtime= 1e7)

    master_list = [list_atom_coords_cooling_tdot_e3,list_atom_coords_cooling_tdot_e4,list_atom_coords_cooling_tdot_e5]
    tdots = [3,4,5]

    #Plot volume per monomer for all plots in one line 
    starttemp = 1.0; endtemp = 0.5;
    temps = np.linspace(starttemp, endtemp, num = 21)
    i = 0
    for atom_coords_tdots in master_list:
        current_tdot = tdots[i]
        list_volumes = []
        for atom_coords in atom_coords_tdots:
            list_volumes.append((atom_coords.volume)/atom_coords.n_atoms)
            print(atom_coords.volume/atom_coords.n_atoms)
        plt.scatter(temps, list_volumes, label = r"$\dot{T} = 10^{-%s}$" %(current_tdot))
        i = i + 1

    plt.title("Volume per monomer vs temperature, different cooling rates")
    plt.xlabel("Temperature [unitless]")
    plt.ylabel("Volume per monomer")
    plt.legend()
    plt.savefig("volume_monomer_different_tdots.pdf")
    plt.show()


def plot_crystallinity_multiple_temps():
    cryst_e3 = np.loadtxt("test_wholebox_frac_cryst_heating_100_tmin_0.5_ttime_10e5.txt")
    cryst_e4 = np.loadtxt("test_wholebox_frac_cryst_heating_100_tmin_0.5_ttime_10e6.txt")
    cryst_e5 = np.loadtxt("test_wholebox_frac_cryst_heating_100_tmin_0.5_ttime_10e7.txt")
    plt.scatter(cryst_e3[:, 0], cryst_e3[:, 1], label = r"$\dot{T} = 10^{-3}$")
    plt.scatter(cryst_e4[:, 0], cryst_e4[:, 1], label = r"$\dot{T} = 10^{-4}$")
    plt.scatter(cryst_e5[:, 0], cryst_e5[:, 1], label = r"$\dot{T} = 10^{-5}$")

    plt.legend()
    plt.title("Crystallinity vs temperature, cooling process, PVA-100")
    plt.xlabel("Temperature [unitless]")
    plt.ylabel("Fraction of crystallinity")

    plt.savefig("crystallinity_multiple_temps.pdf")
    plt.show()


def quick_quench_volume_monomer():
    endtime = 14400000
    no_steps = 13
    t08 = get_list_atom_coords("../../data/pva-100/quick_quench/equil_t_08_tdot_e-3_time", no_steps, endtime= endtime)
    t07 = get_list_atom_coords("../../data/pva-100/quick_quench/equil_t_07_tdot_e-3_time", no_steps, endtime= endtime)
    #list_atom_coords_cooling_tdot_e5 = get_list_atom_coords("../../data/pva-100/cooling_tdot_e-5_time", 21, endtime= 1e7)

    master_list = [t08,t07]
    tdots = [0.8, 0.7]

    #Plot volume per monomer for all plots in one line 
    starttemp = 1.0; endtemp = 0.5;
    time = np.linspace(0, endtime, no_steps)
    #temps = np.linspace(starttemp, endtemp, num = 21)
    i = 0
    for atom_coords_tdots in master_list:
        current_tdot = tdots[i]
        list_volumes = []
        for atom_coords in atom_coords_tdots:
            list_volumes.append((atom_coords.volume)/atom_coords.n_atoms)
            print(atom_coords.volume/atom_coords.n_atoms)
        plt.scatter(time, list_volumes, label = r"${T} = %s$" %(tdots[i]))
        i = i + 1

    plt.title("Volume per monomer after quick quench, time evoluation")
    plt.xlabel("Temperature [unitless]")
    plt.ylabel("Volume per monomer")
    plt.legend()
    plt.savefig("quick_quench_monomer_temp.pdf")
    plt.show()


    

#plot_volume_monomer_multiple_temps()

#plot_crystallinity_multiple_temps()

#quick_quench_volume_monomer()


#const_t08_e3_cooling = get_list_atom_coords("../../data/pva-100/quick_quench/equil_t_08_tdot_e-3_time", 18, endtime= 20400000) #Kept constant at T=0.8 after quench from T = 1.0
#plot_order_param(const_t08_e3_cooling, "Crystallinity vs time, equilibrated system", savestring="equil_t_08_24e7", plot_time_instead= True, starttime = 0, endtime = 20400000, n_samples = 18)


last_timestep_e5 = atom_coords("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
last_timestep_e5.gyration_radius()

# last_timestep_e4 = atom_coords("../../data/pva-100/cooling_tdot_e-4_time_1000000.txt")
# last_timestep_e4.gyration_radius()

#first_timestep_e5 = atom_coords("../../data/pva-100/cooling_tdot_e-5_time_0.txt")
#first_timestep_e5.gyration_radius()
#first_timestep_e5.end_to_end_distance(show_plot= True)