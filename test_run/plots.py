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


    

#plot_volume_monomer_multiple_temps()

plot_crystallinity_multiple_temps()