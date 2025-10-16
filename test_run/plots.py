from analyse5 import *


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
    plt.scatter(temps, list_volumes, label = r"tdot $= 10e^{-%s}$" %(current_tdot))
    i = i + 1

plt.title("Volume per monomer vs temperature, different cooling rates")
plt.xlabel("Temperature [unitless]")
plt.ylabel("Volume per monomer")
plt.legend()
plt.savefig("volume_monomer_different_tdots.pdf")
plt.show()

    