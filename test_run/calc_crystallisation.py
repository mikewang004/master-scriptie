from analyse6 import * 


# Script to calculate crystallisation per data file

data_path = "../../data/pva-100/quick_quench/long_run"
file_name_T07 = "equil_t_07_tdot_e-3"
file_name_T08 = "equil_t_08_tdot_e-3"
file_name_T075 = "equil_t_075_tdot_e-3"
file_name_T085 = "equil_t_085_tdot_e-3"

# Search all 

time_temp_T08 = get_time_temp_from_slurm(data_path + "/slurm-T=0.8.out")
time_temp_T07 = get_time_temp_from_slurm(data_path + "/slurm-T=0.7.out")
time_temp_T085 = get_time_temp_from_slurm(data_path + "/slurm-T=0.85.out")
time_temp_T075 = get_time_temp_from_slurm(data_path + "/slurm-T=0.75.out")
print(time_temp_T07[:, 0])

cryst_T07 = "all_times_cryst_equil_t_07_tdot_e-3.txt"
cryst_T08 = "all_times_cryst_equil_t_08_tdot_e-3.txt"
cryst_T075 = "all_times_cryst_equil_t_075_tdot_e-3.txt"
cryst_T085 = "all_times_cryst_equil_t_085_tdot_e-3.txt"
cryst_list = [cryst_T07, cryst_T08]
temps = [0.7, 0.8]


def calc_crystallisation(file_name_path, times, cryst_file_name):
    """Calculates the crystallisation of all timesteps of a run and saves that to a file"""
    cryst_array_string = "%s/all_times_cryst_%s.txt" %(file_name_path, cryst_file_name)
    with open(cryst_array_string, "w") as file:
        file.write("")
    for time in times:
        current_file_name = "%s/%s_time_%i.txt" %(file_name_path, cryst_file_name, time)
        current_atom_coords = atom_coords(current_file_name)
        frac_cryst = current_atom_coords.get_nematic_vector_4(save_ev= True, save_string = "%s/cryst_%s_time_%i.txt" %(file_name_path, cryst_file_name, time))
        print(current_file_name)
        with open(cryst_array_string, "a") as file:
            file.write(f"{time}, {frac_cryst}\n")


def plot_crystallisation(list_cryst_file_names, temps):
    i = 0
    for item in list_cryst_file_names:
        array = np.loadtxt(item, delimiter="," )
        print(array)
        time = array[:, 0]; cryst = array[:, 1]
        plt.scatter(time, cryst, label = "temperature = %f" %temps[i])
        i = i + 1
    
    plt.title("Crystallisation as function of time")
    plt.xlabel(r"time ($\tau$)")
    plt.ylabel("fraction of crystallinity")
    plt.legend()
    plt.show()



#plot_crystallisation(cryst_list, temps)

#calc_crystallisation(data_path, time_temp_T075[:, 0], file_name_T075)

calc_crystallisation(data_path, time_temp_T085[:, 0], file_name_T085)
