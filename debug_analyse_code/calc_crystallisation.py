from analyse7 import * 


# Script to calculate crystallisation per data file

data_path = "../../data/pva-100/quick_quench/long_run"
# file_name_T07 = "equil_t_07_tdot_e-3"
# file_name_T08 = "equil_t_08_tdot_e-3"
# file_name_T075 = "equil_t_075_tdot_e-3"
# file_name_T085 = "equil_t_085_tdot_e-3"
# file_name_e4_T085 = "equil_t_085_tdot_e-4"
# file_name_e4_T085_run2 = "equil_t_085_tdot_e-4_run2"
# file_name_run2_e3_T07 = "equil_t_07_tdot_e-3_run2"
file_name_run2_e3_T085 = "equil_t_085_tdot_e-3_run2"

# # Search all 

# time_temp_T08 = get_time_temp_from_slurm(data_path + "/slurm-T=0.8.out")
# time_temp_T07 = get_time_temp_from_slurm(data_path + "/slurm-T=0.7.out")
# time_temp_T085 = get_time_temp_from_slurm(data_path + "/slurm-T=0.85.out")
# time_temp_T075 = get_time_temp_from_slurm(data_path + "/slurm-T=0.75.out")
# time_temp_e4_T085 = get_time_temp_from_slurm(data_path + "/slurm-e4-T=0.85.out")
# #time_temp_run2_e3_T07 = get_time_temp_from_slurm(data_path + "/slurm-run2-e3-T=0.7.out")
time_temp_run2_e3_T085 = get_time_temp_from_slurm(data_path + "/slurm-run2-e3-T=0.85.out")
# time_temp_e4_T085_run2 = get_time_temp_from_slurm(data_path + "/equil_t_085_tdot_e-4_run2.out")
#print(time_temp_T07[:, 0])

cryst_path_prefix = "../crystallisation"
cryst_T07 = "%s/all_times_cryst_equil_t_07_tdot_e-3_full_array.txt" %cryst_path_prefix
cryst_T08 = "%s/all_times_cryst_equil_t_08_tdot_e-3.txt" %cryst_path_prefix
cryst_T075 = "%s/all_times_cryst_equil_t_075_tdot_e-3.txt" %cryst_path_prefix
cryst_T085 = "%s/all_times_cryst_equil_t_085_tdot_e-3_full_array.txt" %cryst_path_prefix
cryst_T085_tdot_4 = "%s/all_times_cryst_equil_t_085_tdot_e-4_full_array.txt" %cryst_path_prefix
#e3_run2_cryst_T07 = "all_times_cryst_equil_t_07_tdot_e-3_run2.txt"
e3_run2_cryst_T085 = "all_times_cryst_equil_t_085_tdot_e-3_run2.txt"
#e4_cryst_T085 = "all_times_cryst_equil_t_085_tdot_e-4.txt"
#cryst_list = [cryst_T07, cryst_T075, cryst_T08, cryst_T085]
#temps = [0.7, 0.75, 0.8, 0.85]

cryst_list = [cryst_T07, cryst_T085]
temps = [0.7, 0.85]
#cryst_list = [cryst_T085, cryst_T085_tdot_4]
#temps = [3, 4]


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
        plt.scatter(time, cryst, label = "temperature = %.2f" %temps[i])
        plt.scatter(time, cryst, label = r"cooling rate = $10^{-%i}$" %temps[i])
        i = i + 1
    
    plt.title(r"PVA-100 crystallisation after quick quench at cooling rate $\dot{T} = 10^{-3} \tau^{-1}$")
    #plt.title(r"PVA-100 crystallisation after quick quench at $T = 0.85$")
    plt.xlabel(r"time ($\tau$)")
    plt.ylabel("fraction of crystallinity")
    plt.legend()
    plt.savefig("isothermal_cryst_quick_quench_t07_085.pdf")
    plt.show()


def merge_cryst_files():
    """Help function to merge two or more crystallisation files to one"""
    test = np.loadtxt("../crystallisation/all_times_cryst_equil_t_085_tdot_e-3_run2.txt", delimiter=",")
    print(test[:,0])
    test[:, 0] = test[:,0] + 60000000
    np.savetxt("temp.txt", test, delimiter = ",")


#plot_crystallisation(cryst_list, temps)

#calc_crystallisation(data_path, time_temp_T075[:, 0], file_name_T075)

#calc_crystallisation(data_path, time_temp_T085[:, 0], file_name_T085)

#calc_crystallisation(data_path, time_temp_e4_T085[:, 0], file_name_e4_T085)

calc_crystallisation(data_path, time_temp_run2_e3_T085[:, 0], file_name_run2_e3_T085)

#calc_crystallisation(data_path, time_temp_run2_e3_T07[:, 0], file_name_run2_e3_T07)

#calc_crystallisation(data_path, time_temp_e4_T085_run2[:, 0], file_name_e4_T085_run2)

#merge_cryst_files()


