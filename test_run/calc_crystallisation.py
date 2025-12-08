from analyse6 import * 


# Script to calculate crystallisation per data file

data_path = "../../data/pva-100/quick_quench/long_run"
file_name_T07 = "equil_t_07_tdot_e-3"
file_name_T08 = "equil_t_08_tdot_e-3"

# Search all 

time_temp_T08 = get_time_temp_from_slurm(data_path + "/slurm-T=0.8.out")
time_temp_T07 = get_time_temp_from_slurm(data_path + "/slurm-T=0.7.out")
print(time_temp_T07[:, 0])


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

calc_crystallisation(data_path, time_temp_T07[:, 0], file_name_T07)

calc_crystallisation(data_path, time_temp_T08[:, 0], file_name_T08)
