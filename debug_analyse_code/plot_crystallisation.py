from analyse7 import * 

cryst_path_prefix = "../crystallisation"
cryst_T07 = "%s/all_times_cryst_equil_t_07_tdot_e-3.txt" %cryst_path_prefix
cryst_T08 = "%s/all_times_cryst_equil_t_08_tdot_e-3.txt" %cryst_path_prefix
cryst_T075 = "%s/all_times_cryst_equil_t_075_tdot_e-3.txt" %cryst_path_prefix
cryst_T085 = "%s/all_times_cryst_equil_t_085_tdot_e-3.txt" %cryst_path_prefix

cryst_list = [cryst_T07, cryst_T085]
temps = [0.7, 0.85]


def merge_cryst_runs(name_cryst):
    """Merges crystallisation files to one file."""
    first_file = np.loadtxt(name_cryst)

    

    np.savetxt("%s_full_run.txt", final_array)
    return 0;



def plot_crystallisation(list_cryst_file_names, temps):
    i = 0
    for item in list_cryst_file_names:
        array = np.loadtxt(item, delimiter="," )
        print(array)
        time = array[:, 0]; cryst = array[:, 1]
        plt.scatter(time, cryst, label = "temperature = %.2f" %temps[i])
        i = i + 1
    
    plt.title(r"PVA-100 crystallisation after quick quench at cooling rate $\dot{T} = 10^{-3} \tau^{-1}$")
    plt.xlabel(r"time ($\tau^{-1}$)")
    plt.ylabel("fraction of crystallinity")
    plt.legend()
    plt.savefig("isothermal_cryst_quick_quench_t_07_085.pdf")
    plt.show()






def main():
    plot_crystallisation(cryst_list, temps)

if __name__ == "__main__":
    main()



