from analyse7 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from simulation import Simulation



def main():
    data_path = "../../data/pva-100/quick_quench/long_run"

    # cooling_e3_T07 = Simulation(0.7, -3, "%s%s" %(data_path, "/slurm-T=0.7.out"), "%s/equil_t_07_tdot_e-3_time"%(data_path))
    # cooling_e3_T07.calc_crystallisation("../crystallisation/PVA-100/quick_quench/T07/cryst_equil_T07_e-3.txt", "%s/equil_t_07_tdot_e-3"%(data_path))

    cooling_e3_T085 = Simulation(0.85, -3, "%s%s" %(data_path, "/slurm-e3-T=0.85.out"), "%s/equil_t_085_tdot_e-3_time"%(data_path), no_runs = 3)
    cooling_e3_T085.calc_avg_local_density()
    #cooling_e3_T085.calc_crystallisation("../data_online/PVA-100/quick_quench/T085/cryst_equil_T085_e-3.txt", "../data_online/PVA-100/quick_quench/T085/boxes_eigenvalues/equil_t_085_tdot_e-3")

if __name__== "__main__":
    main()