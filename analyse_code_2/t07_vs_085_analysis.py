from analyse8 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import pandas as pd
from simulation import Simulation



def main():
    PVA_100_T085 = Simulation(100, "../../data/pva-100/quick_quench/e-3_T085", "../data_online/PVA-100/icryst_T085_Tdot_e-3", cooling_rate= -3, target_temp= 0.85)
    PVA_100_T07 = Simulation(100, "../../data/pva-100/quick_quench/e-3_T07", "../data_online/PVA-100/icryst_T07_Tdot_e-3", cooling_rate= -3, target_temp= 0.7)

    PVA_100_T07.domain_analysis.calc_crystallisation()


    print(PVA_100_T07.df_slurm_sim_data)
if __name__== "__main__":
    main()