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
    PVA_100_T075 = Simulation(100, "../../data/pva-100/quick_quench/e-3_T075", "../data_online/PVA-100/icryst_T075_Tdot_e-3", cooling_rate= -3, target_temp= 0.75)
    PVA_100_T08 = Simulation(100, "../../data/pva-100/quick_quench/e-3_T08", "../data_online/PVA-100/icryst_T08_Tdot_e-3", cooling_rate= -3, target_temp= 0.8)
    PVA_100_T08 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3", cooling_rate= -3, target_temp= 0.88)
    simulations = [PVA_100_T07, PVA_100_T085]
    times = np.array([0, 5, 10, 30])*1200000

    for i in range(0, len(times)):
        print(times[i])

    #PVA_100_T07.domain_analysis.calc_crystallisation()


    #print(PVA_100_T07.domain_analysis.read_avg_domain_size())
if __name__== "__main__":
    main()