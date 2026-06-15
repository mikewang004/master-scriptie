from analyse9 import polymer, atom_coords
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import pathlib
import re
from typing import List, Tuple, Optional
import bisect
from simulation import Simulation



def plot_crystallisation_vs_volume():
    pass

def crystallisation_vs_time(simulations: list):
    for simulation in simulations:
        #plt.scatter(simulation.)
        print(simulation.df_slurm_sim_data)


def main():
    PVA_100 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3")
    PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")
    times = np.array([0, 5, 10, 30])*1200000

    for i in range(0, len(times)):
        print(times[i])

if __name__== "__main__":
    main()