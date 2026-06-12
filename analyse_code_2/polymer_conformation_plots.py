from analyse9 import polymer, atom_coords
from simulation import Simulation
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



#TODO: make plot Rg vs time for different chain lengths, show time evolution 
#Make plot Rg, Re for different stages of coarsening for N = 100, 1000



def main():
    PVA_100 = Simulation(100, "../../data/pva-100/quick_quench/equil", "../data_online/PVA-100/icryst_T088_Tdot_e-3")
    PVA_1000 = Simulation(1000, "../../data/PVA-1000/equil", "../data_online/PVA-1000/icryst_T088_Tdot_e-3")
    times = (np.array([0, 5, 10, 30])*1.2e6)

    PVA_100_times = PVA_100.get_mutiple_polymers_by_time(times)
    
    #Calc Rg per time 
    for polymer in PVA_100_times:
        Rg = polymer.gyration_radius()
        print(Rg)

if __name__== "__main__":
    main()