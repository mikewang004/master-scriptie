from analyse7 import polymer, atom_coords
from clibraries.boxAlgorithmsInC import box_algos_lib, hoshen_kopelman_lib

import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt




def main():
    last_timestep_e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
    last_timestep_e5.read_cryst("../test_run/10e5_debug_cryst.txt")
    #last_timestep_e5.read_cryst("../test_run/debug_cryst_analyse_2D.txt") #Should be independent from actual last_timestep used
    last_timestep_e5.merge_boxes()

if __name__ == "__main__":
    main()