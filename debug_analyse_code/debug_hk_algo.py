from analyse7 import polymer, atom_coords
from clibraries.boxAlgorithmsInC import box_algos_lib, hoshen_kopelman_lib

import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt


def validate_hk_matrix(cryst_file, label_matrix, nridges = 33):
    """Checks if the label_matrix satifies following requirements: 
        1. Each point with a label has crystallinity > 0.8;
        2. Neighboring points, if same label, satisfy ev * ev \geq 0.97"""
    cryst = cryst_file
    for i in range(0, nridges):
        for j in range(0, nridges):
            for k in range(0, nridges):
                if label_matrix[i,j,k] != 0:
                    #print(i,j,k)
                    # Find corresponding entry in spreadsheet
                    row_mask = (cryst["xid"] == i) & (cryst["yid"] == j) & (cryst["zid"] == k)
                    current_row = cryst[row_mask]
                    print(current_row)
                   # print(cryst[row_mask])
                    if current_row["cryst_bool"] < 0.8:
                        print(current_row)

                    # Get neighboring rows
                    
                    

def main():
    last_timestep_e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
    #last_timestep_e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
    nridges = 33
    if nridges == 6:
        last_timestep_e5.read_cryst("./debug_cryst_analyse_2D.txt") #Should be independent from actual last_timestep used
    else:
        last_timestep_e5.read_cryst("../test_run/10e5_debug_cryst.txt")
        label_matrix = np.load("hk_label_matrix.npy")
        #last_timestep_e5.read_cryst("../test_run/10e5_T1_debug_cryst.txt")
        ndriges = 33
        validate_hk_matrix(last_timestep_e5.df_cryst, label_matrix)
    #last_timestep_e5.merge_boxes(nridges = nridges)

if __name__ == "__main__":
    main()