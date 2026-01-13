from analyse7 import polymer, atom_coords
from clibraries.boxAlgorithmsInC import box_algos_lib, hoshen_kopelman_lib

import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt


def compare_hk_labels(label_matrix,check_label, nridges, i,j,k):
    indices_list = []
    for i1 in [-1, 1]:
        for j1 in [-1, 1]:
            for k1 in [-1,1]:
                x = (i + i1)% nridges; y = (j + j1)% nridges; z = (k + k1) % nridges
                current_label = label_matrix[x,y,z]
                if check_label == current_label:
                    indices_list.append([x,y,z])
    return indices_list


def validate_hk_matrix(cryst_file, label_matrix, nridges = 33, cryst_cutoff = 0.8, dot_cutoff = 0.97):
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
                    if (current_row["cryst_bool"] > cryst_cutoff).bool() == True:
                        print(i,j,k)
                        indices_list = compare_hk_labels(label_matrix, label_matrix[i,j,k], nridges, i,j,k)
                        #print(indices_list)
                        if not indices_list:
                            # Check if eigenvectors are larger than the cutoff
                            pass
                        else:
                            error_counter = 0;
                            for indices in indices_list: # Used to keep track of evdot < dot_cutoff
                                x,y,z = indices
                                next_row = cryst[(cryst["xid"] == x) & (cryst["yid"] == y) & (cryst["zid"] == z)]
                                current_row_np, next_row_np = current_row.to_numpy()[0], next_row.to_numpy()[0]
                                evdot = np.abs(np.dot(current_row_np[4:], next_row_np[4:]))
                                #evdot = np.abs(current_row_np[4] * next_row_np[4]) + np.abs(current_row_np[5] + next_row_np[5]) + np.abs(current_row_np[6] * next_row_np[6])
                                if evdot < dot_cutoff:
                                    print(evdot)
                                    error_counter = error_counter + 1
                                    # Verify if other neighbours are correct, otherwise throw error 
                                else:
                                    error_counter = error_counter - 1
                            if error_counter > 0:
                                print("current index %i %i %i with label %i" %(i,j,k,label_matrix[i,j,k]))
                                print("with neighbours")
                                print(indices_list)
                                print("having too small evdot to be in a box together")

                    else: 
                        #print(current_row)
                        pass

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