from analyse7 import polymer, atom_coords
from clibraries.boxAlgorithmsInC import box_algos_lib, hoshen_kopelman_lib

import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt


def compare_hk_labels(label_matrix,check_label, nridges, i,j,k):
    indices_list = []
    i1, j1, k1 = 0,0,0
    for i1 in [-1, 1]:
        compare_hk_labels_subroutine(label_matrix,check_label, indices_list, nridges, i,j,k, i1, j1, k1)
    i1, j1, k1 = 0,0,0
    for j1 in [-1, 1]:
        compare_hk_labels_subroutine(label_matrix,check_label, indices_list, nridges, i,j,k, i1, j1, k1)
    i1, j1, k1 = 0,0,0
    for k1 in [-1,1]:
        compare_hk_labels_subroutine(label_matrix,check_label, indices_list, nridges, i,j,k, i1, j1, k1)

    return indices_list


def compare_hk_labels_subroutine(label_matrix,check_label, indices_list, nridges, i,j,k, i1, j1, k1):
    """Repeating part of compare_hk_labels"""
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
    master_list_errorenous_entries = []
    master_list_indices_ev = []
    for i in range(0, nridges):
        for j in range(0, nridges):
            for k in range(0, nridges):
                if label_matrix[i,j,k] != 0:
                    #print(i,j,k)
                    # Find corresponding entry in spreadsheet
                    row_mask = (cryst["xid"] == i) & (cryst["yid"] == j) & (cryst["zid"] == k)
                    current_row = cryst[row_mask]
                    if (current_row["cryst_bool"] > cryst_cutoff).bool() == True:
                        #print(i,j,k)
                        indices_list = compare_hk_labels(label_matrix, label_matrix[i,j,k], nridges, i,j,k)
                        #print(indices_list)
                        ev_dot_list_current = []
                        if not indices_list:
                            # Check if eigenvectors are larger than the cutoff
                            pass
                        else:
                            error_bool = True
                            for indices in indices_list: # Used to keep track of evdot < dot_cutoff
                                x,y,z = indices
                                next_row = cryst[(cryst["xid"] == x) & (cryst["yid"] == y) & (cryst["zid"] == z)]
                                current_row_np, next_row_np = current_row.to_numpy()[0], next_row.to_numpy()[0]
                                #evdot = np.abs(np.dot(current_row_np[4:], next_row_np[4:]))
                                evdot = np.abs(current_row_np[4] * next_row_np[4]) + np.abs(current_row_np[5] * next_row_np[5]) + np.abs(current_row_np[6] * next_row_np[6])
                                ev_dot_list_current.append(evdot)
                                if evdot < dot_cutoff:
                                    #print("Current value under cutoff")
                                    #print(current_row, next_row)
                                    #print(evdot)
                                    pass
                                else:
                                    error_bool = False
                            if error_bool == True: #Only throw error if all neighbours are incorrect
                                # print("current index %i %i %i with label %i" %(i,j,k,label_matrix[i,j,k]))
                                # print("with neighbours")
                                # print(indices_list)
                                # print("and eigenvector-dot values")
                                # print(ev_dot_list_current)
                                # print("having too small evdot to be in a box together")
                                master_list_errorenous_entries.append([i,j,k])
                                indices_ev_dot_array = np.zeros([len(indices_list), 4])
                                for i in range(len(indices_list)):
                                    x,y,z = indices_list[i]
                                    indices_ev_dot_array[i,:3] = x,y,z
                                    indices_ev_dot_array[i,3] = ev_dot_list_current[i]
                                print(indices_ev_dot_array)
                                master_list_indices_ev.append(indices_ev_dot_array)


                    else: 
                        print(current_row)
                        break

                    # Get neighboring rows
    print(master_list_errorenous_entries)
    print(master_list_indices_ev)
                    

def print_hk_matrix(label_matrix, nridges =33):
    """Reads a hk-matrix and prints it in 2d"""
    for i in range(0, nridges):
        for j in range(0, nridges):
            for k in range(0, nridges):
                print("%i %i %i %d " %(i,j,k,label_matrix[i,j,k]))
                    

def main():
    last_timestep_e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
    #last_timestep_e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
    nridges = 33
    if nridges == 6:
        last_timestep_e5.read_cryst("./debug_cryst_analyse_2D.txt") #Should be independent from actual last_timestep used
    else:
        last_timestep_e5.read_cryst("../test_run/10e5_debug_cryst.txt")
        #label_matrix = np.load("hk_label_matrix.npy")
        #last_timestep_e5.read_cryst("../test_run/10e5_T1_debug_cryst.txt")
        ndriges = 33
        #print_hk_matrix(label_matrix, nridges)
        #validate_hk_matrix(last_timestep_e5.df_cryst, label_matrix)
    last_timestep_e5.merge_boxes(nridges = nridges)

if __name__ == "__main__":
    main()