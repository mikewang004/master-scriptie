import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import pandas as pd 



def sort_files_by_atom_id(filename):
    datapd = pd.read_csv(filename, sep = " ", header = None, skiprows = 11)
    datapd.columns = ["atom_id", "mol_id", "xu", "yu", "zu"]
    datapd = datapd.sort_values("atom_id")
    datapd = datapd.set_index("atom_id")
    datapd.to_csv("test.txt", sep = " ")

def sort_files_by_atom1(filename2):
    datapd = pd.read_csv(filename2, sep = " ", header = None, skiprows = 1)
    datapd.columns = ["atom1", "atom2", "bx", "by", "bz", "nx", "ny", "nz"]
    datapd = datapd.drop(["atom2"], axis = 1)
    datapd = datapd.sort_values("atom1")


    datapd = datapd.set_index("atom1")
    print(datapd)
    datapd.to_csv("test.txt", sep = " ")

def main():
    filename = "equil_t_088_tdot_e-3_time_24000000.txt"
    filename2 = "equil_t_088_tdot_e-3_time_24000000.txt_bondvecs.txt"
    sort_files_by_atom1(filename2)
    #sort_files_by_atom_id(filename)

if __name__ == "__main__":
    main()