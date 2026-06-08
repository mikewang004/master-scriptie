import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import pandas as pd




class Polymer:

    def __init__(self, polymer_length, path_to_file):
        self.datapd = self.read_in_file(path_to_file)
        self.datapd = self.replace_mol_id_by_polymer_length()
        self.n_atoms = len(self.datapd.index)
        self.polymer_length = polymer_length
        self.no_polymers = int(self.n_atoms/self.polymer_length)
        self.timestep = self.get_timestep_from_file_name(path_to_file)

        self.volume, self.boxlengths, self.dimensions = self.get_volume_box(path_to_file)


    def read_in_file(self, filename):
        datapd = pd.read_csv(filename, sep = " ", header = None, skiprows = 9, dtype={1: "int32", 2: "float64", 3: "float64", 4: "float64"})
        datapd.columns = ["atom_id", "mol_id", "xu", "yu", "zu"]
        datapd = datapd.set_index("atom_id")
        
        datapd = datapd.sort_values("atom_id")

        return datapd

    def replace_mol_id_by_polymer_length(self):
        self.datapd['mol_id']  = (np.arange(len(self.datapd)) // self.polymer_length)
        return self.datapd


    def get_volume_box(self, path_to_file):
        """Calculates volume of the total (!) simulation box."""
        datapd_first_rows = (pd.read_csv(path_to_file, sep = " ", header = None, skiprows = 5, nrows = 3))
        datapd_first_rows = datapd_first_rows.rename({0: "min", 1: "max"}, axis = 1).rename({0: "x", 1: "y", 2:"z"}, axis = 0)
        length = datapd_first_rows["max"] - datapd_first_rows["min"]
        volume = length.loc["x"]*length.loc["y"]*length.loc["z"]
        return volume, length, datapd_first_rows


    def get_timestep_from_file_name(self, path_to_file):
        """Only reads timestep from the file name."""
        pattern = r"_(\d+)\.txt$"
        timestep = re.search(pattern, path_to_file).group(1)
        return int(timestep)