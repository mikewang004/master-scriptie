import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import pandas as pd
import re 
from numba import jit
from tqdm import tqdm

from nematic_vector import calc_nematic_tensor_2, nematic_vector_loop, nematic_vector_loop_2, compute_Q, orderparameter
from hoshenKopelmanInPython import hk_in_python

def wrap_coordinates(column, minlength, total_length):
    return (column - minlength) % total_length + minlength

def fraction_crystallinity(data, n_ridges_3d, cutoff = 0.8):
    """Data 1d 1d list/array, defined as data > 0.8 -> crystallinity = 1;
    data <= 0.8 -> crystallinity = 0 as in Sommer/Luo Sep 2010"""
    mask = data > cutoff
    fraction = len(data[mask]) / n_ridges_3d
    return fraction, len(data[mask])

class results(object):
    def __init__(self):
        return



class atom_coords:

    """Read in files and do basic data preprocessing"""

    def __init__(self, path_to_file, cell_length = 2, polymer_length = 100):
        #Read in data 
        self.path_to_file = path_to_file
        self.datapd = self.read_in_file(path_to_file)
        self.n_atoms = len(self.datapd.index)
        self.polymer_length = polymer_length
        self.no_polymers = int(self.n_atoms/self.polymer_length)
        self.datapd = self.replace_mol_id_by_polymer_length()
        self.volume, self.boxlengths, self.dimensions = self.get_volume_box(path_to_file)
        self.cell_length = cell_length
        # self.nridges, self.no_nridges_3d = self.calc_nridges()
        self.current_timestep = self.get_timestep_from_file_name(path_to_file) 
        self.wrapped_monomers = self.wrap_coordinates_all_data()
        #Calculate box properties 

        self.bond_vectors = self.calculate_bond_vectors()
        self.bond_vectors, self.nridges, self.no_nridges_3d  = self.make_cell_grid()
        #self.df_cryst = self.get_nematic_vector_5()
        self.results = results()


    def read_in_file(self, filename):
        datapd = pd.read_csv(filename, sep = " ", header = None, skiprows = 9, dtype={1: "int32", 2: "float64", 3: "float64", 4: "float64"})
        datapd.columns = ["atom_id", "mol_id", "xu", "yu", "zu"]
        datapd = datapd.set_index("atom_id")
        
        datapd = datapd.sort_values("atom_id")
        #print(datapd)
        #datapd.to_csv("pva-200-polyno0.txt", sep = " ")
        #datapd = datapd.sort_values(["mol_id", "atom_id"])
        #

        return datapd

    def replace_mol_id_by_polymer_length(self):
        #self.datapd['mol_id'] = (np.arange(len(self.datapd)) // self.polymer_length)
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

    def wrap_coordinates_all_data(self):
        """Converts coordinates to wrapped coordinates. Input can be float or array (float)"""
        wrapped_coords = self.datapd
        wrapped_coords = wrapped_coords.copy()
        wrapped_coords["xu"] = wrap_coordinates(self.datapd["xu"], self.dimensions.loc["x", "min"], self.boxlengths["x"])
        wrapped_coords["yu"] = wrap_coordinates(self.datapd["yu"], self.dimensions.loc["y", "min"], self.boxlengths["y"])
        wrapped_coords["zu"] = wrap_coordinates(self.datapd["zu"], self.dimensions.loc["z", "min"], self.boxlengths["z"])
        return wrapped_coords

    def calculate_bond_vectors(self):

        shifted = self.datapd.groupby('mol_id')[['xu', 'yu', 'zu']].shift(-1)
        shifted = shifted.copy()
        bond_vecs = shifted - self.datapd[['xu', 'yu', 'zu']]

        bond_vecs.columns = ['bx', 'by', 'bz']
        bond_vecs = bond_vecs.dropna()
        norms = np.linalg.norm(bond_vecs[['bx','by','bz']].to_numpy(), axis=1)
        bond_vecs[['bx','by','bz']] = bond_vecs[['bx','by','bz']].div(norms, axis=0)
        bond_vecs = bond_vecs.apply(np.float32)
        return bond_vecs


    def return_positive_midpoint_coords(self):
        pos_coords = self.wrapped_monomers[['mol_id', 'xu', 'yu', 'zu']]# - self.dimensions.iloc[:, 0]
        pos_coords = pos_coords.copy()
        pos_coords["xu"] = pos_coords["xu"] - self.dimensions.loc["x", "min"]
        pos_coords["yu"] = pos_coords["yu"] - self.dimensions.loc["y", "min"]
        pos_coords["zu"] = pos_coords["zu"] - self.dimensions.loc["z", "min"]
        #pos_coords.to_csv("test3.txt", sep = " ")
        # pos_coords.to_csv("pva-100-t088-poly20-pos_coords.txt", sep = " ")
        # #print(self.dimensions.iloc[:, :])

        shifted = pos_coords.groupby('mol_id')[['xu', 'yu', 'zu']].shift(-1).copy() #check if this is valid
        # print(pos_coords)
        # print(shifted)
        xm = (shifted + pos_coords[["xu", "yu", "zu"]])/2.0 #bond vector position
        xm = xm.dropna()



        return xm

    def return_positive_midpoint_coords_c_like(self):
        # xlo,ylo,zlo should match what C reads from the dump
        xlo = self.dimensions.loc["x", "min"]
        ylo = self.dimensions.loc["y", "min"]
        zlo = self.dimensions.loc["z", "min"]

        # use the same base coordinates as C: xu, yu, zu from the dump
        pos_coords = self.datapd[['mol_id', 'xu', 'yu', 'zu']].copy()
        print(pos_coords)
        # shift so that box min is at 0 (like xp = xu - xlo in C)
        pos_coords['xu'] = pos_coords['xu'] - xlo
        pos_coords['yu'] = pos_coords['yu'] - ylo
        pos_coords['zu'] = pos_coords['zu'] - zlo

        # neighbor within each chain (same as C atom1/atom2)
        shifted = pos_coords.groupby('mol_id')[['xu', 'yu', 'zu']].shift(-1)

        # bond midpoints (like xm = (xp1 + xp2)/2)
        xm = (shifted + pos_coords[['xu', 'yu', 'zu']]) / 2.0
        xm = xm.dropna()
        print(xm)
        return xm

    def calc_nridges(self):
        nridges = (self.boxlengths/self.cell_length).astype(int)
        nridges_total_3d = (nridges.iloc[0]*nridges.iloc[1]*nridges.iloc[2]).item()
        return nridges, nridges_total_3d
    
    def make_cell_grid(self, cell_length = 2):
        nridges, nridges_total_3d = self.calc_nridges()
        actual_cell_length = self.boxlengths/nridges
        
        #Shift monomers to start from 0 
        xm = self.return_positive_midpoint_coords().apply(np.float32) #Modify this binning to use (int) xm / lx 
        nx = (xm["xu"] / actual_cell_length["x"]).astype(int)
        ny = (xm["yu"] / actual_cell_length["y"]).astype(int)
        nz = (xm["zu"] / actual_cell_length["z"]).astype(int)
        #xm.to_csv("test.txt", sep = " ")
        self.bond_vectors = self.bond_vectors.assign(
            xm = xm["xu"],
            ym = xm["yu"],
            zm = xm["zu"],
            nx=nx,
            ny=ny,
            nz=nz,
        )
        #self.bond_vectors = self.bond_vectors[(self.bond_vectors["nx"] >= 0) & (self.bond_vectors["ny"] >= 0) & (self.bond_vectors["nz"] >= 0)]
        return self.bond_vectors, nridges, nridges_total_3d

    def get_nematic_vector_5(self, save_string = None, cryst_cutoff = 0.8):
        #self.df_cryst = nematic_vector_loop(data, self.bond_vectors)

        self.df_cryst = (
            self.bond_vectors.groupby(['nx', 'ny', 'nz'])
            #.apply(compute_Q)
            .apply(orderparameter)
            .reset_index()
        )
        print(self.df_cryst)
        self.fraction_crystallinity, _ = fraction_crystallinity(self.df_cryst.iloc[:,3], self.no_nridges_3d, cutoff= cryst_cutoff)
        if isinstance(save_string, str):
            self.df_cryst.to_csv("%s" %save_string, sep = " ", mode = "w")
        print(self.fraction_crystallinity)
        return self.fraction_crystallinity


    def get_closest_monomers(self, atom_id: int, n: int = 10):
        """
        Return:
        mol_id of this atom,
        and a DataFrame with atom_id, mol_id, xu, yu, zu, distance
        for the n closest monomers in the same mol_id.
        """

        # 1. Get central atom row
        datapd = self.datapd
        try:
            center = datapd.loc[atom_id]
        except KeyError:
            raise ValueError(f"atom_id {atom_id} not found in DataFrame index")

        mol_id = center["mol_id"]
        center_pos = center[["xu", "yu", "zu"]].to_numpy(dtype=float)

        # 2. All atoms with same mol_id
        same_mol = datapd[datapd["mol_id"] == mol_id].copy()

        # 3. Compute distances
        coords = same_mol[["xu", "yu", "zu"]].to_numpy(dtype=float)
        diffs = coords - center_pos
        dists = np.sqrt((diffs ** 2).sum(axis=1))
        same_mol["distance"] = dists

        # 4. Remove itself
        same_mol = same_mol[same_mol.index != atom_id]

        # 5. Take n closest
        closest = same_mol.nsmallest(n, "distance")

        # Optionally add atom_id as a column (from index)
        closest = closest.reset_index().rename(columns={"index": "atom_id"})

        return mol_id, closest



class polymer():
    def __init__(self, path_to_file, polymer_length = 100):
        self.atom_coords = atom_coords(path_to_file, polymer_length= polymer_length)
        self.results = results()

    def read_cryst(self, location):
        self.df_cryst = pd.read_csv(location, sep = " ", header = None, skiprows = 1, index_col = False).iloc[:, 1:]
        self.df_cryst.columns = ["xid", "yid", "zid", "cryst_bool", "x_ev", "y_ev", "z_ev"]
        self.frac_cryst, _ = fraction_crystallinity(self.df_cryst.iloc[:, 3], self.atom_coords.no_nridges_3d)
        return 0;

    def merge_boxes_2(self, ndot_cutoff = 0.97,cryst_cutoff = 0.8, save = False, print_results: bool = False):
        label_matrix = hk_in_python(self.df_cryst, ndot_cutoff = ndot_cutoff, nridges = self.atom_coords.nridges, cryst_cutoff = cryst_cutoff)
        total_box_elements = (self.atom_coords.nridges["x"]*self.atom_coords.nridges["y"]*self.atom_coords.nridges["z"]).astype(int)
        unique_values, counts = np.unique(label_matrix, return_counts=True) #Labels and how much each label occurs


        total_number_merged_clusters = counts[counts > 1]
        self.results.total_number_clusters = total_number_merged_clusters.size
        self.results.total_number_independent_clusters = unique_values.size-1
        cluster_counts = counts[1:]
        self.results.mean_cluster_size = np.mean(cluster_counts)
        self.results.total_number_crystalline_grid_elements = np.sum(cluster_counts)
        self.results.fraction_crystallinity, no_crystalline_elements = fraction_crystallinity(self.df_cryst.iloc[:,3], self.atom_coords.no_nridges_3d)
        if print_results == True:
            #print("For time %i"%(self.atom_coords.current_timestep))
            print("total number clusters w/ >= 2 elements: %i" %(self.results.total_number_clusters))
            print("total number independent crystalline domains: %i" %(self.results.total_number_independent_clusters))
            print("average cluster size crystalline domains: %f" %(self.results.mean_cluster_size))
            print("total number crystalline/all grid elements: %i/%i -> cryt_frac = %f" 
                %(self.results.total_number_crystalline_grid_elements,total_box_elements, 
                    self.results.total_number_crystalline_grid_elements/total_box_elements))
            print("earlier calculated frac_cryst = %f" %(self.frac_cryst))
        return label_matrix


    def create_new_polymer_df(self, column_names):
        """Creates an empty dataframe per mol_id with len(column_names) columns. column_names must be a list"""
        new_df = pd.DataFrame(np.zeros([self.atom_coords.no_polymers, len(column_names)]), columns = column_names)
        new_df.index.name = "mol_id"
        new_df.index = new_df.index + 1
        return new_df

    def end_to_end_distance(self, show_plot = False, save_plot_string = None):
        # Calculates end-to-end distance of each polymer 
        #print(self.datapd)
        df_end_end_length = self.create_new_polymer_df(["end_end_length"])
        for i in tqdm(range(0, self.atom_coords.no_polymers)):
            # First calculate end to end distance 
            # defined as r_n - r_i for each position 
            #print(i, (i*self.atom_coords.polymer_length-1), (i+1)*(self.atom_coords.polymer_length-1))
            #subset = self.atom_coords.datapd.iloc[(i*(self.atom_coords.polymer_length-1)):((i+1)*(self.atom_coords.polymer_length-1)), 1:4]
            first_in_chain = self.atom_coords.datapd.iloc[(i*self.atom_coords.polymer_length), 1:4]
            last_in_chain = self.atom_coords.datapd.iloc[(i+1)*(self.atom_coords.polymer_length)-1, 1:4]
            #print(subset)
            dist = last_in_chain - first_in_chain
            df_end_end_length.iloc[i] = (dist.iloc[0]* dist.iloc[0] + dist.iloc[1] * dist.iloc[1] + dist.iloc[2] * dist.iloc[2])
        end_to_end_length = (np.sum(df_end_end_length.to_numpy())/self.atom_coords.no_polymers)
        end_end_distance_normalised = np.sqrt(df_end_end_length.iloc[:])
        self.results.end_to_end_distribution = (end_end_distance_normalised)
        self.results.end_to_end_length = end_to_end_length
        self.results.mean_squared_end_to_end = np.sqrt(end_to_end_length)
        print("mean end-to-end is %f" %self.results.mean_squared_end_to_end)
        return end_to_end_length

    def gyration_radius(self, nridges = 33, show_plot = False):
        """Should confirm to Saras 2018 paper eq. 3"""
        df_gyration_radius = self.create_new_polymer_df(["comx", "comy", "comz", "gyration_radius"])
        counter = 0
        for i in range(0, self.atom_coords.no_polymers):
            subset = self.atom_coords.datapd[(self.atom_coords.datapd["mol_id"] == i+1)].iloc[:, 1:4]

            com = np.mean(subset, axis = 0)

            df_gyration_radius.iloc[i, :3] = com
            # Shift system to have new center of mass as center 
            subset_com = subset - com
            gyration_radius_squared = np.sum((subset_com**2), axis = 1)
            #print(np.mean(gyration_radius_squared))
            df_gyration_radius.iloc[i, 3] = np.mean(gyration_radius_squared) 
        self.results.gyration_radius_distribution = np.sqrt(df_gyration_radius.iloc[:, 3])
        self.results.gyration_radius = df_gyration_radius
        self.results.mean_gyration_radius = np.mean(df_gyration_radius["gyration_radius"]) #Ensemble average
        print("mean gyration is %f" %np.sqrt(self.results.mean_gyration_radius))