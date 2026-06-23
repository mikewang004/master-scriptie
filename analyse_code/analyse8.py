import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import re 
import pandas as pd 
from tqdm import tqdm
import functools
import sys
from clibraries.boxAlgorithmsInC import box_algos_lib, hk_lib
import ctypes
import time
from numba import jit
from nematic_vector import calc_nematic_tensor_2, nematic_vector_loop, nematic_vector_loop_2, compute_Q, orderparameter
from matplotlib.colors import ListedColormap, BoundaryNorm
from hoshenKopelmanInPython import hk_in_python
from matplotlib.patches import Rectangle

def wrap_coordinates(column, minlength, total_length):
    return (column - minlength) % total_length + minlength

def fraction_crystallinity(data, n_ridges_3d, cutoff = 0.8):
    """Data 1d 1d list/array, defined as data > 0.8 -> crystallinity = 1;
    data <= 0.8 -> crystallinity = 0 as in Sommer/Luo Sep 2010"""
    mask = data > cutoff
    fraction = len(data[mask]) / n_ridges_3d
    return fraction, len(data[mask])

def gaussian(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0)**2 / (2 * sigma**2))


def get_time_temp_from_slurm(file_to_path):
    n_columns = 2
    with open(file_to_path, "r") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if "Per MPI rank memory allocation (min/avg/max)" in line:
            start_index = i+2
            break
    collected_rows = []
    for line in lines[start_index:]:
        if "error: *** JOB" in line:
            print(line)
            break
        row = line.strip().split()
        collected_rows.append(row[:n_columns])
    return np.array(collected_rows, dtype = float)


def calc_bond_angle(r1, r2):
    r = np.dot(r1, r2)
    num = (np.linalg.norm(r1) * np.linalg.norm(r2))
    #print(r, num)
    #print(np.degrees(np.arccos(r/ num)))
    return np.degrees(np.arccos(r/ num))

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
        #pos_coords = self.datapd[['mol_id', 'xu', 'yu', 'zu']]

        #pos_coords.to_csv("test2.txt", sep = " ")
        # pos_coords = pos_coords.copy()
        # pos_coords["xu"] = pos_coords["xu"] - self.dimensions.loc["x", "min"]
        # pos_coords["yu"] = pos_coords["yu"] - self.dimensions.loc["y", "min"]
        # pos_coords["zu"] = pos_coords["zu"] - self.dimensions.loc["z", "min"]
        # # #print(self.dimensions.iloc[:, :])


        # shifted = pos_coords.groupby('mol_id')[['xu', 'yu', 'zu']].shift(-1)
        # xm = (shifted + pos_coords[["xu", "yu", "zu"]])/2 #bond vector position

        #pos_coords = self.datapd[['mol_id', 'xu', 'yu', 'zu']]# - self.dimensions.iloc[:, 0]
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
        # xm["xu"] = wrap_coordinates(xm["xu"], self.dimensions.loc["x", "min"], self.boxlengths["x"])
        # xm["yu"] = wrap_coordinates(xm["yu"], self.dimensions.loc["y", "min"], self.boxlengths["y"])
        # xm["zu"] = wrap_coordinates(xm["zu"], self.dimensions.loc["z", "min"], self.boxlengths["z"])


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
        #print(self.bond_vectors)
        xm = self.return_positive_midpoint_coords().apply(np.float32) #Modify this binning to use (int) xm / lx 
        #gridx = np.linspace(0, np.abs(self.dimensions.loc["x", "min"]) +np.abs(self.dimensions.loc["x", "max"]), nridges.loc["x"]+1)
        #gridy = np.linspace(0,np.abs(self.dimensions.loc["y", "min"])+ np.abs(self.dimensions.loc["y", "max"]), nridges.loc["y"]+1)
        #gridz = np.linspace(0,np.abs(self.dimensions.loc["z", "min"])+ np.abs(self.dimensions.loc["z", "max"]), nridges.loc["z"]+1)
        #nx = pd.cut(xm["xu"], bins = gridx, right = False, labels = False)
        #ny = pd.cut(xm["yu"], bins = gridx, right = False, labels = False)
        #nz = pd.cut(xm["zu"], bins = gridx, right = False, labels = False)
        # nx = pd.cut(self.wrapped_monomers["xu"], bins = gridx, right = False, labels = False)
        # ny = pd.cut(self.wrapped_monomers["yu"], bins = gridy, right = False, labels = False)
        # nz = pd.cut(self.wrapped_monomers["zu"], bins = gridz, right = False, labels = False)

        # nx = (xm[["xu"]]/ actual_cell_length[["x"]].to_numpy())
        # ny = (xm[["yu"]]/ actual_cell_length[["y"]].to_numpy())
        # nz = (xm[["zu"]]/ actual_cell_length[["z"]].to_numpy())
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


    def bond_bond_correlation(self, show_plot = False):
        #For PVA-100: 7200 x99 x 99 upper triangular matrix, axis-0: j, axis-1: j + 1
        bond_bond_correlation_array = np.zeros([self.atom_coords.no_polymers, self.atom_coords.polymer_length-1, self.atom_coords.polymer_length-1])
        for i in range(0, self.atom_coords.no_polymers):
        #for i in range(0,1):
            # subset = self.atom_coords.bond_vectors[(self.atom_coords.bond_vectors["mol_id"] == i + 1)].to_numpy()[:, 1:4]
            # #Normalise bond vectors 
            # subset = subset/np.linalg.norm(subset, axis = 1, keepdims = True)

            start = i * self.atom_coords.polymer_length + 1              # 1, L+1, 2L+1, ...
            end   = (i + 1) * self.atom_coords.polymer_length         # L-1, 2L-1, 3L-1, ...

            subset_df = self.atom_coords.bond_vectors.loc[start:end]
            subset = subset_df.to_numpy()[:, 1:4]
            subset = subset / np.linalg.norm(subset, axis=1, keepdims=True)
            #print(subset)
            bond_bond_corr_array_polymer = np.zeros([self.atom_coords.polymer_length-1, self.atom_coords.polymer_length-1])
            box_algos_lib.bond_bond_correlation(subset, bond_bond_corr_array_polymer, subset.shape[0], subset.shape[1])

            bond_bond_correlation_array[i,:,:] = bond_bond_corr_array_polymer
        bond_correlation_average = np.mean(bond_bond_correlation_array, axis = 0)
        # print(bond_correlation_average)
        #np.savetxt("debug_bond_bond_average.txt", bond_correlation_average)
        #bond_correlation_average = np.loadtxt("debug_bond_bond_average.txt")
        diag_means = np.zeros(self.atom_coords.polymer_length-1)
        for i in range(0, self.atom_coords.polymer_length -1):
        #     #print(np.mean(np.diagonal(bond_correlation_average, offset = i)))
            diag_means[i] = np.mean(np.diagonal(bond_correlation_average, offset = i))
        self.results.bond_bond_correlation = diag_means
        self.atom_coords.positions = np.arange(1, 100)
        # plt.scatter(positions, (diag_means))
        # plt.xlabel("n")
        # plt.ylabel(r"cos\theta(n)")
        # plt.title("Distribution of bond-bond correlations")
        #plt.show()

    def bond_distribution(self):
        #self.atom_coords.get_nematic_vector_5()
        counts = (
            self.atom_coords.bond_vectors
            .groupby(["xid", "yid", "zid"])
            .size()                      # count rows per combination
            .reset_index(name="count")   # put counts in a column
        )

        dist = counts["count"].value_counts().sort_index()
        print(dist)

    def atom_distribution(self):
        wrapped_coords = self.atom_coords.wrap_coordinates(self.atom_coords.datapd.iloc[:, 1:4])
        #wrapped_coords = self.atom_coords.datapd.iloc[:, 1:4]
        print(wrapped_coords)

        coords = wrapped_coords[['xu', 'yu', 'zu']].to_numpy()

        # Define 33 bins between min and max for each axis
        ranges = [
            (wrapped_coords['xu'].min(), wrapped_coords['xu'].max()),
            (wrapped_coords['yu'].min(), wrapped_coords['yu'].max()),
            (wrapped_coords['zu'].min(), wrapped_coords['zu'].max()),
        ]

        # 3D histogram: 33 bins along each of x, y, z
        H, edges = np.histogramdd(
            coords,
            bins=(33, 33, 33),
            range=ranges,
        )

        print("Histogram shape:", H.shape)   # (33, 33, 33)
        print("Total count in histogram:", int(H.sum()))




        H_1d = H.ravel()
        unique_counts, freqs = np.unique(H_1d, return_counts=True)

        for c, f in zip(unique_counts, freqs):
            print(f"count={int(c)} in {int(f)} cells")



    def get_density_dist_per_box(self, nridges = 33):
        """Uses new method with combinations to calculate local density (i.e. density per box)"""
        data = self.atom_coords.datapd
        data_indexed = data.reset_index().rename({"index": "atom_id"}).set_index(["xid", "yid", "zid"]).sort_index()
        #print(data_indexed)
        data_grouped = data_indexed.groupby(level=["xid", "yid", "zid"])
        no_points_per_box = data_grouped.size()
        print(no_points_per_box)
        local_densities = no_points_per_box / self.atom_coords.local_volume
        print(data.shape[0]/self.atom_coords.volume)
        #no_points_per_box.to_csv("local_densities.txt")
        self.results.avg_local_density = np.mean(local_densities)
        return local_densities

    def get_density_dist(self, nridges = None):
        if nridges is None:
            nridges = self.atom_coords.nridges
        _, _, boxes = self.atom_coords.assign_center_of_mass(nridges = nridges)
        boxes = boxes.set_index(["xid", "yid", "zid"]).sort_index()
        data_grouped = boxes.groupby(level=["xid", "yid", "zid"])
        no_points_per_box = data_grouped.size()
        return no_points_per_box/self.atom_coords.local_volume


    def get_kuhn_length(self):
        """See eq. 20 and 21 in Saras paper 2018"""
        pp_avg_bond_length = np.zeros(self.atom_coords.no_polymers)
        N_minus_1 = self.atom_coords.polymer_length - 1
        for i in tqdm(range(0, self.atom_coords.no_polymers)):
        #for i in range(0, 5):
            # Calculate |r_{i+1}-r_i|
            subset = self.atom_coords.datapd.iloc[(i*(self.atom_coords.polymer_length)):((i+1)*(self.atom_coords.polymer_length)), 1:4]
            #print(i*(self.atom_coords.polymer_length), (i+1)*(self.atom_coords.polymer_length)-1)
            #print(subset)
            point_to_point_dist = np.diff(subset, axis = 0)
            #print(point_to_point_dist)
            abs_point_to_point_dist = np.sqrt(point_to_point_dist[:, 0]**2 + point_to_point_dist[:, 1]**2 + point_to_point_dist[:, 2]**2)
            #print(abs_point_to_point_dist)
            #print(abs_point_to_point_dist)
            #pp_avg_bond_length[i] = 1/(N_minus_1) * ((np.sum(abs_point_to_point_dist)))
            pp_avg_bond_length[i] = abs_point_to_point_dist.mean()
            #print(pp_avg_bond_length[i])
            #print(point_to_point_dist)
        R_e2 = self.end_to_end_distance()
        N_e = R_e2/(N_minus_1*np.mean(pp_avg_bond_length))
        print(np.mean(N_e))
        return N_e


    def get_entanglement_length(self, bond_cutoff: int = 175):
        #print(self.atom_coords.datapd)
        #print(self.atom_coords.datapd.columns)
        poly = self.atom_coords.datapd.groupby(["mol_id"])
        self.atom_coords.ppa = self.atom_coords.datapd.copy()
        self.atom_coords.ppa = self.atom_coords.ppa.drop(columns = ["xu", "yu", "zu"])
        self.atom_coords.ppa["bond_angle"] = np.nan
        self.atom_coords.ppa["ppa_index"] = np.nan
        self.atom_coords.ppa["ppa_length"] = np.nan
        #print(poly)
        ppa_index = 1
        ppa_length = 0
        #for mol_id, block in poly:
        mean_ppa_list = np.zeros(self.atom_coords.no_polymers); j = 0;
        for mol_id, block in tqdm(poly, total=poly.ngroups, desc="Processing mol_id groups"):
            atom_ids = block.index.to_list()

            # bonds between consecutive atoms in this monomer
            block_bv = self.atom_coords.bond_vectors.loc[atom_ids[:-1]]
            # PPA info for *internal* atoms (where angles are defined)
            block_ppa = self.atom_coords.ppa.loc[atom_ids[1:-1]].copy()

            # number of angles = number of internal atoms
            n_angles = len(block_ppa)          # == len(block_bv) - 1

            ppa_start_loc = 0
            ppa_length = 0
            # assume these exist as columns in block_ppa
            col_angle = block_ppa.columns.get_loc("bond_angle")
            col_index = block_ppa.columns.get_loc("ppa_index")
            col_length = block_ppa.columns.get_loc("ppa_length")

            for i in range(n_angles):          # i = 0 .. n_angles-1
                r1 = -1 * block_bv.iloc[i, :3].to_numpy()
                r2 =      block_bv.iloc[i+1, :3].to_numpy()
                bond_angle = calc_bond_angle(r1, r2)

                block_ppa.iat[i, col_angle] = bond_angle
                block_ppa.iat[i, col_index] = ppa_index
                ppa_length += 1

                if bond_angle > bond_cutoff:
                    # close current PPA stretch up to and including i
                    block_ppa.iloc[ppa_start_loc:(i+1), col_length] = ppa_length
                    ppa_start_loc = i + 1
                    ppa_length = 0
                    ppa_index += 1

                # last angle in this block: close final stretch
                elif i == n_angles - 1:
                    block_ppa.iloc[ppa_start_loc:(i+1), col_length] = ppa_length
                    ppa_start_loc = i + 1
                    ppa_length = 0
                    ppa_index += 1

            #print(block_ppa)
            #print(block_ppa["ppa_length"].mean())
            mean_ppa_list[j] = block_ppa["ppa_length"].mean()
            j = j + 1
        return mean_ppa_list


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
            print("total number clusters w/ >= 2 elements: %i" %(self.results.total_number_clusters))
            print("total number independent crystalline domains: %i" %(self.results.total_number_independent_clusters))
            print("average cluster size crystalline domains: %f" %(self.results.mean_cluster_size))
            print("total number crystalline/all grid elements: %i/%i -> cryt_frac = %f" 
                %(self.results.total_number_crystalline_grid_elements,total_box_elements, 
                    self.results.total_number_crystalline_grid_elements/total_box_elements))
            print("earlier calculated frac_cryst = %f" %(self.frac_cryst))
        return label_matrix

    def bin_label_matrix(self):#bins: np.array):
        """Function to bin the cluster size into a histogram"""
        label_matrix = self.merge_boxes_2()
        result = cluster_pdfs_from_label_matrix(label_matrix, vol_per_site= self.atom_coords.volume/self.atom_coords.n_atoms)
        unique_values, counts = np.unique(label_matrix, return_counts=True)
        mask = unique_values > 0
        cluster_sizes = counts[mask]
        labelcount = len(cluster_sizes)   
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

        ax1.bar(result["bin_mid"], result["pdf1"]/labelcount,
                width=result["bin_high"] - result["bin_low"],
                align="center", edgecolor="k", alpha=0.7)
        ax1.set_xscale("log"); ax1.set_yscale("log")
        ax1.set_xlabel("Cluster size s"); ax1.set_ylabel("  [amount of clusters]")
        ax1.set_title("Cluster size distribution")

        ax2.bar(result["bin_mid"], result["pdf2"],
                width=result["bin_high"] - result["bin_low"],
                align="center", edgecolor="k", color="orange", alpha=0.7)
        ax2.set_xscale("log"); ax2.set_yscale("log")
        ax2.set_xlabel("Cluster size s"); ax2.set_ylabel("  [volume fraction]")
        ax2.set_title("Volume-weighted distribution ")
        plt.suptitle(r"Cluster distributuion, PVA-%i, $\phi = 0.442$" %(self.atom_coords.polymer_length))
        plt.tight_layout()
        plt.show()

    #hist, bin_edges = np.histogram(counts[1:], bins)


def make_ncluster_from_labels(unique_values, counts):
    """
    Build ncluster[size] = number of clusters of that size.
    Assumes label 0 is background and is ignored.
    """
    mask = unique_values > 0
    cluster_sizes = counts[mask]

    imax1 = int(cluster_sizes.max())
    jmax  = imax1 + 1000

    ncluster = np.zeros(jmax + 1, dtype=int)
    size_hist = np.bincount(cluster_sizes)
    ncluster[:len(size_hist)] = size_hist

    return ncluster, imax1, jmax


def bin_ncluster(ncluster, imax1, jmax):
    """
    Apply the non-uniform binning from the C++ protocol.
    Returns:
      bin_low   : lower size (inclusive) of each non-empty bin
      bin_high  : upper size (inclusive) of each non-empty bin
      bin_count : total number of clusters in that bin
    """
    bin_low, bin_high, bin_count = [], [], []

    def add_bin(low, high):
        if low > imax1:
            return
        high_eff = min(high, imax1)
        total = ncluster[low:high_eff + 1].sum()
        if total == 0:
            return                  # skip empty bins (mirrors C++ if bincluster[i]>0)
        bin_low.append(low)
        bin_high.append(high_eff)
        bin_count.append(total)

    # width-1 bins: sizes 1–5
    for s in range(1, 6):
        add_bin(s, s)

    # width-2 bin: 6–7
    add_bin(6, 7)

    # width-3 bin: 8–10
    add_bin(8, 10)

    # width-10 bins: 11–20, 21–30, ..., 61–70
    for high in range(20, 71, 10):
        add_bin(high - 9, high)

    # width-30 bin: 71–100
    add_bin(71, 100)

    # width-100 bins: 101–200, ..., 601–700
    for high in range(200, 701, 100):
        add_bin(high - 99, high)

    # width-300 bin: 701–1000
    add_bin(701, 1000)

    # width-1000 bins: 1001–2000, 2001–3000, ...
    for high in range(2000, jmax, 1000):
        add_bin(high - 999, high)

    return np.array(bin_low), np.array(bin_high), np.array(bin_count)


def compute_pdfs(bin_low, bin_high, bin_count, nemcount, labelcount):
    """
    Compute pdf1 and pdf2 for each bin, mirroring the C++ code.

    pdf1[i] = bin_count / bin_width
            = number of clusters per unit cluster size
            → the cluster NUMBER distribution (normalised for bin width)

    pdf2[i] = bin_midpoint * pdf1[i] / nemcount
            = (n1+n2)/2 * pdf1 / nemcount
            → the VOLUME-WEIGHTED distribution per unit cluster size,
              normalised by total number of particles in clusters

    Parameters
    ----------
    bin_low, bin_high : arrays of bin edges (inclusive)
    bin_count         : number of clusters in each bin
    nemcount          : total number of particles/sites belonging to clusters
                        (= sum over all cluster sizes s of s * ncluster[s])

    Returns
    -------
    bin_mid  : bin midpoint (n1+n2)/2  — use as x-axis
    pdf1     : number PDF per unit cluster size
    pdf2     : volume-weighted PDF per unit cluster size, normalised by nemcount
    """
    bin_width = bin_high - bin_low + 1          # width of each bin in size units
    bin_mid   = (bin_low + bin_high) / 2.0      # midpoint, same as (n1+n2)/2 in C++

    # C++: pdf1[i] = bincluster[i] / (n2 - n1)
    # Note: C++ uses n2-n1 where n1 is the PREVIOUS bin's upper edge,
    # which equals bin_width here (since bins are contiguous).
    pdf1 = bin_count / (bin_width.astype(float))

    # C++: pdf2[i] = (n1+n2) * pdf1[i] * 0.5 / nemcount
    pdf2 = bin_mid * pdf1 / nemcount

    return bin_mid, pdf1, pdf2


def averages_from_ncluster(ncluster, imax1, vol, nemcount, labelcount):
    """
    Compute average cluster size and average volume directly from ncluster,
    mirroring the C++ averages loop.

    C++:
      avencluster += ncluster[i] * i        → mean cluster size
      avevolume   += ncluster[i] * i * i    → second moment → mean volume

    Parameters
    ----------
    vol        : volume per particle/site
    nemcount   : total number of particles in clusters
    labelcount : total number of clusters
    """
    sizes = np.arange(1, imax1 + 1)
    nc    = ncluster[1:imax1 + 1]

    avencluster = np.sum(nc * sizes) / labelcount
    avevolume   = np.sum(nc * sizes * sizes) * vol / nemcount

    return avencluster, avevolume


# ── Full pipeline ──────────────────────────────────────────────────────────────

def cluster_pdfs_from_label_matrix(label_matrix, vol_per_site=1.0):
    """
    Full pipeline from a 3D label_matrix to pdf1 and pdf2.

    Parameters
    ----------
    label_matrix : 3D numpy array of integer cluster labels (0 = background)
    vol_per_site : volume of a single lattice site (default 1)

    Returns
    -------
    dict with keys:
      bin_low, bin_high, bin_mid, bin_count
      pdf1        : number PDF per unit cluster size
      pdf2        : volume-weighted PDF per unit cluster size / nemcount
      avencluster : mean cluster size
      avevolume   : mean cluster volume
    """
    unique_values, counts = np.unique(label_matrix, return_counts=True)

    # Build ncluster histogram
    ncluster, imax1, jmax = make_ncluster_from_labels(unique_values, counts)

    # Derived counts
    mask       = unique_values > 0
    cluster_sizes = counts[mask]
    labelcount = len(cluster_sizes)               # total number of clusters
    nemcount   = int(cluster_sizes.sum())         # total particles in clusters

    # Binning
    bin_low, bin_high, bin_count = bin_ncluster(ncluster, imax1, jmax)

    # PDFs
    bin_mid, pdf1, pdf2 = compute_pdfs(bin_low, bin_high, bin_count, nemcount, labelcount)

    # Averages
    avencluster, avevolume = averages_from_ncluster(
        ncluster, imax1, vol_per_site, nemcount, labelcount
    )

    print(f"Total clusters (labelcount) : {labelcount}")
    print(f"Total particles in clusters : {nemcount}")
    print(f"Largest cluster size        : {imax1}")
    print(f"Average cluster size        : {avencluster:.4f}")
    print(f"Average cluster volume      : {avevolume:.4f}")

    return {
        "bin_low"    : bin_low,
        "bin_high"   : bin_high,
        "bin_mid"    : bin_mid,
        "bin_count"  : bin_count,
        "pdf1"       : pdf1,
        "pdf2"       : pdf2,
        "avencluster": avencluster,
        "avevolume"  : avevolume,
    }
    unique_values, counts = np.unique(label_matrix, return_counts=True)
    
    ncluster, imax1, jmax = make_ncluster_from_labels(unique_values, counts)
    bin_low, bin_high, bin_count = bin_ncluster(ncluster, imax1, jmax)
    
    return {
        "bin_low": bin_low,      # lower size bound of each bin (inclusive)
        "bin_high": bin_high,    # upper size bound of each bin (inclusive)
        "bin_count": bin_count,  # number of clusters in that bin
        "imax1": imax1,
        "jmax": jmax,
    }



def check_labels(label_matrix, nridges, i,j,k):
    """Helper functuion to check if label_matrix is equal to one of its neighbours. If not, print current value of label_matrix and all neighbours"""
    if label_matrix[i,j,k] != label_matrix[(i -1)% nridges,j,k] or label_matrix[i,j,k] != label_matrix[(i+1)%nridges, j, k]:
        if label_matrix[i,j,k] != label_matrix[i,(j-1)% nridges,k] or label_matrix[i,j,k] != label_matrix[i, (j+1) % nridges, k]:
            if label_matrix[i,j,k] != label_matrix[i,j,(k-1) % nridges] or label_matrix[i,j,k] != label_matrix[i, j, (k+1) % nridges]:
                print("position %i %i %i, current label %i \nneighbours in -x %i and +x %i" %(i,j,k,label_matrix[i,j,k], label_matrix[(i -1)% nridges, j,k], label_matrix[(i +1)% nridges, j,k]))
                print("neighbours in -y %i, +y %i" %(label_matrix[i, (j - 1) %nridges, k], label_matrix[i, (j + 1) %nridges, k]))
                print("neighbours in -z direction %i, +z %i \n" %(label_matrix[i, j, (k - 1)%nridges], label_matrix[i, j, (k + 1)%nridges]))


    # if label_matrix[i,j,k] != label_matrix[(i-1)% nridges,j,k]:
    #     if label_matrix[i,j,k] != label_matrix[(i+1)%nridges, j, k]:  
    #         if label_matrix[i,j,k] != label_matrix[i,(j-1)% nridges,k]: 
    #             if label_matrix[i,j,k] != label_matrix[i, (j+1) % nridges, k]:
    #                 if label_matrix[i,j,k] != label_matrix[i,j,(k-1) % nridges]:
    #                     if label_matrix[i,j,k] != label_matrix[i, j, (k+1) % nridges]:

class results(object):
    def __init__(self):
        return

    # Include plotting functions here 


def plot_hk_matrix_2d(polymer, ndot_cutoff=0.98, cryst_threshold=0.8):
    # 3D array of cluster labels
    data = polymer.merge_boxes_2(ndot_cutoff = ndot_cutoff, print_results=True)
    Nx, Ny, Nz = data.shape

    dfc = polymer.df_cryst.copy()

    # Integer indices for the grid
    ix = dfc["xid"].astype(int).to_numpy()
    iy = dfc["yid"].astype(int).to_numpy()
    iz = dfc["zid"].astype(int).to_numpy()

    # ---- Build cryst_bool, x_ev, y_ev, z_ev arrays aligned with data ----
    cryst_vals = np.full((Nx, Ny, Nz), np.nan, dtype=float)
    xev_vals   = np.full((Nx, Ny, Nz), np.nan, dtype=float)
    yev_vals   = np.full((Nx, Ny, Nz), np.nan, dtype=float)
    zev_vals   = np.full((Nx, Ny, Nz), np.nan, dtype=float)

    cryst_vals[ix, iy, iz] = dfc["cryst_bool"].to_numpy()
    xev_vals[ix, iy, iz]   = dfc["x_ev"].to_numpy()
    yev_vals[ix, iy, iz]   = dfc["y_ev"].to_numpy()
    zev_vals[ix, iy, iz]   = dfc["z_ev"].to_numpy()

    # ---- Global colour scale for cluster labels ----
    global_max = int(data.max())
    n_labels = max(global_max, 0)

    base_colors = plt.cm.viridis(np.linspace(0, 1, n_labels + 1))
    base_colors[0] = np.array([1.0, 0.0, 0.0, 1.0])  # RGBA
    cmap = ListedColormap(base_colors)
    bounds = np.arange(-0.5, n_labels + 1.5, 1)
    norm = BoundaryNorm(bounds, cmap.N)

    for k in range(Nz):
        slice_k = data[:, :, k]
        cryst_k = cryst_vals[:, :, k]

        cell_size = 0.8
        fig_width = max(4, Ny * cell_size)
        fig_height = max(4, Nx * cell_size)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        im = ax.imshow(slice_k, cmap=cmap, norm=norm, origin='lower')

        ax.set_title(
            f"Cluster labels, PVA-100, T = 0.88, Tdot = 10**(-3), z = {k}"
        )
        ax.set_xlabel("y")
        ax.set_ylabel("x")
        ax.set_xlim(-0.5, Ny - 0.5)
        ax.set_ylim(-0.5, Nx - 0.5)

        # --- Cell labels: cluster + cryst_bool, plus z-edge values above/below ---
        label_only_nonzero = False
        step_x = 1
        step_y = 1

        for x in range(0, Nx, step_x):
            for y in range(0, Ny, step_y):
                val = slice_k[x, y]
                if label_only_nonzero and val == 0:
                    continue

                cryst_val = cryst_k[x, y]

                # Central cell text: cluster + cryst_bool
                lines = [f"{int(val)}"]
                if not np.isnan(cryst_val):
                    lines.append(f"{cryst_val:.3f}")
                text_str = "\n".join(lines)

                ax.text(
                    y, x, text_str,
                    ha='center', va='center',
                    color='white' if val < global_max / 2 else 'black',
                    fontsize=5,
                )

                # If this cell is "red" (cryst_bool < threshold), skip z-val_edge
                if np.isnan(cryst_val) or cryst_val < cryst_threshold:
                    continue

                # Periodic neighbors in z
                k_above = (k + 1) % Nz   # border "above" current slice
                k_below = (k - 1) % Nz   # border "below" current slice

                # z-edge above: between (x,y,k) and (x,y,k_above)
                val_above = None
                if (not np.isnan(xev_vals[x, y, k]) and
                    not np.isnan(xev_vals[x, y, k_above]) and
                    not np.isnan(yev_vals[x, y, k]) and
                    not np.isnan(yev_vals[x, y, k_above]) and
                    not np.isnan(zev_vals[x, y, k]) and
                    not np.isnan(zev_vals[x, y, k_above])):

                    P2Q_above = (
                        (xev_vals[x, y, k] * xev_vals[x, y, k_above]) +
                        (yev_vals[x, y, k] * yev_vals[x, y, k_above]) +
                        (zev_vals[x, y, k] * zev_vals[x, y, k_above])
                    )
                    val_above = 1.5 * P2Q_above**2 - 0.5

                # z-edge below: between (x,y,k_below) and (x,y,k)
                val_below = None
                if (not np.isnan(xev_vals[x, y, k_below]) and
                    not np.isnan(xev_vals[x, y, k]) and
                    not np.isnan(yev_vals[x, y, k_below]) and
                    not np.isnan(yev_vals[x, y, k]) and
                    not np.isnan(zev_vals[x, y, k_below]) and
                    not np.isnan(zev_vals[x, y, k])):

                    P2Q_below = (
                        (xev_vals[x, y, k_below] * xev_vals[x, y, k]) +
                        (yev_vals[x, y, k_below] * yev_vals[x, y, k]) +
                        (zev_vals[x, y, k_below] * zev_vals[x, y, k])
                    )
                    val_below = 1.5 * P2Q_below**2 - 0.5

                # Z-edge values: separate texts, colored by ndot_cutoff, with z+/z- prefix
                if val_above is not None:
                    color_above = 'blue' if val_above >= ndot_cutoff else 'red'
                    ax.text(
                        y, x + 0.35, f"z+ {val_above:.2f}",
                        ha='center', va='center',
                        color=color_above,
                        fontsize=4,
                    )

                if val_below is not None:
                    color_below = 'blue' if val_below >= ndot_cutoff else 'red'
                    ax.text(
                        y, x - 0.35, f"z- {val_below:.2f}",
                        ha='center', va='center',
                        color=color_below,
                        fontsize=4,
                    )

        # --- Edge labels in x,y directions with periodic BC, colored by ndot_cutoff ---

        edge_step_x = 1
        edge_step_y = 1

        # "Horizontal" borders (between rows): (x,y)–(x_next,y), x_next periodic
        for x in range(0, Nx, edge_step_x):
            x_next = (x + 1) % Nx
            for y in range(0, Ny, edge_step_y):
                c1 = cryst_vals[x, y, k]
                # only require THIS cell not red
                if np.isnan(c1) or c1 < cryst_threshold:
                    continue

                if (np.isnan(xev_vals[x,     y, k]) or np.isnan(xev_vals[x_next, y, k]) or
                    np.isnan(yev_vals[x,     y, k]) or np.isnan(yev_vals[x_next, y, k]) or
                    np.isnan(zev_vals[x,     y, k]) or np.isnan(zev_vals[x_next, y, k])):
                    continue

                P2Q = (
                    (xev_vals[x,     y, k] * xev_vals[x_next, y, k]) +
                    (yev_vals[x,     y, k] * yev_vals[x_next, y, k]) +
                    (zev_vals[x,     y, k] * zev_vals[x_next, y, k])
                )
                val_edge = 1.5 * P2Q**2 - 0.5
                edge_color = 'blue' if val_edge >= ndot_cutoff else 'red'

                # Midpoint in x with special handling for wrap edge (Nx-1 ↔ 0)
                if x == Nx - 1 and x_next == 0:
                    mid_x = Nx - 0.5  # bottom boundary
                else:
                    mid_x = x + 0.5

                mid_y = y
                ax.text(
                    mid_y, mid_x, f"{val_edge:.2f}",
                    ha='center', va='center',
                    color=edge_color,
                    fontsize=4,
                )

        # "Vertical" borders (between columns): (x,y)–(x,y_next), y_next periodic
        for x in range(0, Nx, edge_step_x):
            for y in range(0, Ny, edge_step_y):
                y_next = (y + 1) % Ny

                c1 = cryst_vals[x, y, k]
                if np.isnan(c1) or c1 < cryst_threshold:
                    continue

                if (np.isnan(xev_vals[x, y,     k]) or np.isnan(xev_vals[x, y_next, k]) or
                    np.isnan(yev_vals[x, y,     k]) or np.isnan(yev_vals[x, y_next, k]) or
                    np.isnan(zev_vals[x, y,     k]) or np.isnan(zev_vals[x, y_next, k])):
                    continue

                P2Q = (
                    (xev_vals[x, y,     k] * xev_vals[x, y_next, k]) +
                    (yev_vals[x, y,     k] * yev_vals[x, y_next, k]) +
                    (zev_vals[x, y,     k] * zev_vals[x, y_next, k])
                )
                val_edge = 1.5 * P2Q**2 - 0.5
                edge_color = 'blue' if val_edge >= ndot_cutoff else 'red'

                # Midpoint in y with special handling for wrap edge (Ny-1 ↔ 0)
                if y == Ny - 1 and y_next == 0:
                    mid_y = Ny - 0.5  # right boundary
                else:
                    mid_y = y + 0.5

                mid_x = x
                ax.text(
                    mid_y, mid_x, f"{val_edge:.2f}",
                    ha='center', va='center',
                    color=edge_color,
                    fontsize=4,
                )

        cbar = fig.colorbar(im, ax=ax, label='cluster label')
        cbar.set_ticks(np.arange(0, n_labels + 1))
        fig.tight_layout()
        plt.savefig(f"hk_debug/PVA_1000_polynumber_20_T088_zlayer_{k}.pdf")
        plt.close()


def plot_hk_matrix(polymer):
    label_matrix = polymer.merge_boxes()


    Nx, Ny, Nz = label_matrix.shape

    # Create coordinate grids for x, y, z (index positions)
    x_idx, y_idx, z_idx = np.indices((Nx, Ny, Nz))

    # Flatten to 1D
    x = x_idx.ravel()
    y = y_idx.ravel()
    z = z_idx.ravel()
    l = label_matrix.ravel()      # values

    mask = l > 0
    x = x[mask]
    y = y[mask]
    z = z[mask]
    l = l[mask]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(x, y, z, c=l, cmap='viridis', s=10)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    fig.colorbar(sc, ax=ax, label='l')

    plt.show()



def main():
    """Testing function"""

    #first_timestep_e5 = polymer("../../data/pva-100/quick_quench/equil_t_085_tdot_e-3_run2_goodrun_time_0.txt")
    slow_quench_e5_time_100e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
    # slow_quench_e5_time_95e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_9500000.txt")
    # slow_quench_e5_time_90e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_9000000.txt")
    # slow_quench_e5_time_85e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_8500000.txt")
    # slow_quench_e5_time_80e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_8000000.txt")


    slow_quench_e5_time_100e5.atom_coords.get_nematic_vector_5("10e5_T05_last_timestep_boxes_ev.txt")
    # slow_quench_e5_time_95e5.atom_coords.get_nematic_vector_5("10e5_T05_timestep_boxes_ev_time_95e5.txt")
    # slow_quench_e5_time_90e5.atom_coords.get_nematic_vector_5("10e5_T05_timestep_boxes_ev_time_90e5.txt")
    # slow_quench_e5_time_85e5.atom_coords.get_nematic_vector_5("10e5_T05_timestep_boxes_ev_time_85e5.txt")
    # slow_quench_e5_time_80e5.atom_coords.get_nematic_vector_5("10e5_T05_timestep_boxes_ev_time_80e5.txt")
    slow_quench_e5_time_100e5.read_cryst("hk_debug/10e5_T05_last_timestep_boxes_ev.txt")
    # slow_quench_e5_time_95e5.read_cryst("hk_debug/10e5_T05_timestep_boxes_ev_time_95e5.txt")
    # slow_quench_e5_time_90e5.read_cryst("hk_debug/10e5_T05_timestep_boxes_ev_time_90e5.txt")
    # slow_quench_e5_time_85e5.read_cryst("hk_debug/10e5_T05_timestep_boxes_ev_time_85e5.txt")
    # slow_quench_e5_time_80e5.read_cryst("hk_debug/10e5_T05_timestep_boxes_ev_time_80e5.txt")
    #print(first_timestep_e5.atom_coords.combinations)
    #last_timestep_e5.get_density_dist()
    label_matrix = slow_quench_e5_time_100e5.merge_boxes()
    # slow_quench_e5_time_95e5.merge_boxes()
    # slow_quench_e5_time_90e5.merge_boxes()
    # slow_quench_e5_time_85e5.merge_boxes()
    # slow_quench_e5_time_80e5.merge_boxes()
    bins = np.array([1,2,3,4,5,7,10,20,30,40,50,60,70])
    bin_label_matrix(label_matrix, bins)


    # plot_hk_matrix_2d(slow_quench_e5_time_100e5)
    
    
if __name__ == "__main__":
    main()