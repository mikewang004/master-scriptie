import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import pandas as pd
from numba import jit, njit
import matplotlib.cm as cm
import matplotlib.colors as colors

def wrap_coordinates(column, minlength, total_length):
    return (column - minlength) % total_length + minlength

@njit
def norm_keep_dims_1d(A):
    n, m = A.shape
    result = np.empty((n, 1), dtype = A.dtype)
    for i in range(n):
        s = 0.0
        for j in range(m):
            s += A[i, j] * A[i, j]
        result[i, 0] = np.sqrt(s)
    return result


def nematic_vector_loop_2(bond_vectors: pd.DataFrame) -> pd.DataFrame:
    """
    Uses bond vector cell alignment to calculate crystallinity per cell.
    Expects columns: bx, by, bz, nx, ny, nz, index = atom_id.
    """
    # bring atom_id out of the index (if needed)
    bond_vectors_indexed = bond_vectors.reset_index()  # atom_id becomes a column

    # group by grid cell
    g = bond_vectors_indexed.groupby(["nx", "ny", "nz"])

    # apply nematic tensor calculation per cell
    df_cryst = g.apply(calc_nematic_tensor_pandas)

    # tidy up
    df_cryst = df_cryst.reset_index()
    df_cryst = df_cryst.sort_values(by=["nz", "nx", "ny"]).reset_index(drop=True)
    return df_cryst

def calc_nematic_tensor_pandas(block: pd.DataFrame) -> pd.Series:
    # use bx, by, bz from the block (one cell)
    array = block[["bx", "by", "bz"]].to_numpy()
    array_length = array.shape[0]

    if array_length == 0:
        return pd.Series(
            {"cryst_bool": 0.0, "x_ev": 0.0, "y_ev": 0.0, "z_ev": 0.0}
        )

    # normalize each bond vector; norm_keep_dims_1d should handle row-wise norms
    array = array / norm_keep_dims_1d(array)

    # nematic tensor Q
    outer = (array.T @ array) / array_length
    Q = 1.5 * outer - 0.5 * np.eye(3)

    # eigenvalues / eigenvectors
    labda, ev = np.linalg.eigh(Q)
    idx = np.argmax(labda)
    max_labda = labda[idx]
    max_ev = ev[:, idx]

    # enforce +z convention
    if max_ev[2] < 0:
        max_ev = -max_ev

    return pd.Series(
        {"cryst_bool": max_labda, "x_ev": max_ev[0], "y_ev": max_ev[1], "z_ev": max_ev[2]}
    )

def plot_cryst_bool(df):
    df_plot = df.copy()
    df_plot['nx'] = df_plot['nx'].astype(int)
    df_plot['ny'] = df_plot['ny'].astype(int)
    df_plot['nz'] = df_plot['nz'].astype(int)

    # Global normalization for viridis, based only on values > 0.8
    mask_good = df_plot['cryst_bool'] > 0.8
    if mask_good.any():
        vmin = df_plot.loc[mask_good, 'cryst_bool'].min()
        vmax = df_plot.loc[mask_good, 'cryst_bool'].max()
    else:
        # fallback in case everything is <= 0.8
        vmin, vmax = 0.8, 1.0

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap('viridis')

    for z in sorted(df_plot['nz'].unique()):
        df_z = df_plot[df_plot['nz'] == z]

        # Pivot to 2D grid: rows = ny, columns = nx
        grid = df_z.pivot_table(
            index='ny',
            columns='nx',
            values='cryst_bool',
            aggfunc='mean'  # in case of duplicates
        )

        # Sort axes to get consistent layout
        grid = grid.sort_index(axis=0).sort_index(axis=1)

        y_vals = grid.index.to_numpy()
        x_vals = grid.columns.to_numpy()
        Z = grid.to_numpy()

        # Build RGBA image using viridis for >0.8, red for <=0.8
        rgba = cmap(norm(Z))  # shape (ny, nx, 4)

        # NaNs should be transparent
        nan_mask = np.isnan(Z)
        rgba[nan_mask, 3] = 0.0

        # Threshold mask: red for cryst_bool <= 0.8
        low_mask = (Z <= 0.8) & ~nan_mask
        rgba[low_mask] = np.array([1.0, 0.0, 0.0, 1.0])  # pure red, fully opaque

        # Plot
        fig, ax = plt.subplots(figsize=(6, 6))

        # Define extent so that cell centers are at integer (nx, ny)
        x_min, x_max = x_vals.min(), x_vals.max()
        y_min, y_max = y_vals.min(), y_vals.max()
        extent = [x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5]

        im = ax.imshow(
            rgba,
            origin='lower',
            extent=extent,
            aspect='equal'
        )

        # Annotate each cell with cryst_bool
        for iy, y in enumerate(y_vals):
            for ix, x in enumerate(x_vals):
                val = Z[iy, ix]
                if not np.isnan(val):
                    ax.text(
                        x, y,
                        f"{val:.2f}",
                        ha='center',
                        va='center',
                        color='black',
                        fontsize=6
                    )

        # Axes/ticks
        ax.set_xlabel("nx")
        ax.set_ylabel("ny")
        ax.set_title(f"Crystallinity, nz = {z}")

        ax.set_xticks(x_vals)
        ax.set_yticks(y_vals)
        ax.set_xticklabels(x_vals)
        ax.set_yticklabels(y_vals)

        # Optional: grid lines
        ax.set_xticks(x_vals + 0.5, minor=True)
        ax.set_yticks(y_vals + 0.5, minor=True)
        ax.grid(which='minor', color='k', linestyle='-', linewidth=0.3)
        ax.tick_params(which='minor', length=0)

        # Save one PDF per z level
        outfile = f"crystallinity_z{z}.pdf"
        plt.savefig(outfile, bbox_inches='tight')
        plt.close(fig)



class atom_coords():

    def __init__(self, path_to_file, nridges = 33):
        self.datapd = self.read_in_file(path_to_file)
        self.volume, self.boxlengths, self.dimensions = self.get_volume_box(path_to_file)
        self.wrapped_monomers = self.wrap_coordinates_all_data()
        self.bond_vectors = self.calculate_bond_vectors()
        self.bond_vectors = self.make_cell_grid()
        self.df_cryst = self.calc_cryst()


    def read_in_file(self, filename, ovito_single_z = False):
        if ovito_single_z == True:
            skiprows = 6
        else:
            skiprows = 9
        datapd = pd.read_csv(filename, sep = " ", header = None, skiprows = skiprows)
        datapd.columns = ["atom_id", "mol_id", "xu", "yu", "zu"]
        datapd = datapd.sort_values("atom_id")
        datapd = datapd.set_index("atom_id")

        return datapd

    def get_volume_box(self, path_to_file, ovito_single_z = False):
        """Calculates volume of the total (!) simulation box."""
        if ovito_single_z == True:
            skiprows = 2
        else:
            skiprows = 5
        datapd_first_rows = (pd.read_csv(path_to_file, sep = " ", header = None, skiprows = skiprows, nrows = 3))
        datapd_first_rows = datapd_first_rows.rename({0: "min", 1: "max"}, axis = 1).rename({0: "x", 1: "y", 2:"z"}, axis = 0)
        length = datapd_first_rows.iloc[:, 1] - datapd_first_rows.iloc[:, 0]
        volume = length.loc["x"]*length.loc["y"]*length.loc["z"]

        return volume, length, datapd_first_rows


    def wrap_coordinates_all_data(self):
        """Converts coordinates to wrapped coordinates. Input can be float or array (float)"""
        wrapped_coords = self.datapd
        wrapped_coords["xu"] = wrap_coordinates(self.datapd["xu"], self.dimensions.loc["x", "min"], self.boxlengths["x"])
        wrapped_coords["yu"] = wrap_coordinates(self.datapd["yu"], self.dimensions.loc["y", "min"], self.boxlengths["y"])
        wrapped_coords["zu"] = wrap_coordinates(self.datapd["zu"], self.dimensions.loc["z", "min"], self.boxlengths["z"])
        return wrapped_coords

    def calculate_bond_vectors(self):
        print(self.datapd)

        shifted = self.datapd.groupby('mol_id')[['xu', 'yu', 'zu']].shift(-1)
        bond_vecs = shifted - self.datapd[['xu', 'yu', 'zu']]

        bond_vecs.columns = ['bx', 'by', 'bz']
        bond_vecs = bond_vecs.dropna()
        return bond_vecs


    def make_cell_grid(self, cell_length = 2):

        nridges = (self.boxlengths/cell_length).astype(int)
        actual_cell_length = self.boxlengths/nridges
        print(self.dimensions)
        print(nridges)
        gridx = np.linspace(self.dimensions.loc["x", "min"], self.dimensions.loc["x", "max"], nridges.loc["x"]+1)
        gridy = np.linspace(self.dimensions.loc["y", "min"], self.dimensions.loc["y", "max"], nridges.loc["y"]+1)
        gridz = np.linspace(self.dimensions.loc["z", "min"], self.dimensions.loc["z", "max"], nridges.loc["z"]+1)
        print(gridx)
        #Shift monomers to start from 0 
        #print(self.bond_vectors)
        nx = pd.cut(self.wrapped_monomers["xu"] + self.bond_vectors["bx"]/2, bins = gridx, right = False, labels = False)
        ny = pd.cut(self.wrapped_monomers["yu"] + self.bond_vectors["by"]/2, bins = gridy, right = False, labels = False)
        nz = pd.cut(self.wrapped_monomers["zu"] + self.bond_vectors["bz"]/2, bins = gridz, right = False, labels = False)
        self.bond_vectors = self.bond_vectors.assign(
            nx=nx,
            ny=ny,
            nz=nz,
        )
        print(self.bond_vectors)

        return self.bond_vectors


    def calc_cryst(self):
        df_cryst = nematic_vector_loop_2(self.bond_vectors)
        print(df_cryst.loc)
        return df_cryst


def main():
    filename = "equil_t_088_tdot_e-3_time_24000000.txt"
    pva_100_tnumber_20 = atom_coords(filename)
    df_cryst = pva_100_tnumber_20.df_cryst
    plot_cryst_bool(df_cryst)
    #print(pva_100_tnumber_20.datapd)

if __name__ == "__main__":
    main()