from analyse7 import polymer, atom_coords
from clibraries.boxAlgorithmsInC import box_algos_lib, hoshen_kopelman_lib

import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # may be needed for some MPL versions
import os 
from matplotlib.colors import ListedColormap, BoundaryNorm


def return_cell(df, x,y,z):
    return df[(df["xid"] == x) & (df["yid"] == y) & (df["zid"] == z)]


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import numpy as np

def plot_3d_with_nematic_cell_center(atoms_df, nematic_df, xid, yid, zid,
                                     rel_length=0.5, abs_scale=None):
    """
    3D scatter of atoms at (xu, yu, zu) colored by mol_id,
    plus nematic vector (x_ev, y_ev, z_ev) for the given cell (xid, yid, zid),
    starting at the cell center.

    Arrow length control:
    - If abs_scale is given: arrow length = abs_scale (in xu/yu/zu units).
    - Else: arrow length = rel_length * (cell bounding-box diagonal).

    Parameters
    ----------
    atoms_df : pandas.DataFrame
        Columns: 'xu', 'yu', 'zu', 'mol_id', 'xid', 'yid', 'zid'.
    nematic_df : pandas.DataFrame
        Columns: 'xid', 'yid', 'zid', 'x_ev', 'y_ev', 'z_ev'.
    xid, yid, zid : float or int
        Cell indices to select.
    rel_length : float
        Fraction of cell diagonal to use as arrow length (if abs_scale is None).
    abs_scale : float or None
        If not None, use this absolute length for the nematic vector.
    """

    # Atoms in this cell
    cell_atoms = atoms_df[
        (atoms_df['xid'] == xid) &
        (atoms_df['yid'] == yid) &
        (atoms_df['zid'] == zid)
    ]
    if cell_atoms.empty:
        print("No atoms in cell (xid, yid, zid) =", (xid, yid, zid))
        return

    # Nematic data for this cell
    cell_nem = nematic_df[
        (nematic_df['xid'] == xid) &
        (nematic_df['yid'] == yid) &
        (nematic_df['zid'] == zid)
    ]
    if cell_nem.empty:
        print("No nematic data for cell (xid, yid, zid) =", (xid, yid, zid))
        return

    nem = cell_nem.iloc[0]
    vx, vy, vz = nem['x_ev'], nem['y_ev'], nem['z_ev']

    # Normalize eigenvector (just in case)
    vnorm = np.sqrt(vx**2 + vy**2 + vz**2)
    if vnorm > 0:
        vx, vy, vz = vx / vnorm, vy / vnorm, vz / vnorm

    # Cell center (bounding-box center)
    xmin, xmax = cell_atoms['xu'].min(), cell_atoms['xu'].max()
    ymin, ymax = cell_atoms['yu'].min(), cell_atoms['yu'].max()
    zmin, zmax = cell_atoms['zu'].min(), cell_atoms['zu'].max()

    ox = 0.5 * (xmin + xmax)
    oy = 0.5 * (ymin + ymax)
    oz = 0.5 * (zmin + zmax)

    # Cell diagonal length
    cell_diag = np.sqrt((xmax - xmin)**2 + (ymax - ymin)**2 + (zmax - zmin)**2)

    # Arrow length
    if abs_scale is not None:
        L = abs_scale
    else:
        L = rel_length * cell_diag

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Color atoms by mol_id
    mol_cat = cell_atoms['mol_id'].astype('category')
    mol_codes = mol_cat.cat.codes

    sc = ax.scatter(
        cell_atoms['xu'], cell_atoms['yu'], cell_atoms['zu'],
        c=mol_codes,
        cmap='tab20',
        s=15,
        alpha=0.9
    )

    # Nematic vector as a short arrow
    ax.quiver(
        ox, oy, oz,
        vx * L, vy * L, vz * L,
        color='black',
        linewidth=2,
        arrow_length_ratio=0.15
    )

    # Optional line
    ax.plot(
        [ox, ox + vx * L],
        [oy, oy + vy * L],
        [oz, oz + vz * L],
        color='black',
        linestyle='--',
        linewidth=1
    )

    cbar = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.1)
    cbar.set_label('mol_id')
    cbar.set_ticks(range(len(mol_cat.cat.categories)))
    cbar.set_ticklabels(mol_cat.cat.categories)

    ax.set_xlabel('xu')
    ax.set_ylabel('yu')
    ax.set_zlabel('zu')
    ax.set_title(r'Positions & nematic vector, $\dot{T}=10^{-5}, T = 0.5$, cell %i,%i,%i' %(xid,yid,zid))

    plt.tight_layout()
    plt.show()

def plot_3d_with_nematic_cell_center_bond_vector(
    atoms_df,
    nematic_df,
    bond_df,
    xid,
    yid,
    zid,
    rel_length=0.5,
    abs_scale=None,
    bond_scale=1.0,
):
    """
    3D scatter of atoms at (xu, yu, zu) colored by atom_id,
    plus:
      - nematic vector (x_ev, y_ev, z_ev) from cell center
      - bond vector from each atom position.

    Arrow length for nematic:
      - If abs_scale is given: arrow length = abs_scale.
      - Else: rel_length * (cell bounding-box diagonal).

    Bond vectors:
      - Drawn from each atom position using (bondx, bondy, bondz) * bond_scale.
    """

    # ---- Select atoms in this cell ----
    cell_atoms = atoms_df[
        (atoms_df['xid'] == xid) &
        (atoms_df['yid'] == yid) &
        (atoms_df['zid'] == zid)
    ]
    if cell_atoms.empty:
        print("No atoms in cell (xid, yid, zid) =", (xid, yid, zid))
        return

    # ---- Nematic data for this cell ----
    cell_nem = nematic_df[
        (nematic_df['xid'] == xid) &
        (nematic_df['yid'] == yid) &
        (nematic_df['zid'] == zid)
    ]
    if cell_nem.empty:
        print("No nematic data for cell (xid, yid, zid) =", (xid, yid, zid))
        return

    nem = cell_nem.iloc[0]
    vx, vy, vz = nem['x_ev'], nem['y_ev'], nem['z_ev']

    # Normalize nematic eigenvector
    vnorm = np.sqrt(vx**2 + vy**2 + vz**2)
    if vnorm > 0:
        vx, vy, vz = vx / vnorm, vy / vnorm, vz / vnorm

    # ---- Cell center from bounding box of atoms ----
    xmin, xmax = cell_atoms['xu'].min(), cell_atoms['xu'].max()
    ymin, ymax = cell_atoms['yu'].min(), cell_atoms['yu'].max()
    zmin, zmax = cell_atoms['zu'].min(), cell_atoms['zu'].max()

    ox = 0.5 * (xmin + xmax)
    oy = 0.5 * (ymin + ymax)
    oz = 0.5 * (zmin + zmax)

    # Cell diagonal and nematic arrow length
    cell_diag = np.sqrt((xmax - xmin)**2 + (ymax - ymin)**2 + (zmax - zmin)**2)
    if abs_scale is not None:
        L = abs_scale
    else:
        L = rel_length * cell_diag

    # ---- Prepare bond vectors for atoms in this cell ----
    # Join bond_df on atom_id index with cell_atoms
    # This keeps only atoms that have bond vectors
    cell_with_bonds = cell_atoms.join(
        bond_df[['bondx', 'bondy', 'bondz']],
        how='inner'  # only those atom_ids present in bond_df
    )

    # ---- Plotting ----
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Color atoms by atom_id (index)
    atom_ids = cell_atoms.index.astype('category')
    atom_codes = atom_ids.codes

    sc = ax.scatter(
        cell_atoms['xu'], cell_atoms['yu'], cell_atoms['zu'],
        c=atom_codes,
        cmap='tab20',
        s=15,
        alpha=0.9
    )

    # Nematic vector from cell center
    ax.quiver(
        ox, oy, oz,
        vx * L, vy * L, vz * L,
        color='black',
        linewidth=2,
        arrow_length_ratio=0.15
    )
    ax.plot(
        [ox, ox + vx * L],
        [oy, oy + vy * L],
        [oz, oz + vz * L],
        color='black',
        linestyle='--',
        linewidth=1
    )

    # Bond vectors from each atom position (if available)
    if not cell_with_bonds.empty:
        ax.quiver(
            cell_with_bonds['xu'],
            cell_with_bonds['yu'],
            cell_with_bonds['zu'],
            cell_with_bonds['bondx'] * bond_scale,
            cell_with_bonds['bondy'] * bond_scale,
            cell_with_bonds['bondz'] * bond_scale,
            color='red',
            linewidth=1,
            arrow_length_ratio=0.2
        )

    # Colorbar labeled by atom_id
    cbar = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.1)
    cbar.set_label('atom_id')
    # Be careful: many atoms -> crowded; for small cells it's OK
    cbar.set_ticks(range(len(atom_ids.categories)))
    cbar.set_ticklabels(atom_ids.categories)

    ax.set_xlabel('xu')
    ax.set_ylabel('yu')
    ax.set_zlabel('zu')
    ax.set_title(r'Positions & nematic vector, $\dot{T}=10^{-5}, T = 0.5$, cell %i,%i,%i' %(xid,yid,zid))

    plt.tight_layout()
    plt.show()


def load_ncolor_3d(path):
    # Load all numeric rows; comment lines starting with '#' are skipped
    data2d = np.loadtxt(path, comments='#')
    # data2d.shape = ((nz+1)*(ny+1), nx+1)

    # Count how many z-slices there are from the headers
    with open(path) as f:
        nz_plus_1 = sum(
            1 for line in f
            if line.lstrip().startswith("# z =")
        )

    total_rows, nx_plus_1 = data2d.shape
    ny_plus_1 = total_rows // nz_plus_1

    # Reshape to [z, y, x]
    arr_zyx = data2d.reshape(nz_plus_1, ny_plus_1, nx_plus_1)

    # Transpose to [x, y, z]
    arr_xyz = np.transpose(arr_zyx, (2, 1, 0))

    # Make sure it's integer type
    return arr_xyz.astype(int)

def plot_hk_matrix_2d_from_array(path,
    title_prefix="PVA-100, T = 0.88, clusters via nematic20cc",
    output_prefix="hk_debug/pva-100_T088_nematic20cc_npolymer20",
    label_cells=True,
):
    """
    Plot z-slices of a 3D integer array 'data' (shape: Nx, Ny, Nz).

    Parameters
    ----------
    data : np.ndarray
        3D array of cluster labels, shape (Nx, Ny, Nz).
    title_prefix : str
        Text prefix for the plot title (before ', z = k').
    output_prefix : str
        Prefix for output PDF filenames. Files will be:
        f"{output_prefix}_{k}.pdf"
    label_cells : bool
        If True, write the integer cluster label in each cell.
    """
    # Ensure output directory exists
    data = load_ncolor_3d(path)
    data = data.astype(float)
    data[data > 35937] = np.nan
    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    # Shape of the data
    Nx, Ny, Nz = data.shape

    # ---- Global colour scale for cluster labels ----
    global_max = int(np.nanmax(data))
    print(global_max)
    n_labels = max(global_max, 0)

    # Build colormap: label 0 is red, the rest use viridis
    base_colors = plt.cm.viridis(np.linspace(0, 1, n_labels + 1))
    base_colors[0] = np.array([1.0, 0.0, 0.0, 1.0])  # label 0 = red
    cmap = ListedColormap(base_colors)
    cmap.set_bad(color="black")  # how NaNs will appear
    bounds = np.arange(-0.5, n_labels + 1.5, 1)
    norm = BoundaryNorm(bounds, cmap.N)
    for k in range(Nz):
        slice_k = data[:, :, k]

        cell_size = 0.8
        fig_width = max(4, Ny * cell_size)
        fig_height = max(4, Nx * cell_size)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        im = ax.imshow(slice_k, cmap=cmap, norm=norm, origin='lower')

        ax.set_title(f"{title_prefix}, z = {k}")
        ax.set_xlabel("y")
        ax.set_ylabel("x")
        ax.set_xlim(-0.5, Ny - 0.5)
        ax.set_ylim(-0.5, Nx - 0.5)

        # --- Optional cell labels: just the cluster ID ---
        if label_cells:
            for x in range(Nx):
                for y in range(Ny):
                    print(x,y)
                    val = slice_k[x, y]
                    ax.text(
                        y, x, f"{int(val)}",
                        ha='center', va='center',
                        color='white' if val < global_max / 2 else 'black',
                        fontsize=5,
                    )
        print("test")
        # Colorbar
        # cbar = fig.colorbar(im, ax=ax, label='cluster label')

        # max_ticks = 20
        # if n_labels <= max_ticks:
        #     cbar.set_ticks(np.arange(0, n_labels + 1))
        # elif n_labels > 0:
        #     step = max(1, n_labels // max_ticks)
        #     ticks = np.arange(0, n_labels + 1, step)
        #     cbar.set_ticks(ticks)
        fig.tight_layout()
        out_name = f"{output_prefix}_{k}.pdf"
        plt.savefig(out_name)
        plt.close()
        print(f"Saved {out_name}")


def main():
    # last_timestep_e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
    # last_timestep_e5.atom_coords.get_nematic_vector_5()
    # #print(last_timestep_e5.atom_coords.bond_vectors)
    # x,y,z = 6,11,25
    # #x,y,z = 6,6,25
    # cryst_row = return_cell(last_timestep_e5.atom_coords.df_cryst, x,y,z)
    # single_cell = return_cell(last_timestep_e5.atom_coords.datapd, x,y,z)
    # bond_with_cells = last_timestep_e5.atom_coords.bond_vectors.join(last_timestep_e5.atom_coords.datapd[['xid','yid','zid']], how = "right")
    # single_bonds = return_cell(bond_with_cells,x,y,z)
    # print(single_cell, single_bonds)
    # plot_3d_with_nematic_cell_center_bond_vector(single_cell, cryst_row,single_bonds,  x,y,z)
    path = "clibraries/nematic20/pva-100_comparison/equil_t_088_tdot_e-3_time_24000000.txt_ncolourprint.txt"
    plot_hk_matrix_2d_from_array(path)
    

if __name__ == "__main__":
    main()