from analyse7 import polymer, atom_coords
from clibraries.boxAlgorithmsInC import box_algos_lib, hoshen_kopelman_lib

import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # may be needed for some MPL versions



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

def main():
    last_timestep_e5 = polymer("../../data/pva-100/cooling_tdot_e-5_time_10000000.txt")
    last_timestep_e5.atom_coords.get_nematic_vector_5()
    #print(last_timestep_e5.atom_coords.bond_vectors)
    x,y,z = 6,11,25
    #x,y,z = 6,6,25
    cryst_row = return_cell(last_timestep_e5.atom_coords.df_cryst, x,y,z)
    single_cell = return_cell(last_timestep_e5.atom_coords.datapd, x,y,z)
    bond_with_cells = last_timestep_e5.atom_coords.bond_vectors.join(last_timestep_e5.atom_coords.datapd[['xid','yid','zid']], how = "right")
    single_bonds = return_cell(bond_with_cells,x,y,z)
    print(single_cell, single_bonds)
    plot_3d_with_nematic_cell_center_bond_vector(single_cell, cryst_row,single_bonds,  x,y,z)

if __name__ == "__main__":
    main()