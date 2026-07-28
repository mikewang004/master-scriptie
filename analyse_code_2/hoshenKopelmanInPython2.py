import numpy as np
import pandas as pd


# ── helpers ──────────────────────────────────────────────────────────────────

def find(labels, x):
    """Path-compressed find for the union-find structure."""
    while labels[x] != x:
        labels[x] = labels[labels[x]]   # path compression (halving)
        x = labels[x]
    return x


def union(labels, a, b):
    """Union by root; returns the new common root."""
    ra, rb = find(labels, a), find(labels, b)
    if ra == rb:
        return ra
    # merge the larger root into the smaller index (keeps roots small)
    if ra < rb:
        labels[rb] = ra
        return ra
    else:
        labels[ra] = rb
        return rb


# ── main function ─────────────────────────────────────────────────────────────

def hoshen_kopelman_domains(
    df: pd.DataFrame,
    s_threshold: float = 0.8,
    dot_threshold: float = 0.97,
    label_matrix = None
) -> np.ndarray:
    """
    Cluster crystalline lattice cells into nematic domains using the
    Hoshen-Kopelman algorithm with periodic boundary conditions.

    Parameters
    ----------
    df            : DataFrame with columns nx, ny, nz, cryst_bool (S value),
                    x_ev, y_ev, z_ev.
    s_threshold   : minimum S value for a cell to be considered crystalline.
    dot_threshold : minimum |n1·n2| for two cells to be considered aligned.

    Returns
    -------
    domain_labels : int32 ndarray of shape (Nx, Ny, Nz).
                    0  → not crystalline / no domain
                    k  → domain index k  (k ≥ 1)
    """

    # ── grid dimensions ───────────────────────────────────────────────────────
    Nx = df["nx"].max() + 1
    Ny = df["ny"].max() + 1
    Nz = df["nz"].max() + 1
    N  = Nx * Ny * Nz

    # ── build lookup arrays indexed by flat index ─────────────────────────────
    cryst = np.zeros(N, dtype=bool)
    evecs = np.zeros((N, 3), dtype=float)

    for _, row in df.iterrows():
        idx = int(row["nx"]) * (Ny * Nz) + int(row["ny"]) * Nz + int(row["nz"])
        cryst[idx] = row["cryst_bool"] > s_threshold
        evecs[idx] = [row["x_ev"], row["y_ev"], row["z_ev"]]

    # ── union-find initialisation ─────────────────────────────────────────────
    # labels[i] == i means i is its own root; 0 is reserved for "no domain"
    # We use indices 1..N for cells (shift by 1 so 0 stays free).
    uf = np.arange(N + 1, dtype=np.intp)   # uf[0] unused

    def cell_id(ix, iy, iz):
        """Flat 1-based cell id."""
        return ix * (Ny * Nz) + iy * Nz + iz + 1   # +1: 1-based

    def flat0(ix, iy, iz):
        """0-based flat index into cryst / evecs arrays."""
        return ix * (Ny * Nz) + iy * Nz + iz

    def aligned(i0, j0):
        """True if both crystalline and eigenvectors are aligned."""
        if not (cryst[i0] and cryst[j0]):
            return False
        dot = abs(np.dot(evecs[i0], evecs[j0]))
        return dot >= dot_threshold

    # ── first pass: raster scan with look-back neighbours (+ PBC wrap) ───────
    # For each cell we check the three "previous" neighbours along x, y, z
    # (with periodic wrap-around handled explicitly).
    for ix in range(Nx):
        for iy in range(Ny):
            for iz in range(Nz):
                i0 = flat0(ix, iy, iz)
                if not cryst[i0]:
                    continue

                cid = cell_id(ix, iy, iz)

                # neighbours: -x, -y, -z  (with PBC)
                neighbours = [
                    (ix - 1) % Nx, iy,           iz,
                    ix,           (iy - 1) % Ny,  iz,
                    ix,            iy,            (iz - 1) % Nz,
                ]
                for k in range(3):
                    nx_, ny_, nz_ = neighbours[3*k], neighbours[3*k+1], neighbours[3*k+2]
                    j0 = flat0(nx_, ny_, nz_)
                    if aligned(i0, j0):
                        union(uf, cid, cell_id(nx_, ny_, nz_))

    # ── second pass: flatten all roots ───────────────────────────────────────
    for i in range(1, N + 1):
        find(uf, i)   # side-effect: path compression

    # ── relabel roots to consecutive integers starting from 1 ────────────────
    root_map: dict[int, int] = {}
    label_counter = 0

    domain_labels = np.zeros((Nx, Ny, Nz), dtype=np.int32)

    for ix in range(Nx):
        for iy in range(Ny):
            for iz in range(Nz):
                i0 = flat0(ix, iy, iz)
                if not cryst[i0]:
                    domain_labels[ix, iy, iz] = 0
                    continue
                root = find(uf, cell_id(ix, iy, iz))
                if root not in root_map:
                    label_counter += 1
                    root_map[root] = label_counter
                domain_labels[ix, iy, iz] = root_map[root]

    return domain_labels


# ── example usage ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Replace with your actual df_cryst
    path = "%s/equil_t_088_tdot_e-3_cryst_time_12000000.txt" %("../data_online/PVA-100/icryst_T088_Tdot_e-3/crystallisation")
    df_cryst = pd.read_csv(path, sep = " ")
    labels = hoshen_kopelman_domains(df_cryst)

    Nx = df_cryst["nx"].max() + 1
    Ny = df_cryst["ny"].max() + 1
    Nz = df_cryst["nz"].max() + 1

    print(f"Grid shape : {labels.shape}")
    print(f"Unique domains (excl. 0): {np.unique(labels[labels > 0])}")
    print(f"Number of domains        : {labels.max()}")
    print(f"Crystalline cells        : {(labels > 0).sum()} / {Nx*Ny*Nz}")
