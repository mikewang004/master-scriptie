import pandas as pd


def hk_in_python(df_cryst, ndot_cutoff = 0.97, nridges = 33, cryst_cutoff = 0.8):
    # Set (xid, yid, zid) as an index for fast lookup
    df_cryst = df_cryst.set_index(["xid", "yid", "zid"])

    # Iterate over all rows
    for (xid, yid, zid), row in df_cryst.iterrows():
        # Define neighbor coordinates
        neighbours = {
            "x_minus": (apply_pbc(xid, nridges= nridges), yid, zid),
            "y_minus": (xid, apply_pbc(yid, nridges= nridges), zid),
            "z_minus": (xid, yid,apply_pbc(zid, nridges= nridges)),
        }

        # Retrieve neighbors if they exist (skip boundaries)
        neighbours_rows = {}
        for name, coord in neighbours.items():
            try:
                neighbours_rows[name] = df_cryst.loc[coord]
            except KeyError:
                neighbours_rows[name] = None  # no neighbor (e.g. xid == 0)



def apply_pbc(x, pbc = False, nridges = 33):
    if pbc == False: 
        return x - 1 
    return (x - 1 + nridges) % nridges