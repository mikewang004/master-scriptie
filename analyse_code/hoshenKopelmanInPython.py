import pandas as pd
import numpy as np

class unionFind:
    """
    Union-Find (Disjoint Set Union) with the same semantics as the C version:

    - labels[x] is an alias for label x.
    - labels[0] stores the highest label already used.
    - Labels for actual sets start at 1.
    """

    def __init__(self, max_labels):
        """
        Initialize the data structure.

        max_labels: size of the labels array (indices 0 .. max_labels-1)
        """
        self.n_labels = max_labels
        # equivalent to calloc: zero-initialized
        self.labels = [0] * self.n_labels
        self.labels[0] = 0  # no labels used yet

    def find(self, x):
        """
        Return the canonical label (root) for the equivalence class containing x.
        Includes path compression.
        """
        y = x
        # follow chain to root
        while self.labels[y] != y:
            y = self.labels[y]

        # path compression
        while self.labels[x] != x:
            z = self.labels[x]
            self.labels[x] = y
            x = z

        return y

    def union(self, x, y):
        """
        Join the sets containing x and y.
        Return the canonical label of the resulting class.
        """
        root_y = self.find(y)
        self.labels[self.find(x)] = root_y
        return root_y

    def make_set(self):
        """
        Create a new equivalence class and return its label.
        """
        self.labels[0] += 1
        assert self.labels[0] < self.n_labels, "Out of label space"
        self.labels[self.labels[0]] = self.labels[0]
        return self.labels[0]


def hk_in_python(df_cryst, nridges, ndot_cutoff = 0.97, cryst_cutoff = 0.8):
    # Set (xid, yid, zid) as an index for fast lookup
    df_cryst = df_cryst.set_index(["xid", "yid", "zid"]).sort_index()
    nridgex = nridges[["x"]].item(); nridgey = nridges[["y"]].item(); nridgez = nridges[["z"]].item()
    label_array = np.zeros((nridgex, nridgey, nridgez), dtype = int)
    max_labels = int(33**5)
    new_labels = np.zeros(int(nridgex * nridgey * nridgez), dtype = int)
    hk = unionFind(max_labels)
    print(ndot_cutoff)
    # Iterate over all rows
    i = 0
    for m in range(0,2):
        if m == 0: 
            apply_pbc_bool = False
        else:
            apply_pbc_bool = True
        for (xid, yid, zid), current_row in df_cryst.iterrows():

            # Only attempt if larger than cryst cutoff
            if current_row.loc["cryst_bool"] > cryst_cutoff:

                xm = apply_pbc(xid, pbc = apply_pbc_bool, nridges = nridgex)
                ym = apply_pbc(yid, pbc = apply_pbc_bool, nridges = nridgey)
                zm = apply_pbc(zid, pbc = apply_pbc_bool, nridges = nridgez)
            # Define neighbor coordinates
                x_minus = df_cryst.loc[(xm, yid, zid)] if (xm, yid, zid) in df_cryst.index else None
                y_minus = df_cryst.loc[(xid, ym, zid)] if (xid, ym, zid) in df_cryst.index else None
                z_minus = df_cryst.loc[(xid, yid, zm)] if (xid, yid, zm) in df_cryst.index else None
            
            # See if neighbouring coordinates can be merged
                xmin_check = check_ndot(current_row, x_minus, ndot_cutoff= ndot_cutoff, cryst_cutoff= cryst_cutoff)
                ymin_check = check_ndot(current_row, y_minus, ndot_cutoff= ndot_cutoff, cryst_cutoff= cryst_cutoff)
                zmin_check = check_ndot(current_row, z_minus, ndot_cutoff= ndot_cutoff, cryst_cutoff= cryst_cutoff)

                check_total = xmin_check + ymin_check + zmin_check
                ix, iy, iz = int(xid), int(yid), int(zid)   # ensure ints

                if check_total == 0: #Start new cluster
                    #print(ix, iy, iz)
                    label_array[ix, iy, iz] = hk.make_set()
                    #print(label_array[ix, iy, iz])
                # elif check_total == 1: #Merge with x/y/z neighbour
                #     if xmin_check == 1:
                #         label_array[ix, iy, iz] = label_array[xm, iy, iz]
                #     elif ymin_check == 1:
                #         label_array[ix, iy, iz] = label_array[ix, ym, iz]
                #     else: #zmin_check == 1
                #         label_array[ix, iy, iz] = label_array[ix, iy, zm]
                # elif check_total == 2: #Merge two clusters 
                #     if xmin_check == 1:
                #         if ymin_check == 1:
                #             label_array[ix, iy, iz] = hk.union(label_array[xm, iy, iz], label_array[ix, ym, iz])
                #         else: #zmin_check ==1 
                #             label_array[ix, iy, iz] = hk.union(label_array[xm, iy, iz], label_array[ix, iy, zm])
                #     else: #ymin_check == 1 and zmin_check == 1
                #         label_array[ix, iy, iz] = hk.union(label_array[ix, ym, iz], label_array[ix, iy, zm])
                # else: #check_total == 3
                #     label_array[ix, iy, iz] = hk.union(label_array[xm, iy, iz], label_array[ix, ym, iz])
                #     label_array[ix, iy, iz] = hk.union(label_array[xm, iy, iz], label_array[ix, iy, zm])
                #     label_array[ix, iy, iz] = hk.union(label_array[ix, ym, iz], label_array[ix, iy, zm])

                elif check_total == 1:
                    if xmin_check == 1:
                        label_array[ix, iy, iz] = hk.find(label_array[xm, iy, iz]) 
                    elif ymin_check == 1:
                        label_array[ix, iy, iz] = hk.find(label_array[ix, ym, iz])  
                    else:
                        label_array[ix, iy, iz] = hk.find(label_array[ix, iy, zm])  

                elif check_total == 2:
                    if xmin_check == 1:
                        if ymin_check == 1:
                            label_array[ix, iy, iz] = hk.union(
                                hk.find(label_array[xm, iy, iz]),
                                hk.find(label_array[ix, ym, iz])
                            )
                        else:
                            label_array[ix, iy, iz] = hk.union(
                                hk.find(label_array[xm, iy, iz]),
                                hk.find(label_array[ix, iy, zm])
                            )
                    else:
                        label_array[ix, iy, iz] = hk.union(
                            hk.find(label_array[ix, ym, iz]),
                            hk.find(label_array[ix, iy, zm])
                        )
                else:  # check_total == 3
                    r1 = hk.find(label_array[xm, iy, iz])
                    r2 = hk.find(label_array[ix, ym, iz])
                    r3 = hk.find(label_array[ix, iy, zm])
                    root = hk.union(r1, r2)
                    root = hk.union(root, r3)
                    label_array[ix, iy, iz] = root
    #np.savetxt("hk_old_labels.txt", hk.labels)

    for x in range(0, nridgex):
        for y in range(0, nridgey):
            for z in range(0, nridgez):
                lab = (label_array[x, y, z])
                if lab == 0:
                    continue  # ignore background/unlabeled
                root = hk.find(lab)
                # optional compaction of labels:
                if new_labels[root] == 0:
                    new_labels[0] += 1
                    new_labels[root] = new_labels[0]
                label_array[x, y, z] = new_labels[root]
    #np.savetxt("hk_new_labels.txt", new_labels)

    # with open("labels_all_sheets.txt", "w") as f:
    #     nx, ny, nz = label_array.shape
    #     for z in range(nz):
    #         f.write(f"# z = {z}\n")
    #         np.savetxt(f, label_array[:, :, z], fmt="%d")
    #         f.write("\n")
    return label_array

        


# def apply_pbc(x, pbc = False, nridges = 33):
#     if pbc == False: 
#         return int(x - 1)
#     return int((x - 1 + nridges) % nridges)

def apply_pbc(x, pbc=False, nridges=33):
    """
    x is already 0-based (0..nridges-1).
    Returns the backward neighbor index (0-based).
    """
    if not pbc:
        return int(x) - 1              # plain backward: may be -1 for x=0 (→ None lookup)
    return (int(x) - 1) % nridges      # PBC backward: wraps 0 → nridges-1

def check_ndot(current_row, neighbour_row, cryst_cutoff = 0.8, ndot_cutoff = 0.97):
    """Returns 1 if 1.5 * (n_i dot n_j) - 0.5 >= ndot_cutoff, otherwise 0 """
    if neighbour_row is None:
        return 0;
    elif neighbour_row.loc["cryst_bool"] < cryst_cutoff:
        return 0;

    inner_product = (current_row.loc["x_ev"] * neighbour_row.loc["x_ev"] + current_row.loc["y_ev"] * neighbour_row.loc["y_ev"] + \
        current_row.loc["z_ev"] * neighbour_row.loc["z_ev"]) 
    ndot_check = 1.5 * inner_product* inner_product - 0.5
    # print(current_row)
    # print(neighbour_row)
    # print(inner_product)
    # print(ndot_check)
    # print()

    #ndot_check = inner_product
    if ndot_check >= ndot_cutoff: 
        return 1; 
    return 0;