import numpy as np 
import scipy as sp 
import pandas as pd

class coarse_grained_hist():

    def __init__(self, label_matrix, total_volume, local_volume, nemcount):
        self.label_matrix = label_matrix
        self.local_volume = local_volume
        self.total_volume = total_volume
        self.nemcount = nemcount
        self.unique_values, self.counts = np.unique(label_matrix, return_counts = True)



    def build_bin_edges(self):
        """
        Reproduces the manual bin boundaries from the C++ code.
        Returns an array of bin edges (left-exclusive, right-inclusive).
        """
        edges = []
        edges += list(range(0, 6))          # bins [0,1], [1,2], ..., [4,5]  → width 1
        edges += [7]                         # [5, 7]                         → width 2
        edges += [10]                        # [7, 10]                        → width 3
        edges += list(range(20, 80, 10))     # [10,20], [20,30], ..., [60,70] → width 10
        edges += [100]                       # [70, 100]                      → width 30
        edges += list(range(200, 800, 100))  # [100,200], ..., [600,700]      → width 100
        edges += [1000]                      # [700, 1000]                    → width 300
        edges += list(range(2000, 10001, 1000))  # [1000,2000], ...           → width 1000
        return np.array(edges)


    def coarse_bin(self, cluster_sizes: np.ndarray):
        """
        Parameters
        ----------
        cluster_sizes : 1D array of cluster sizes (one entry per cluster)
        nemcount      : total number of particles
        vol           : volume per particle (cell volume)

        Returns
        -------
        bin_centers : midpoint of each non-empty bin
        pdf1        : count density (clusters per unit size)
        pdf2        : volume-weighted PDF (fraction of particles per unit size)
        bin_counts  : raw cluster counts per bin
        """
        edges = self.build_bin_edges()

        # np.histogram does all the aggregation in one call
        bin_counts, _ = np.histogram(cluster_sizes, bins=edges)

        # Keep only non-empty bins
        mask        = bin_counts > 0
        counts      = bin_counts[mask]
        left_edges  = edges[:-1][mask]
        right_edges = edges[1:][mask]

        widths      = right_edges - left_edges          # bin widths
        bin_centers = 0.5 * (left_edges + right_edges)  # bin midpoints

        pdf1 = counts / widths                          # count density
        pdf2 = bin_centers * pdf1 / self.nemcount            # volume-weighted, normalised

        return bin_centers, pdf1, pdf2, counts



    def build_hist(self, crystallinity, save_string = None):
        centers, pdf1, pdf2, counts = self.coarse_bin(self.counts)
        df = pd.DataFrame({
            "clustersize": centers * 2,          # bin label ≈ right edge (i in C++)
            "volume":      centers * 2 * self.local_volume,
            "vol_Vtot":    centers * 2 * self.local_volume / self.total_volume,
            "Ncluster":    counts,
            "pdf1":        pdf1 / len(self.unique_values),    # labelcount = total clusters
            "volume_pdf":  pdf2,
            "volume_pdf_x_cryst": pdf2 * crystallinity,
        })
        if save_string == None:
            print(df)
        else:
            df.to_csv("%s" %save_string, index=False, sep=" ")

        return df