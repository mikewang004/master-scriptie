import numpy as np 
import scipy as sp 
import pandas as pd

class coarse_grained_hist():

    def __init__(self, label_matrix, total_volume,local_volume, n_atoms, crystallinity):
        self.label_matrix = label_matrix
        self.local_volume = local_volume
        self.total_volume = total_volume
        self.n_atoms = n_atoms
        self.crystallinity = crystallinity
        self.unique_values, self.counts = np.unique(label_matrix, return_counts = True)



#TODO look at how pdf2/volume pdf is calculated

    def build_bin_edges(self):
        labels = (
            list(range(1, 8))           # 1,2,3,4,5,6,7
            + [10]                      # 10
            + list(range(20, 71, 10))   # 20,30,40,50,60,70
            + [100]                     # 100
            + list(range(200, 701, 100))# 200,300,400,500,600,700
            + [1000]                    # 1000
            + list(range(2000, 10001, 1000))  # 2000,3000,...,10000
        )

        # edges: prepend 0 so the first bin is [0, 1) → size 1 only
        # each bin is [labels[i-1], labels[i]) → contains sizes labels[i-1]..labels[i]-1
        edges = [0] + labels

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
        pdf2 = bin_centers * pdf1 / (np.sum(self.counts[1:]))            # volume-weighted, normalised

        return bin_centers, pdf1, pdf2, counts



    def build_hist(self, crystallinity, save_string = None):
        mask = self.unique_values > 0
        cluster_sizes = self.counts[mask]
        centers, pdf1, pdf2, counts = self.coarse_bin(cluster_sizes)
        df = pd.DataFrame({
            "clustersize": centers * 2,          
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