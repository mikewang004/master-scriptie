import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
from experiment2 import quench_PVA, calculate_ppa_from_file
from analyse8 import polymer



N_values = np.array([50, 200, 300, 500, 1000])
Ne_values_quench_start = np.array([31.139, 42.884, 44.490, 46.659, 48.111])
Ne_values_quench_end = np.array([30.788, 42.392, 44.232, 46.324, 47.898])


plt.scatter(N_values, Ne_values_quench_start, label = "T = 1.0")
plt.scatter(N_values, Ne_values_quench_end, label = "T = 0.88")
plt.title("Ne vs N")
plt.ylabel("Ne")
plt.xlabel("N")
plt.savefig("Ne_vs_N.pdf")
plt.show()