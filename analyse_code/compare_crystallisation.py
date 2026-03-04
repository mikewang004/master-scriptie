import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 


times = np.array([6000000, 12000000, 24000000])
cryst_sara_e3 = np.array([0.302891, 0.389404, 0.432256])
cryst_sara_e4 = np.array([0.306425, 0.388235, 0.436959])

cryst_mike_e3 = np.array([0.281242, 0.372067, 0.4216267])
cryst_mike_e4 = np.array([0.277704, 0.372902, 0.427832])


plt.scatter(times, cryst_sara_e3, marker = "1", c = "red", label = r"$\dot{T}=10^{-3}$, nem20")
plt.scatter(times, cryst_sara_e4, marker = "1", c = "b", label = r"$\dot{T}=10^{-4}$, nem20")
plt.scatter(times, cryst_mike_e3, marker = "2", c = "r", label = r"$\dot{T}=10^{-3}$, python")
plt.scatter(times, cryst_mike_e4, marker = "2", c = "b", label = r"$\dot{T}=10^{-4}$, python")
plt.title(r"Comparison crystallisation, PVA-100, T = 0.85")
plt.xlabel("time")
plt.legend()
plt.ylabel(r"$\phi(t)$")
plt.savefig("plots/phit_selected_values_comparison.pdf")
plt.show()