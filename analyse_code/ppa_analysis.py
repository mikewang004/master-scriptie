import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 

def calc_normalised_entanglement_length(ppa, poly_length):
    #ppa_mean = np.mean(ppa)
    denom = poly_length/ppa
    print(np.std(ppa))
    #ne = np.mean(poly_length / (denom-1))
    ne = np.mean(poly_length/ppa)
    err = (poly_length**2/((poly_length - np.mean(ppa))**2) * np.std(ppa))
    print(err)
    return ne, err

ppa_50 = np.loadtxt("ppa_length_pva_50.txt")
ppa_100 = np.loadtxt("ppa_length_pva_100.txt")
ppa_200 = np.loadtxt("ppa_length_pva_200.txt")
ppa_300 = np.loadtxt("ppa_length_pva_300.txt")
ppa_500 = np.loadtxt("ppa_length_pva_500.txt")
ppa_1000 = np.loadtxt("ppa_length_pva_1000.txt")

ne_50, ne_err_50 = calc_normalised_entanglement_length(ppa_50, 50)
ne_100, ne_err_100 = calc_normalised_entanglement_length(ppa_100, 100)
ne_200, ne_err_200 = calc_normalised_entanglement_length(ppa_200, 200)
ne_300, ne_err_300 = calc_normalised_entanglement_length(ppa_300, 300)
ne_500, ne_err_500 = calc_normalised_entanglement_length(ppa_500, 500)
ne_1000, ne_err_1000 = calc_normalised_entanglement_length(ppa_1000, 1000)


polymer_lengths = np.array([50, 100, 200, 300, 500, 1000])

#plt.scatter(polymer_lengths, np.array([np.mean(ppa_50), np.mean(ppa_100),np.mean(ppa_200), np.mean(ppa_300),np.mean(ppa_500), np.mean(ppa_1000)]))
plt.scatter(polymer_lengths, np.array([ne_50, ne_100, ne_200, ne_300, ne_500, ne_1000])),# yerr = np.array([ne_err_50,ne_err_100,ne_err_200,ne_err_300,ne_err_500,ne_err_1000]), fmt = ".")
plt.title("Entanglement lengths")
plt.ylabel("N_e")
plt.xlabel("polymer length")
plt.savefig("PPA_t=18_different_polymers_unnormalised.pdf")
plt.show()