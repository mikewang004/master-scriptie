import numpy as np 

lpp = np.loadtxt("Lpp_values.dat")

print(lpp)
N = 500;
N_minus_1 = N - 1
mean_lpp = np.mean(lpp)
re_sqrt = 26.631469

lbpp = mean_lpp/N_minus_1
print(lbpp)

lk = re_sqrt**2/(N_minus_1 * lbpp)

ne = lk/lbpp

print(ne)