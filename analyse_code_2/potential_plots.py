import numpy as np
import matplotlib.pyplot as plt

# Parameters
k_bond = 2704
b0 = 0.5
eps = 1.511
sigma0 = 0.89

# --- Bond potential ---
r_bond = np.linspace(0.01, 2.0, 500)
U_bond = 0.5 * k_bond * (r_bond - b0)**2

# --- LJ potential with cutoff at 1.6 ---
r_lj = np.linspace(0.4, 1.6, 500)
U_lj = eps * ((sigma0 / r_lj)**9 - (sigma0 / r_lj)**6)

# --- Plotting ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Bond potential
axes[0].plot(r_bond, U_bond, color='steelblue', lw=2)
axes[0].axvline(b0, color='gray', linestyle='--', label=f'$b_0 = {b0}$')
axes[0].set_xlabel('$r$', fontsize=13)
axes[0].set_ylabel('$U_{\\mathrm{bond}}(r)$', fontsize=13)
axes[0].set_title('Harmonic Bond Potential', fontsize=14)
axes[0].legend(fontsize=11)
axes[0].set_xlim(0, 2)
axes[0].set_ylim(0, None)
axes[0].grid(True, alpha=0.3)

# LJ potential
axes[1].plot(r_lj, U_lj, color='tomato', lw=2)
axes[1].axhline(0, color='black', linestyle='-', lw=0.8)
axes[1].axvline(1.6, color='gray', linestyle='--', label='cutoff $r = 1.6$')
axes[1].set_xlabel('$r$', fontsize=13)
axes[1].set_ylabel('$U_{\\mathrm{LJ}}(r)$', fontsize=13)
axes[1].set_title('Lennard-Jones Potential (9-6)', fontsize=14)
axes[1].legend(fontsize=11)
axes[1].set_xlim(0.4, 2)
axes[1].set_ylim(-2, 5)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()