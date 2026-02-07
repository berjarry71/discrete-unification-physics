import numpy as np
import matplotlib.pyplot as plt

# Parameters
lambda_vals = np.linspace(0.1, 2.0, 30)
v = 1.0
g = 1.0

# Effective Higgs condensate
phi_sq_avg = v**2 * np.ones_like(lambda_vals)

# Gauge boson mass squared
mA_sq = g**2 * phi_sq_avg

# Plot
plt.figure(figsize=(7,5))
plt.plot(lambda_vals, mA_sq, 'o-', label=r"$m_A^2 \propto g^2\langle|\phi|^2\rangle$")

plt.xlabel(r"Higgs self-coupling $\lambda$")
plt.ylabel(r"Gauge boson mass squared $m_A^2$")
plt.title("Gauge boson mass generation from lattice Higgs field")

plt.legend()
plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()
