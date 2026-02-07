import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 1.0
alpha = 4.0

# Momentum range (Brillouin zone)
k = np.linspace(-np.pi/a, np.pi/a, 600)

# Naive fermion dispersion
E_naive = np.abs(np.sin(k * a)) / a

# Z3-projected dispersion (effective suppression of doublers)
suppression = np.exp(-alpha * np.sin(3 * k * a / 2)**2)
E_proj = E_naive * suppression

# Plot
plt.figure(figsize=(7,5))

plt.plot(k, E_naive, '--', label="Naive fermion dispersion")
plt.plot(k, E_proj, label=r"$\mathbb{Z}_3$-projected dispersion", linewidth=2)

plt.xlabel(r"Momentum $k$")
plt.ylabel(r"Energy $E(k)$")
plt.title("Fermion doubling and $\\mathbb{Z}_3$ projection")

plt.legend()
plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()
