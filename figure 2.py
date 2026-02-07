import numpy as np
import matplotlib.pyplot as plt

# Physical parameters (dimensionless units)
hbar = 1.0
m = 1.0
a = 1.0

# Momentum range (first Brillouin zone)
k = np.linspace(0, np.pi / a, 400)

# Lattice dispersion (nearest-neighbor Laplacian)
omega_lattice = (2 / (m * a**2)) * (1 - np.cos(k * a))

# Continuum Schrödinger dispersion
omega_cont = (hbar * k**2) / (2 * m)

# Plot
plt.figure(figsize=(7, 5))

plt.plot(k, omega_lattice, label="Lattice dispersion", linewidth=2)
plt.plot(k, omega_cont, '--', label="Continuum Schrödinger limit", linewidth=2)

plt.xlabel(r"Wave number $k$")
plt.ylabel(r"Frequency $\omega(k)$")
plt.title("Emergence of the Schrödinger dispersion relation")

plt.legend()
plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()
