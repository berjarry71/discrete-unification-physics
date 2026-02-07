import numpy as np
import matplotlib.pyplot as plt

# Energy range (GeV)
E = np.logspace(0, 19, 500)

# Discreteness scale (Planck scale)
Lambda = 1e19  # GeV
eta = 1.0

# Velocity deviation
delta_v = eta * (E / Lambda)**2

# Plot
plt.figure(figsize=(7,5))
plt.loglog(E, delta_v, linewidth=2)

plt.xlabel("Energy $E$ [GeV]")
plt.ylabel(r"Velocity deviation $\Delta v$")
plt.title("Universal dispersion corrections from spacetime discreteness")

plt.grid(which='both', alpha=0.4)
plt.tight_layout()
plt.show()


