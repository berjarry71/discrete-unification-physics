import numpy as np
import matplotlib.pyplot as plt

# Discrete scale factor
a_vals = np.linspace(0.01, 2.0, 400)

# Effective discrete action (toy quantum gravity)
lambda_c = 1.0
S = (1.0 / a_vals) + lambda_c * a_vals**2

# Quantum probability amplitude
P = np.exp(-S)

# Normalize
P /= np.trapz(P, a_vals)

# Plot
plt.figure(figsize=(7,5))
plt.plot(a_vals, P, linewidth=2)
plt.axvline(0, color='k', linestyle='--', alpha=0.5)

plt.xlabel("Scale factor $a$")
plt.ylabel("Quantum probability density $P(a)$")
plt.title("Discrete gravitational path integral and singularity resolution")

plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()
