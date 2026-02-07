import numpy as np
import matplotlib.pyplot as plt

# Parameters
N = 90              # lattice sites in compact direction
M = N // 3          # Z3 shift
n = np.arange(N)

# Z3 eigenmodes
psi0 = np.ones(N) / np.sqrt(N)
psi1 = np.exp(2j * np.pi * n / 3) / np.sqrt(N)
psi2 = np.exp(4j * np.pi * n / 3) / np.sqrt(N)

# Plot
plt.figure(figsize=(8,5))

plt.plot(n, np.real(psi0), label=r"Sector $\lambda = 1$")
plt.plot(n, np.real(psi1), label=r"Sector $\lambda = e^{2\pi i/3}$")
plt.plot(n, np.real(psi2), label=r"Sector $\lambda = e^{4\pi i/3}$")

plt.xlabel("Compact lattice coordinate")
plt.ylabel("Mode amplitude (real part)")
plt.title(r"$\mathbb{Z}_3$ orbifold projection and invariant sectors")

plt.legend()
plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()
