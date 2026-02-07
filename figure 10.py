import numpy as np
import matplotlib.pyplot as plt

# Time grid
t = np.linspace(-10, 10, 1000)
dt = t[1] - t[0]

# Parameters
alpha = 0.02
beta = 0.001

# Initial conditions
a = np.zeros_like(t)
v = np.zeros_like(t)

a[0] = 1.0
v[0] = 0.0

# Effective potential force
def force(a):
    return 2 / a**3 - 2 * alpha * a

# Integrate equation of motion
for i in range(len(t) - 1):
    v[i+1] = v[i] + force(a[i]) * dt
    a[i+1] = a[i] + v[i+1] * dt

# Plot
plt.figure(figsize=(7,5))
plt.plot(t, a, linewidth=2)

plt.xlabel("Cosmic time $t$")
plt.ylabel("Scale factor $a(t)$")
plt.title("Emergent cosmological evolution from lattice dynamics")

plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()

