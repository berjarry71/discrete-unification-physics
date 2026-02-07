import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 1.0
epsilon = 0.15
N = 20

# Regular 2D lattice
x, y = np.meshgrid(np.arange(N), np.arange(N))
points = np.column_stack((x.flatten(), y.flatten()))

# Introduce random link-length fluctuations
np.random.seed(0)
height = epsilon * np.random.randn(N, N)

# Compute discrete curvature proxy (Laplacian of height field)
curvature = (
    -4 * height
    + np.roll(height, 1, axis=0)
    + np.roll(height, -1, axis=0)
    + np.roll(height, 1, axis=1)
    + np.roll(height, -1, axis=1)
)

# Plot
plt.figure(figsize=(7,5))
plt.imshow(curvature, cmap='RdBu')
plt.colorbar(label="Discrete curvature (deficit proxy)")
plt.title("Emergent curvature from lattice link-length fluctuations")
plt.axis('off')
plt.tight_layout()
plt.show()
