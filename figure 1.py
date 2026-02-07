import matplotlib.pyplot as plt
import numpy as np
from itertools import product

# Lattice size
N = 4
a = 1.0

# Generate lattice points (3D projection of 4D lattice)
points = np.array(list(product(range(N), repeat=3)))

# Nearest-neighbor connections
edges = []
for p in points:
    for d in range(3):
        q = p.copy()
        q[d] += 1
        if q[d] < N:
            edges.append((p, q))

# Choose central site
center = np.array([1, 1, 1])
neighbors = []

for p in points:
    if np.sum(np.abs(p - center)) == 1:
        neighbors.append(p)

# Plot
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(projection='3d')

# Plot edges
for p, q in edges:
    ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]],
            color='gray', alpha=0.6)

# Plot lattice points
ax.scatter(points[:,0], points[:,1], points[:,2],
           s=30, color='black')

# Highlight center and neighbors
ax.scatter(*center, s=120, color='red', label='Site $x$')
for n in neighbors:
    ax.scatter(*n, s=80, color='blue')

ax.set_title("Discrete spacetime lattice with nearest-neighbor coupling")
ax.set_axis_off()
ax.legend()

plt.tight_layout()
plt.show()
