import numpy as np
import matplotlib.pyplot as plt

# Physical constants (dimensionless)
hbar = 1.0
m = 1.0

# Lattice parameters
L = 200
sigma = 10.0
k0 = 1.0

a_values = np.linspace(0.2, 2.0, 20)
uncertainty_products = []

for a in a_values:
    x = np.arange(-L//2, L//2) * a

    # Discrete Gaussian wave packet
    psi = np.exp(-x**2 / (2 * sigma**2)) * np.exp(1j * k0 * x)
    psi /= np.sqrt(np.sum(np.abs(psi)**2))

    # Position uncertainty
    x_mean = np.sum(np.abs(psi)**2 * x)
    dx = np.sqrt(np.sum(np.abs(psi)**2 * (x - x_mean)**2))

    # Momentum space (discrete Fourier transform)
    p = np.fft.fftfreq(len(x), d=a) * 2 * np.pi
    phi = np.fft.fft(psi)
    phi /= np.sqrt(np.sum(np.abs(phi)**2))

    p_mean = np.sum(np.abs(phi)**2 * p)
    dp = np.sqrt(np.sum(np.abs(phi)**2 * (p - p_mean)**2))

    uncertainty_products.append(dx * dp)

# Plot
plt.figure(figsize=(7,5))
plt.plot(a_values, uncertainty_products, 'o-', label=r"$\Delta x \Delta p$")
plt.axhline(hbar/2, color='r', linestyle='--', label=r"$\hbar/2$")

plt.xlabel("Lattice spacing $a$")
plt.ylabel(r"$\Delta x \Delta p$")
plt.title("Uncertainty relation from discrete lattice dynamics")
plt.legend()
plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()

