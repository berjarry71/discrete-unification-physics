#!/usr/bin/env python3
"""
dispersion_naive_vs_z3.py

Figure 5 (illustrative): Naive vs Z3-projected fermion dispersion on a lattice.

Goal:
- Show the naive lattice Dirac dispersion with multiple low-energy "doubler" points.
- Show a Z3-projected effective dispersion where doubler branches are strongly suppressed.

This is a minimal, transparent toy model designed for reproducibility and pedagogy.
It is not a full lattice-QFT simulation.

Dependencies: numpy, matplotlib

Author: (your name)
"""

from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt


def naive_dirac_energy_1d(k: np.ndarray, a: float, m: float = 0.0, c: float = 1.0) -> np.ndarray:
    """
    Minimal 1D naive lattice Dirac-like dispersion (toy):
        p_lat(k) = (1/a) sin(k a)
        E(k) = sqrt( (c p_lat)^2 + (m c^2)^2 )

This captures the key feature: sin(k a) produces extra near-zeros at k ~ 0, ±π/a, ...
In higher dimensions, doubling becomes more pronounced, but 1D is enough to illustrate the mechanism.
"""
    p_lat = (1.0 / a) * np.sin(k * a)
    return np.sqrt((c * p_lat) ** 2 + (m * c**2) ** 2)


def z3_projector_weight(k: np.ndarray, a: float, sigma: float = 0.25) -> np.ndarray:
    """
    A simple smooth "Z3 projection" weight in momentum space (toy).

Idea:
- Keep the physical branch around k ~ 0.
- Strongly suppress contributions from the doubler region near |k| ~ π/a.

We implement a smooth damping centered at the doubler points:
    w(k) ~ 1 - exp(-((|k| - π/a)^2)/(2 sigma^2 (π/a)^2))
so that w ~ 1 near k=0 and w ~ 0 near k=±π/a.

sigma controls the width of the suppression region (dimensionless relative to π/a).
"""
    kmax = np.pi / a
    # distance from doubler shell
    d = np.abs(np.abs(k) - kmax) / kmax
    # near d=0 -> suppress; far from d=0 -> keep
    suppress = np.exp(-(d**2) / (2.0 * sigma**2))
    w = 1.0 - suppress
    return np.clip(w, 0.0, 1.0)


def projected_effective_energy(E_naive: np.ndarray, w: np.ndarray, floor: float = 1e-6) -> np.ndarray:
    """
    Build an illustrative "projected" effective dispersion.

We interpret the projection as suppressing the doubler branch contribution.
A simple way to visualize this is to keep energies where weight is high and push
energies where weight is low upward (or effectively suppress them).

Here we define:
    E_proj = E_naive / max(w, floor)

So near k=0: w~1 => E_proj ~ E_naive
Near k~π/a: w~0 => E_proj becomes large => doubler low-energy branch disappears visually.

This is purely illustrative but gives a clean, full-page plot for the paper.
"""
    w_eff = np.maximum(w, floor)
    return E_naive / w_eff


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Illustrative comparison: naive lattice fermion dispersion vs Z3-projected effective dispersion."
    )
    parser.add_argument("--a", type=float, default=1.0, help="Lattice spacing a (default: 1.0)")
    parser.add_argument("--m", type=float, default=0.05, help="Mass parameter m (default: 0.05)")
    parser.add_argument("--c", type=float, default=1.0, help="Speed parameter c (default: 1.0)")
    parser.add_argument(
        "--n",
        type=int,
        default=2000,
        help="Number of k-samples across the Brillouin zone (default: 2000)",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=0.22,
        help="Suppression width (dimensionless relative to π/a) (default: 0.22)",
    )
    parser.add_argument(
        "--save",
        type=str,
        default="fig5_naive_vs_z3_dispersion.png",
        help="Output filename (default: fig5_naive_vs_z3_dispersion.png)",
    )
    args = parser.parse_args()

    a = float(args.a)
    m = float(args.m)
    c = float(args.c)

    # k in the first Brillouin zone: [-π/a, π/a]
    kmax = np.pi / a
    k = np.linspace(-kmax, kmax, args.n)

    E_naive = naive_dirac_energy_1d(k, a=a, m=m, c=c)
    w = z3_projector_weight(k, a=a, sigma=args.sigma)
    E_proj = projected_effective_energy(E_naive, w=w)

    # Normalize for clean plotting (optional)
    # We normalize by the minimum of E_naive near k~0 to show suppression clearly
    E0 = float(np.min(E_naive))
    if E0 <= 0:
        E0 = 1.0
    E_naive_n = E_naive / E0
    E_proj_n = E_proj / E0

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(k * a, E_naive_n, label=r"Naive $E(k)$")
    plt.plot(k * a, E_proj_n, label=r"Z$_3$-projected $E_{\mathrm{proj}}(k)$")
    plt.xlabel(r"Dimensionless momentum $ka$  (Brillouin zone: $[-\pi, \pi]$)")
    plt.ylabel(r"Normalized energy $E/E_{\min}$")
    plt.title("Naive vs Z$_3$-Projected Fermion Dispersion (Illustrative)")
    plt.ylim(0, np.percentile(E_proj_n, 95))  # keep plot readable (clips very large suppressed region)
    plt.legend()
    plt.tight_layout()

    plt.savefig(args.save, dpi=300)
    print(f"Saved Figure 5 plot to: {args.save}")


if __name__ == "__main__":
    main()
