#!/usr/bin/env python3

plaquette_holonomy_u1.py

Figure 7 (illustrative): Gauge fields from lattice links and plaquettes (U(1) case).

We build a 2D square lattice, define link variables
    Ux(x,y) = exp(i a A_x(x,y))
    Uy(x,y) = exp(i a A_y(x,y))

For a constant magnetic field B (Landau gauge):
    A_x = 0
    A_y = B x

The plaquette holonomy is:
    U_p(x,y) = Ux(x,y) * Uy(x+1,y) * conj(Ux(x,y+1)) * conj(Uy(x,y))

The plaquette angle theta_p = arg(U_p) satisfies:
    theta_p ≈ a^2 F_xy = a^2 B   (for smooth fields, small a)

We demonstrate convergence by computing <theta_p>/a^2 for decreasing a.

Dependencies: numpy, matplotlib
Author: (your name)
"""

from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt


def build_u1_links(N: int, a: float, B: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Build U(1) link variables Ux, Uy on an NxN periodic lattice.
    Using Landau gauge for constant B:
        A_x = 0
        A_y = B x
    """
    x = (np.arange(N) - N // 2) * a  # centered coordinates
    y = (np.arange(N) - N // 2) * a

    X, Y = np.meshgrid(x, y, indexing="ij")

    Ax = np.zeros_like(X)
    Ay = B * X

    Ux = np.exp(1j * a * Ax)
    Uy = np.exp(1j * a * Ay)
    return Ux, Uy


def plaquette_holonomy(Ux: np.ndarray, Uy: np.ndarray) -> np.ndarray:
    """
    Compute plaquette holonomy U_p(x,y) with periodic boundary conditions:
        U_p = Ux(x,y) * Uy(x+1,y) * conj(Ux(x,y+1)) * conj(Uy(x,y))
    """
    Ux_xy = Ux
    Uy_x1y = np.roll(Uy, shift=-1, axis=0)           # x+1
    Ux_xy1 = np.roll(Ux, shift=-1, axis=1)           # y+1
    Uy_xy = Uy

    Up = Ux_xy * Uy_x1y * np.conjugate(Ux_xy1) * np.conjugate(Uy_xy)
    return Up


def wrapped_phase(z: np.ndarray) -> np.ndarray:
    """Return principal value phase in (-pi, pi]."""
    return np.angle(z)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="U(1) lattice links + plaquette holonomy; show theta/a^2 -> B as a -> 0."
    )
    parser.add_argument("--N", type=int, default=128, help="Lattice size N (default: 128)")
    parser.add_argument("--B", type=float, default=0.02, help="Constant field strength B (default: 0.02)")
    parser.add_argument(
        "--a_values",
        type=float,
        nargs="+",
        default=[2.0, 1.0, 0.5, 0.25, 0.125],
        help="Lattice spacings to test (default: 2 1 0.5 0.25 0.125)",
    )
    parser.add_argument(
        "--save",
        type=str,
        default="fig7_plaquette_convergence.png",
        help="Output filename (default: fig7_plaquette_convergence.png)",
    )
    args = parser.parse_args()

    N = int(args.N)
    B = float(args.B)

    a_vals = []
    est_vals = []
    ratios = []

    print("\nPlaquette holonomy convergence (U(1), constant B)")
    print(f"Parameters: N={N}, B={B}\n")
    header = f"{'a':>10} {'<theta>':>14} {'<theta>/a^2':>14} {'ratio to B':>14}"
    print(header)
    print("-" * len(header))

    for a in args.a_values:
        a = float(a)

        Ux, Uy = build_u1_links(N=N, a=a, B=B)
        Up = plaquette_holonomy(Ux, Uy)
        theta = wrapped_phase(Up)

        # Use mean over lattice as estimator (constant field => uniform theta)
        theta_mean = float(np.mean(theta))
        est = theta_mean / (a**2)

        a_vals.append(a)
        est_vals.append(est)
        ratios.append(est / B if B != 0 else np.nan)

        print(f"{a:>10.6g} {theta_mean:>14.6g} {est:>14.6g} {ratios[-1]:>14.6g}")

    a_vals = np.array(a_vals, dtype=float)
    est_vals = np.array(est_vals, dtype=float)
    ratios = np.array(ratios, dtype=float)

    # Plot: convergence of estimator toward B
    plt.figure(figsize=(10, 6))
    plt.plot(a_vals, est_vals, marker="o", label=r"Estimator $\langle\theta_\square\rangle/a^2$")
    plt.axhline(B, linestyle="--", label=r"Target $F_{xy}=B$")
    plt.gca().invert_xaxis()
    plt.xlabel("Lattice spacing a (smaller → smoother continuum limit)")
    plt.ylabel(r"Field strength estimate")
    plt.title("Plaquette Holonomy: $\u27E8\\theta_{\\square}\\u27E9/a^2 \\to F_{xy}$ as $a\\to 0$")
    plt.legend()
    plt.tight_layout()

    plt.savefig(args.save, dpi=300)
    print(f"\nSaved Figure 7 plot to: {args.save}")


if __name__ == "__main__":
    main()
