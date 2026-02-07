#!/usr/bin/env python3
"""
quantum_gravity_lattice_bounce.py

Figure 10 (illustrative): Quantum gravity from lattice configurations.

We model a simplified discrete path integral over geometric configurations
characterized by a scale factor a_n evolving in discrete "time" steps.

Key idea:
- Classical evolution leads to a -> 0 singularity
- Quantum (lattice) corrections generate a repulsive term
- The singularity is replaced by a quantum bounce

This is a conceptual illustration, not a full quantum gravity simulation.

Dependencies: numpy, matplotlib
Author: (your name)
"""

from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt


def classical_friedmann(a: float, rho: float) -> float:
    """
    Classical Friedmann-like contraction term (schematic):
        (da/dt)^2 ~ rho a^2
    """
    return -rho * a


def quantum_correction(a: float, alpha: float) -> float:
    """
    Lattice-induced quantum correction (repulsive at small a):
        ~ + alpha / a^3
    """
    return alpha / max(a**3, 1e-8)


def evolve_scale_factor(
    a0: float,
    rho: float,
    alpha: float,
    dt: float,
    n_steps: int,
) -> np.ndarray:
    """
    Discrete evolution of the scale factor including quantum correction.
    """
    a = np.zeros(n_steps)
    a[0] = a0

    for n in range(1, n_steps):
        force = classical_friedmann(a[n - 1], rho) + quantum_correction(a[n - 1], alpha)
        a[n] = max(a[n - 1] + dt * force, 1e-6)

    return a


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Illustrative quantum gravity lattice bounce."
    )
    parser.add_argument("--a0", type=float, default=1.0, help="Initial scale factor")
    parser.add_argument("--rho", type=float, default=0.4, help="Matter density parameter")
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.02,
        help="Quantum geometry strength (lattice correction)",
    )
    parser.add_argument("--dt", type=float, default=0.01, help="Discrete time step")
    parser.add_argument("--steps", type=int, default=2000, help="Number of steps")
    parser.add_argument(
        "--save",
        type=str,
        default="fig10_quantum_gravity_bounce.png",
        help="Output filename",
    )
    args = parser.parse_args()

    a_quantum = evolve_scale_factor(
        a0=args.a0,
        rho=args.rho,
        alpha=args.alpha,
        dt=args.dt,
        n_steps=args.steps,
    )

    # Classical reference (no quantum correction)
    a_classical = evolve_scale_factor(
        a0=args.a0,
        rho=args.rho,
        alpha=0.0,
        dt=args.dt,
        n_steps=args.steps,
    )

    t = np.arange(args.steps) * args.dt

    plt.figure(figsize=(10, 6))
    plt.plot(t, a_classical, linestyle="--", label="Classical evolution (singularity)")
    plt.plot(t, a_quantum, label="Quantum lattice evolution (bounce)")
    plt.xlabel("Discrete time")
    plt.ylabel("Scale factor a(t)")
    plt.title("Quantum Gravity from Lattice Configurations: Singularity vs Bounce")
    plt.legend()
    plt.tight_layout()

    plt.savefig(args.save, dpi=300)
    print(f"Saved Figure 10 plot to: {args.save}")


if __name__ == "__main__":
    main()
