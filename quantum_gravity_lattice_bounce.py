#!/usr/bin/env python3
"""
cosmology_emergent_evolution.py

Figure 11 (illustrative): Emergent cosmological evolution from the lattice.

We build a simple effective evolution for the scale factor a(t) that captures:
1) A quantum bounce (no singularity) via a repulsive term at small a
2) Early accelerated expansion driven by quantum geometry (strong near bounce)
3) Late-time acceleration driven by a residual lattice tension (small constant push)

This is a conceptual toy model meant for visualization and intuition.
It is not a precision cosmology fit.

Dependencies: numpy, matplotlib
Author: (your name)
"""

from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt


def H_effective(a: float, rho0: float, a0: float, lam_eff: float) -> float:
    """
    A schematic Hubble-like term from matter + effective late-time tension.
    rho(a) ~ rho0 (a0/a)^3  (dust-like)
    plus a small constant term lam_eff (lattice residual tension analogue).
    """
    rho = rho0 * (a0 / max(a, 1e-12)) ** 3
    return np.sqrt(max(rho + lam_eff, 0.0))


def quantum_repulsion(a: float, alpha: float) -> float:
    """
    Repulsive quantum-geometry correction dominant at small a,
    producing a bounce. Chosen as ~ alpha / a^3.
    """
    return alpha / max(a**3, 1e-12)


def early_acceleration_boost(a: float, beta: float, a_star: float) -> float:
    """
    Early-time acceleration term (dominant near bounce) that decays as a grows.
    We model it as:
        beta * exp(-(a/a_star)^2)
    """
    return beta * np.exp(-((a / max(a_star, 1e-12)) ** 2))


def evolve(a_init: float, dt: float, steps: int, params: dict) -> tuple[np.ndarray, np.ndarray]:
    """
    Discrete-time evolution for a(t):

    We integrate a first-order effective equation:
        da/dt = a * H(a) + (early boost) + (quantum repulsion)

    - a*H(a): standard expansion term
    - early boost: drives early accelerated expansion near bounce
    - quantum repulsion: prevents a -> 0 and creates a bounce

    We simulate both backward and forward in time around the bounce by running
    forward from a_init and also generating a mirrored branch for illustration.
    """
    rho0 = params["rho0"]
    a0 = params["a0"]
    lam_eff = params["lam_eff"]
    alpha = params["alpha"]
    beta = params["beta"]
    a_star = params["a_star"]

    a = np.zeros(steps)
    t = np.arange(steps) * dt
    a[0] = a_init

    for n in range(1, steps):
        H = H_effective(a[n - 1], rho0=rho0, a0=a0, lam_eff=lam_eff)
        boost = early_acceleration_boost(a[n - 1], beta=beta, a_star=a_star)
        rep = quantum_repulsion(a[n - 1], alpha=alpha)

        da = dt * (a[n - 1] * H + boost + rep)
        a[n] = max(a[n - 1] + da, 1e-10)

    # Build a "pre-bounce" branch by reversing time (visual symmetry)
    # This is illustrative: we mirror around t=0 with same |t|.
    t_full = np.concatenate((-t[::-1], t[1:]))
    a_full = np.concatenate((a[::-1], a[1:]))

    return t_full, a_full


def finite_diff_second(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Second derivative via finite differences (for acceleration sign)."""
    dy = np.gradient(y, t)
    d2y = np.gradient(dy, t)
    return d2y


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Illustrative emergent cosmological evolution: bounce + early and late acceleration."
    )
    parser.add_argument("--a_init", type=float, default=0.08, help="Initial scale factor near bounce (default: 0.08)")
    parser.add_argument("--dt", type=float, default=0.01, help="Time step (default: 0.01)")
    parser.add_argument("--steps", type=int, default=2500, help="Number of forward steps (default: 2500)")

    # Matter + late-time term
    parser.add_argument("--rho0", type=float, default=0.30, help="Matter density scale (default: 0.30)")
    parser.add_argument("--a0", type=float, default=1.0, help="Reference scale factor (default: 1.0)")
    parser.add_argument("--lam_eff", type=float, default=0.003, help="Residual lattice tension (default: 0.003)")

    # Quantum + early terms
    parser.add_argument("--alpha", type=float, default=0.0008, help="Quantum repulsion strength (default: 0.0008)")
    parser.add_argument("--beta", type=float, default=0.05, help="Early acceleration boost (default: 0.05)")
    parser.add_argument("--a_star", type=float, default=0.25, help="Decay scale for early boost (default: 0.25)")

    parser.add_argument("--save", type=str, default="fig11_emergent_cosmology.png", help="Output filename")
    args = parser.parse_args()

    params = {
        "rho0": args.rho0,
        "a0": args.a0,
        "lam_eff": args.lam_eff,
        "alpha": args.alpha,
        "beta": args.beta,
        "a_star": args.a_star,
    }

    t, a = evolve(a_init=args.a_init, dt=args.dt, steps=args.steps, params=params)

    # Compute acceleration indicator: second derivative of a(t)
    acc = finite_diff_second(t, a)

    # Determine regions of acceleration (acc > 0)
    # We'll mark early and late acceleration zones by simple thresholds.
    early_zone = (np.abs(t) < 2.0)  # near bounce
    late_zone = (np.abs(t) > 12.0)  # far future/past (illustrative)

    plt.figure(figsize=(10, 6))
    plt.plot(t, a, label="Scale factor a(t) with bounce")

    # Mark early acceleration region
    plt.fill_between(
        t,
        a.min(),
        a.max(),
        where=(early_zone & (acc > 0)),
        alpha=0.2,
        label="Early accelerated expansion (quantum geometry)",
    )

    # Mark late-time acceleration region
    plt.fill_between(
        t,
        a.min(),
        a.max(),
        where=(late_zone & (acc > 0)),
        alpha=0.2,
        label="Late-time acceleration (residual lattice tension)",
    )

    plt.xlabel("Cosmic time (schematic)")
    plt.ylabel("Scale factor a(t)")
    plt.title("Emergent Cosmological Evolution from the Lattice (Illustrative)")
    plt.legend()
    plt.tight_layout()

    plt.savefig(args.save, dpi=300)
    print(f"Saved Figure 11 plot to: {args.save}")


if __name__ == "__main__":
    main()
