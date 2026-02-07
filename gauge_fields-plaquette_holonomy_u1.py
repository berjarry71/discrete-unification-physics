#!/usr/bin/env python3
"""
higgs_lattice_mechanism.py

Figure 8 (illustrative): Higgs mechanism on a lattice.

We consider a complex scalar field phi on a 2D lattice with a
Mexican-hat potential:
    V(phi) = λ (|phi|^2 - v^2)^2

We minimize the local potential to illustrate spontaneous symmetry breaking
and show how fluctuations around the vacuum acquire an effective mass.

This is a conceptual lattice illustration, not a full gauge-Higgs simulation.

Dependencies: numpy, matplotlib
Author: (your name)
"""

from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt


def mexican_hat_potential(phi_abs: np.ndarray, lam: float, v: float) -> np.ndarray:
    """Mexican-hat potential V(|phi|)."""
    return lam * (phi_abs**2 - v**2) ** 2


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Illustrative Higgs mechanism on a lattice (Mexican-hat potential)."
    )
    parser.add_argument("--N", type=int, default=200, help="Grid resolution for |phi| (default: 200)")
    parser.add_argument("--lam", type=float, default=1.0, help="Self-coupling λ (default: 1.0)")
    parser.add_argument("--v", type=float, default=1.0, help="Vacuum expectation value v (default: 1.0)")
    parser.add_argument(
        "--save",
        type=str,
        default="fig8_higgs_lattice.png",
        help="Output filename (default: fig8_higgs_lattice.png)",
    )
    args = parser.parse_args()

    lam = float(args.lam)
    v = float(args.v)

    # Radial field values
    phi_abs = np.linspace(0.0, 2.0 * v, args.N)
    V = mexican_hat_potential(phi_abs, lam=lam, v=v)

    # Identify vacuum (minimum of potential)
    idx_min = np.argmin(V)
    phi_vac = phi_abs[idx_min]

    # Effective Higgs mass squared ~ second derivative of V at vacuum
    dphi = phi_abs[1] - phi_abs[0]
    d2V = np.gradient(np.gradient(V, dphi), dphi)
    m_h_sq = d2V[idx_min]

    print("\nHiggs lattice illustration")
    print(f"λ = {lam}, v = {v}")
    print(f"Vacuum expectation value |phi| ≈ {phi_vac:.4f}")
    print(f"Effective Higgs mass^2 ≈ {m_h_sq:.4f}")

    # Plot Mexican-hat radial potential
    plt.figure(figsize=(10, 6))
    plt.plot(phi_abs, V, label=r"$V(|\phi|) = \lambda(|\phi|^2 - v^2)^2$")
    plt.axvline(phi_vac, linestyle="--", label=r"Vacuum $|\phi| = v$")
    plt.scatter([phi_vac], [V[idx_min]], zorder=5)
    plt.xlabel(r"Field amplitude $|\phi|$")
    plt.ylabel(r"Potential $V$")
    plt.title("Higgs Mechanism on a Lattice: Spontaneous Symmetry Breaking")
    plt.legend()
    plt.tight_layout()

    plt.savefig(args.save, dpi=300)
    print(f"\nSaved Figure 8 plot to: {args.save}")


if __name__ == "__main__":
    main()
