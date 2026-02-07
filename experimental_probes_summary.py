#!/usr/bin/env python3
"""
experimental_probes_scale_chart.py

Figure 12 (paper-ready style): Experimental probes of a discrete lattice framework
shown on a logarithmic length scale with an approximate corresponding energy scale.

- Bottom axis: log10(Length [m])
- Top axis: log10(Energy [eV]) using E ≈ ħ c / L (order-of-magnitude mapping)

Bands indicate qualitative sensitivity regimes for different experimental channels.
This is an illustrative summary chart, not a constraint plot.

Dependencies: numpy, matplotlib
Author: (your name)
"""

from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt


# Physical constants for mapping E ≈ ħ c / L (approximate, order-of-magnitude)
HBAR = 1.054_571_817e-34  # J*s
C = 299_792_458.0         # m/s
EV = 1.602_176_634e-19    # J/eV


def logE_from_logL(logL: np.ndarray) -> np.ndarray:
    """
    Map log10(L[m]) -> log10(E[eV]) using E ≈ ħ c / L.
    """
    L = 10.0 ** logL
    E_eV = (HBAR * C) / (L * EV)
    return np.log10(E_eV)


def logL_from_logE(logE: np.ndarray) -> np.ndarray:
    """
    Inverse map log10(E[eV]) -> log10(L[m]) using L ≈ ħ c / E.
    """
    E = 10.0 ** logE
    L_m = (HBAR * C) / (E * EV)
    return np.log10(L_m)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a log-scale chart of experimental probes vs length/energy."
    )
    parser.add_argument(
        "--save",
        type=str,
        default="fig12_experimental_probes_scale_chart.png",
        help="Output filename (default: fig12_experimental_probes_scale_chart.png)",
    )
    args = parser.parse_args()

    # Length range: from subnuclear to cosmological scales
    logL_min, logL_max = -20.0, 26.0  # meters

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_axes([0.10, 0.16, 0.84, 0.72])

    ax.set_xlim(logL_min, logL_max)
    ax.set_ylim(0.0, 6.0)
    ax.set_yticks([])

    ax.set_xlabel(r"$\log_{10}(\mathrm{Length\ [m]})$")
    ax.set_title("Experimental Probes of the Discrete Lattice Framework (Illustrative)")

    # Top axis for energy (approximate mapping)
    ax_top = ax.twiny()
    ax_top.set_xlim(logE_from_logL(np.array([logL_min, logL_max])))
    ax_top.set_xlabel(r"$\log_{10}(\mathrm{Energy\ [eV]})$  (approx.\ via $E \simeq \hbar c/L$)")

    # Define bands: (name, y_center, y_height, logL_left, logL_right, notes)
    bands = [
        ("Particle physics (colliders)", 5.1, 0.7, -19.0, -15.0,
         "High-energy scattering\npossible lattice remnants,\nfermion sector signatures"),
        ("Precision tests (atomic/quantum)", 4.0, 0.7, -12.0, -6.0,
         "Interferometry,\nclock comparisons,\nphase-sensitive tests"),
        ("Astrophysics (photons/neutrinos)", 2.9, 0.7, -18.0, 6.0,
         "Time-of-flight dispersion,\nthreshold anomalies,\npropagation effects"),
        ("Gravitational waves", 1.8, 0.7, -3.0, 12.0,
         "Propagation / polarization,\nhigh-frequency deviations"),
        ("Cosmology (CMB/LSS/BBN)", 0.7, 0.7, 12.0, 26.0,
         "Bounce imprints,\nprimordial correlations,\nlate-time acceleration"),
    ]

    # Draw bands as rectangles with text
    for name, y, h, x0, x1, note in bands:
        rect = plt.Rectangle((x0, y - h/2), x1 - x0, h, fill=False, linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x0 + 0.2, y + 0.18, name, ha="left", va="center", fontsize=10, weight="bold")
        ax.text(x0 + 0.2, y - 0.18, note, ha="left", va="center", fontsize=9)

    # Add vertical guide lines for reference scales (optional but useful)
    ref_scales = [
        ("Proton\n~1 fm", -15.0),
        ("Atom\n~1 Å", -10.0),
        ("Optical\n~500 nm", -6.3),
        ("Human\n~1 m", 0.0),
        ("Earth\n~10^7 m", 7.0),
        ("Hubble\n~10^26 m", 26.0),
    ]
    for label, x in ref_scales:
        ax.vlines(x, 0.0, 6.0, linestyles="dotted", linewidth=1.0)
        ax.text(x, 5.85, label, rotation=90, ha="right", va="top", fontsize=8)

    # Add a few "signature" callouts as markers
    signatures = [
        ("Modified\ndispersion", -16.5, 3.3),
        ("Lorentz tests", -17.5, 5.3),
        ("Plaquette / gauge\ncontinuum limit", -16.0, 5.0),
        ("Quantum bounce\nimprints", 18.0, 0.9),
        ("Residual tension\n(late accel.)", 22.0, 0.5),
    ]
    for txt, x, y in signatures:
        ax.plot([x], [y], marker="o")
        ax.text(x + 0.4, y, txt, ha="left", va="center", fontsize=8)

    # Styling tweaks for print friendliness
    for spine in ["left", "right", "top"]:
        ax.spines[spine].set_visible(False)

    ax.grid(axis="x", linestyle=":", linewidth=0.8)

    # Footer note
    fig.text(
        0.10,
        0.06,
        "Note: bands indicate qualitative sensitivity regimes. The top axis uses an order-of-magnitude mapping E ≈ ħc/L.",
        fontsize=9,
    )

    plt.savefig(args.save, dpi=300)
    print(f"Saved Figure 12 chart to: {args.save}")


if __name__ == "__main__":
    main()
