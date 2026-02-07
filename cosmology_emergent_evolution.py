#!/usr/bin/env python3
"""
experimental_probes_summary.py

Figure 12 (illustrative): Experimental probes of a discrete lattice framework.

Generates a full-page summary diagram showing multiple experimental/observational
channels and the approximate scales (energy / length) they probe.

This is an illustrative infographic, not a data-driven constraint plot.

Dependencies: numpy, matplotlib
Author: (your name)
"""

from __future__ import annotations

import argparse
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate an infographic: experimental probes of a discrete lattice."
    )
    parser.add_argument(
        "--save",
        type=str,
        default="fig12_experimental_probes.png",
        help="Output filename (default: fig12_experimental_probes.png)",
    )
    args = parser.parse_args()

    # We build a diagram with a vertical axis representing "scale" (qualitative),
    # and horizontal categories for channels.
    # No hard numbers are claimed; we use qualitative bands with example labels.

    channels = [
        "Particle Physics\n(Colliders)",
        "Precision Tests\n(Atomic/Quantum)",
        "Astrophysics\n(Photons/Neutrinos)",
        "Gravitational Waves",
        "Cosmology\n(CMB/LSS/BBN)",
    ]

    # Qualitative vertical "scale" levels (top = smallest length / highest energy)
    y_levels = [
        ("Planck / UV frontier\n(very speculative)", 0.92),
        ("Ultra-High Energy\nAstroparticles", 0.78),
        ("TeV–PeV Scale\n(High-energy tests)", 0.64),
        ("GeV–TeV Scale\n(Precision SM)", 0.50),
        ("Low-energy / IR\n(Astrophysical & cosmological)", 0.36),
        ("Largest scales\n(Structure formation)", 0.20),
    ]

    # Each entry: (channel_index, y_center, height, label)
    # These are illustrative "bands" indicating where each probe is most sensitive.
    bands = [
        (0, 0.60, 0.22, "Lorentz tests,\nmodified dispersion,\nrare processes"),
        (0, 0.48, 0.16, "High-energy scattering,\nfermion sector structure"),
        (1, 0.46, 0.18, "Precision clocks,\ninterferometry,\nquantum phases"),
        (2, 0.70, 0.20, "Time-of-flight delays,\nenergy-dependent dispersion,\nthreshold anomalies"),
        (2, 0.56, 0.16, "Neutrino propagation,\nflavor effects"),
        (3, 0.42, 0.22, "Propagation effects,\npolarization,\nhigh-frequency deviations"),
        (4, 0.30, 0.24, "Primordial signatures,\nbounce footprints,\nCMB correlations"),
        (4, 0.22, 0.18, "Late-time acceleration,\nresidual lattice tension"),
    ]

    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_axes([0.08, 0.08, 0.84, 0.84])

    ax.set_xlim(-0.6, len(channels) - 0.4)
    ax.set_ylim(0.0, 1.0)

    # Remove axes for infographic style
    ax.set_xticks(range(len(channels)))
    ax.set_xticklabels(channels, fontsize=10)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title("Experimental Probes of the Discrete Lattice Framework (Illustrative)", fontsize=14, pad=20)

    # Draw horizontal scale separators + labels
    for label, y in y_levels:
        ax.hlines(y, -0.6, len(channels) - 0.4, linestyles="dotted", linewidth=1)
        ax.text(-0.62, y, label, va="center", ha="right", fontsize=9)

    # Draw channel vertical guides
    for i in range(len(channels)):
        ax.vlines(i, 0.08, 0.92, linestyles="dashed", linewidth=0.8)

    # Draw bands as rectangles with text
    for (ci, y, h, txt) in bands:
        left = ci - 0.35
        width = 0.70
        bottom = y - h / 2.0

        rect = plt.Rectangle((left, bottom), width, h, fill=False, linewidth=1.5)
        ax.add_patch(rect)
        ax.text(ci, y, txt, ha="center", va="center", fontsize=9)

    # Add a small legend note
    ax.text(
        0.5,
        0.05,
        "Bands indicate qualitative sensitivity regions; this figure summarizes targets rather than presenting constraints.",
        ha="center",
        va="center",
        fontsize=9,
        transform=ax.transAxes,
    )

    plt.savefig(args.save, dpi=300)
    print(f"Saved Figure 12 plot to: {args.save}")


if __name__ == "__main__":
    main()
