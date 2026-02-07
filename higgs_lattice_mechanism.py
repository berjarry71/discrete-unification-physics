#!/usr/bin/env python3
"""
regge_emergent_geometry.py

Figure 9 (illustrative): Emergent spacetime geometry from lattice dynamics.

We construct a simple 2D triangulated lattice with variable edge lengths.
Curvature is encoded via deficit angles at vertices, following Regge calculus.

Goal:
- Show how local fluctuations of link lengths generate curvature.
- Illustrate the emergence of curved geometry from discrete structure.

This is a conceptual illustration, not a full numerical GR solver.

Dependencies: numpy, matplotlib
Author: (your name)
"""

from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt


def triangle_angle(a: float, b: float, c: float) -> float:
    """
    Return the angle opposite side a in a triangle with sides (a, b, c),
    using the law of cosines.
    """
    cosA = (b**2 + c**2 - a**2) / (2.0 * b * c)
    cosA = np.clip(cosA, -1.0, 1.0)
    return np.arccos(cosA)


def deficit_angle(edge_lengths: list[tuple[float, float, float]]) -> float:
    """
    Compute the deficit angle around a vertex given the list of triangles
    sharing that vertex.

    edge_lengths: list of (a, b, c) where 'a' is the edge opposite the vertex.
    """
    angle_sum = 0.0
    for a, b, c in edge_lengths:
        angle_sum += triangle_angle(a, b, c)
    return 2.0 * np.pi - angle_sum


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Illustrative Regge calculus: curvature from deficit angles."
    )
    parser.add_argument(
        "--noise",
        type=float,
        default=0.15,
        help="Relative fluctuation amplitude of edge lengths (default: 0.15)",
    )
    parser.add_argument(
        "--save",
        type=str,
        default="fig9_regge_curvature.png",
        help="Output filename (default: fig9_regge_curvature.png)",
    )
    args = parser.parse_args()

    np.random.seed(1)

    # Reference equilateral triangle side length
    a0 = 1.0

    # Build a simple vertex with 6 surrounding triangles (flat case)
    n_tri = 6

    # Introduce fluctuations in edge lengths
    triangles = []
    for _ in range(n_tri):
        a = a0 * (1.0 + args.noise * np.random.randn())
        b = a0 * (1.0 + args.noise * np.random.randn())
        c = a0 * (1.0 + args.noise * np.random.randn())
        triangles.append((a, b, c))

    delta = deficit_angle(triangles)

    print("\nRegge curvature illustration")
    print(f"Number of triangles around vertex: {n_tri}")
    print(f"Relative edge-length noise: {args.noise}")
    print(f"Deficit angle Δ = {delta:.6f} rad")
    print(f"Effective curvature sign: {'positive' if delta > 0 else 'negative'}")

    # Visualization: schematic curvature indicator
    theta = np.linspace(0, 2 * np.pi, 400)
    r_flat = np.ones_like(theta)
    r_curved = r_flat * (1.0 + 0.5 * delta / np.pi)

    plt.figure(figsize=(6, 6))
    plt.plot(r_flat * np.cos(theta), r_flat * np.sin(theta), linestyle="--", label="Flat geometry")
    plt.plot(r_curved * np.cos(theta), r_curved * np.sin(theta), label="Emergent curved geometry")
    plt.gca().set_aspect("equal")
    plt.axis("off")
    plt.title("Emergent Curvature from Lattice Deficit Angles")
    plt.legend()
    plt.tight_layout()

    plt.savefig(args.save, dpi=300)
    print(f"\nSaved Figure 9 plot to: {args.save}")


if __name__ == "__main__":
    main()
