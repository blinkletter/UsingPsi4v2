"""Helper functions for the UsingPsi4 notebooks.

Usage:
    from psi4_helpers import show_molecule, bond_length, bond_angle, plot_history, search_file
"""

import numpy as np
import matplotlib.pyplot as plt
import psi4
import py3Dmol


def show_molecule(mol):
    """Display a Psi4 molecule as an interactive 3D ball-and-stick model."""
    view = py3Dmol.view(width=400, height=300)
    view.addModel(mol.save_string_xyz_file(), "xyz")
    view.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
    view.setBackgroundColor("#f2f2f2")
    view.zoomTo()
    view.show()


def bond_length(mol, i, j):
    """Return the distance between atoms i and j in angstroms (0-based indices)."""
    xyz = mol.geometry().np * psi4.constants.bohr2angstroms
    return np.linalg.norm(xyz[i] - xyz[j])


def bond_angle(mol, i, j, k):
    """Return the i-j-k angle in degrees, with j as the central atom (0-based indices)."""
    xyz = mol.geometry().np
    a = xyz[i] - xyz[j]
    b = xyz[k] - xyz[j]
    cos_angle = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return np.degrees(np.arccos(cos_angle))


def plot_history(history, title):
    """Plot energy relative to the final step (kJ/mol) for a psi4.optimize history."""
    energies = np.array(history["energy"])
    delta_E = (energies - energies[-1]) * psi4.constants.hartree2kJmol
    steps = np.arange(1, len(energies) + 1)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))
    ax.plot(steps, delta_E, color="black", linewidth=0.7, marker="o",
            markerfacecolor="white", markeredgecolor="black", markeredgewidth=0.5)
    ax.set_xticks(steps)
    ax.set_xlabel("Optimization step")
    ax.set_ylabel(r"$E - E_{final}$ (kJ/mol)")
    ax.set_title(title)
    plt.show()


def search_file(filename, text):
    """Print every line in filename containing text (case-insensitive), with line numbers."""
    with open(filename) as f:
        for number, line in enumerate(f, start=1):
            if text.lower() in line.lower():
                print(f"{number:6d}: {line.rstrip()}")
