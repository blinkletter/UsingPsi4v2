"""Helper functions for the UsingPsi4 notebooks.

Usage:
    import psi4_helpers as ph
    ph.show_molecule(mol, interactive=INTERACTIVE)

Every drawing helper takes ``interactive``:
    True  -> a rotatable py3Dmol viewer (needs a live notebook; large when saved)
    False -> a static matplotlib picture (small; works when the book is built)
Static isosurfaces use scikit-image (``skimage.measure.marching_cubes``).
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import psi4
import py3Dmol

BOHR2ANG = psi4.constants.bohr2angstroms

# Colours and covalent radii (Å) used for static ball-and-stick pictures
ELEMENTS = {
    "H": ("white", 0.31), "C": ("#555555", 0.76), "N": ("#3050F8", 0.71),
    "O": ("#FF0D0D", 0.66), "F": ("#90E050", 0.57), "P": ("#FF8000", 1.07),
    "S": ("#FFFF30", 1.05), "Cl": ("#1FF01F", 1.02), "Br": ("#A62929", 1.20),
}
Z2SYMBOL = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl", 35: "Br"}


# ------------------------------------------------------------------ static drawing tools

def _new_3d_axes(ax, figsize=(4, 4)):
    """Return (ax, created) — a new 3D axes if none was given."""
    if ax is not None:
        return ax, False
    fig = plt.figure(figsize=figsize)
    return fig.add_subplot(projection="3d"), True


def _finish_axes(ax, points, pad=0.3, elev=20, azim=35, title="", zoom=1.4):
    """Equal-aspect limits that just enclose points (Å), set the view, hide the axes."""
    centre = (points.max(0) + points.min(0)) / 2
    half = np.abs(points - centre).max() + pad
    for set_lim, c in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), centre):
        set_lim(c - half, c + half)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=9)


def _draw_atoms(ax, symbols, xyz, scale=1.0):
    """Ball-and-stick model (xyz in Å) on a 3D axes."""
    cols = [ELEMENTS.get(s, ("pink", 0.8))[0] for s in symbols]
    radii = np.array([ELEMENTS.get(s, ("pink", 0.8))[1] for s in symbols])
    for i in range(len(symbols)):
        for j in range(i + 1, len(symbols)):
            if np.linalg.norm(xyz[i] - xyz[j]) < 1.2 * (radii[i] + radii[j]):
                ax.plot(*xyz[[i, j]].T, color="#444444", linewidth=2.5 * scale)
    ax.scatter(*xyz.T, s=900 * scale * radii**2, c=cols, edgecolors="black",
               linewidths=0.5, depthshade=False)


def _shaded(tris, colors, light=(1.0, -0.6, 1.2)):
    """Darken each triangle's colour by how far it faces away from the light."""
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    light = np.asarray(light) / np.linalg.norm(light)
    k = 0.45 + 0.55 * np.abs(n @ light)
    rgba = np.array(mcolors.to_rgba_array(colors), float)
    if len(rgba) == 1:
        rgba = np.repeat(rgba, len(tris), 0)
    rgba[:, :3] *= k[:, None]
    return rgba


def _opacity(transparency, default):
    """Opacity (0-1) from a transparency fraction (0 = solid, 1 = invisible); None keeps the default."""
    return default if transparency is None else 1.0 - transparency


def _draw_surface(ax, tris, colors, alpha=0.7):
    if len(tris):
        ax.add_collection3d(Poly3DCollection(tris, facecolors=_shaded(tris, colors),
                                             edgecolors="none", linewidths=0,      # no mesh lines
                                             antialiased=True, alpha=alpha))


def _molecule_atoms(mol):
    """Element symbols and coordinates (Å) of a psi4 Molecule."""
    return [mol.symbol(i).capitalize() for i in range(mol.natom())], mol.geometry().np * BOHR2ANG


def _xyz_frames(xyz_text):
    """Split a (multi-frame) XYZ string into a list of (comment, symbols, xyz Å)."""
    lines = xyz_text.strip().splitlines()
    frames, i = [], 0
    while i < len(lines):
        n = int(lines[i].split()[0])
        atoms = [l.split() for l in lines[i + 2:i + 2 + n]]
        frames.append((lines[i + 1].strip(), [a[0].capitalize() for a in atoms],
                       np.array([a[1:4] for a in atoms], float)))
        i += n + 2
    return frames


# ------------------------------------------------------------------ cube files

def read_cube(path):
    """Read a Gaussian cube file into a dict: data grid, origin and axes (Å), symbols, xyz (Å)."""
    lines = Path(path).read_text().splitlines()
    natoms, *origin = lines[2].split()
    natoms = int(natoms)
    shape, axes = [], []
    for line in lines[3:6]:
        n, *vec = line.split()
        shape.append(int(n))
        axes.append(np.array(vec[:3], float))
    atom_lines = [l.split() for l in lines[6:6 + abs(natoms)]]
    start = 6 + abs(natoms) + (1 if natoms < 0 else 0)            # skip the MO-index line if present
    data = np.array(" ".join(lines[start:]).split(), float).reshape(shape)
    return {"data": data,
            "origin": np.array(origin[:3], float) * BOHR2ANG,
            "axes": np.array(axes) * BOHR2ANG,
            "symbols": [Z2SYMBOL.get(int(a[0]), "X") for a in atom_lines],
            "xyz": np.array([a[2:5] for a in atom_lines], float) * BOHR2ANG}


def isosurface(cube, isoval):
    """Triangles (n, 3, 3) in Å where the cube data equals isoval (scikit-image marching cubes)."""
    from skimage.measure import marching_cubes
    data = cube["data"]
    if not data.min() < isoval < data.max():
        return np.empty((0, 3, 3))
    verts, faces, _, _ = marching_cubes(data, level=isoval)       # verts in grid-index units
    verts = cube["origin"] + verts @ cube["axes"]
    return verts[faces]


# ------------------------------------------------------------------ public helpers

def show_molecule(mol, interactive=True, ax=None, numbers=False, **view):
    """Ball-and-stick model of a psi4 Molecule: py3Dmol viewer, or a static picture.

    numbers=True writes each atom's number (counting from 1) beside it.
    """
    symbols, xyz = _molecule_atoms(mol)
    if interactive:
        v = py3Dmol.view(width=400, height=300)
        v.addModel(mol.save_string_xyz_file(), "xyz")
        v.setStyle({"stick": {"radius": 0.05}, "sphere": {"scale": 0.15}})     # thinner bonds (default is 0.25)
#        v.setStyle({"stick": {}, "sphere": {"scale": 0.25}}) # traditional settings
#        v.setViewStyle({"style": "outline", "color": "black", "width": 0.03}) # set outline
        if numbers:
            for n, (x, y, z) in enumerate(xyz, start=1):
                v.addLabel(str(n), {"position": {"x": x, "y": y, "z": z}, "fontColor": "black",
                                    "fontSize": 14, "backgroundOpacity": 0.0})
        v.setBackgroundColor("#f2f2f2")
        v.zoomTo()
        v.show()
        return
    ax, created = _new_3d_axes(ax, figsize=(3, 3))
    _draw_atoms(ax, symbols, xyz)
    if numbers:
        for n, (x, y, z) in enumerate(xyz, start=1):
            ax.text(x, y, z, f" {n}", fontsize=8, zorder=10)
    _finish_axes(ax, xyz, pad=0.6, **view)
    if created:
        plt.show()


def show_orbital(cube_path, mol=None, isoval=0.05, interactive=True, ax=None,
                 transparency=None, **view):
    """Orbital from a cube file, positive lobe blue and negative lobe red.

    transparency: 0 = solid lobes, 0.75 = 75 % see-through, 1 = invisible.
    """
    if interactive:
        cube_text = Path(cube_path).read_text()
        v = py3Dmol.view(width=250, height=250)
        v.addModel(mol.save_string_xyz_file() if mol is not None else cube_text,
                   "xyz" if mol is not None else "cube")
        v.setStyle({"stick": {"radius": 0.05}, "sphere": {"scale": 0.15}})     # thinner bonds (default is 0.25)
#        v.setStyle({"stick": {}, "sphere": {"scale": 0.25}}) # traditional seetings
#        v.setViewStyle({"style": "outline", "color": "black", "width": 0.03}) # set outline
        v.addVolumetricData(cube_text, "cube", {"isoval": isoval, "color": "blue",
                                                 "opacity": _opacity(transparency, 0.7)})
        v.addVolumetricData(cube_text, "cube", {"isoval": -isoval, "color": "red",
                                                 "opacity": _opacity(transparency, 0.7)})
        v.setBackgroundColor("#f2f2f2")
        v.zoomTo()
        if view.get("title"):
            print(view["title"])
        v.show()
        return
    cube = read_cube(cube_path)
    pos = isosurface(cube, isoval)
    neg = isosurface(dict(cube, data=-cube["data"]), isoval)        # negative lobe
    ax, created = _new_3d_axes(ax)
    _draw_atoms(ax, cube["symbols"], cube["xyz"], scale=0.5)
    _draw_surface(ax, pos, "#4C72B0", alpha=_opacity(transparency, 0.7))
    _draw_surface(ax, neg, "#DD4444", alpha=_opacity(transparency, 0.7))
    _finish_axes(ax, np.vstack([cube["xyz"], pos.reshape(-1, 3), neg.reshape(-1, 3)]), **view)
    if created:
        plt.show()


def show_orbitals(cube_paths, titles=None, mol=None, isoval=0.05, interactive=True,
                  transparency=None, **view):
    """Several orbitals: one static figure with a panel for each, or one viewer each."""
    cube_paths = list(cube_paths)
    titles = titles if titles is not None else [Path(p).stem for p in cube_paths]
    if interactive:
        for path, title in zip(cube_paths, titles):
            show_orbital(path, mol=mol, isoval=isoval, interactive=True,
                         transparency=transparency, title=title)
        return
    n = len(cube_paths)
    fig, axs = plt.subplots(nrows=1, ncols=n, figsize=(2.4 * n, 2.8),
                            subplot_kw={"projection": "3d"}, squeeze=False)
    for ax, path, title in zip(axs[0], cube_paths, titles):
        show_orbital(path, isoval=isoval, interactive=False, ax=ax,
                     transparency=transparency, title=title, **view)
    plt.show()


def show_density(cube_path, isoval=0.001, esp_cube_path=None, esp_range=0.05,
                 mol=None, interactive=True, ax=None, transparency=None,
                 esp_min=None, esp_max=None, **view):
    """Density isosurface; with esp_cube_path, coloured by the ESP (red negative, blue positive).

    The colour scale runs from -esp_range to +esp_range; give esp_min and esp_max
    instead for a scale that is not centred on zero (e.g. for an ion).
    transparency: 0 = solid surface, 0.75 = 75 % see-through, 1 = invisible.
    """
    esp_min = -esp_range if esp_min is None else esp_min
    esp_max = esp_range if esp_max is None else esp_max
    if interactive:
        dens = Path(cube_path).read_text()
        v = py3Dmol.view(width=400, height=350)
        v.addModel(mol.save_string_xyz_file() if mol is not None else dens,
                   "xyz" if mol is not None else "cube")
        v.setStyle({"stick": {"radius": 0.05}, "sphere": {"scale": 0.15}})     # thinner bonds (default is 0.25)
#        v.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
#        v.setViewStyle({"style": "outline", "color": "black", "width": 0.03}) # set outline
        opts = {"isoval": isoval, "color": "lightblue", "opacity": _opacity(transparency, 0.6)}
        if esp_cube_path is not None:
            opts = {"isoval": isoval, "opacity": _opacity(transparency, 0.9),
                    "voldata": Path(esp_cube_path).read_text(), "volformat": "cube",
                    "volscheme": {"gradient": "rwb", "min": esp_min, "max": esp_max}}
        v.addVolumetricData(dens, "cube", opts)
        v.setBackgroundColor("#f2f2f2")
        v.zoomTo()
        v.show()
        return
    from scipy.ndimage import map_coordinates
    cube = read_cube(cube_path)
    tris = isosurface(cube, isoval)
    ax, created = _new_3d_axes(ax)
    _draw_atoms(ax, cube["symbols"], cube["xyz"], scale=0.5)
    if esp_cube_path is None:
        _draw_surface(ax, tris, "lightblue", alpha=_opacity(transparency, 0.5))
    else:
        esp = read_cube(esp_cube_path)
        centres = tris.mean(1)
        grid_index = np.linalg.solve(esp["axes"].T, (centres - esp["origin"]).T)   # Å -> grid
        values = map_coordinates(esp["data"], grid_index, order=1)
        colours = plt.get_cmap("RdBu")(mcolors.Normalize(esp_min, esp_max)(values))
        _draw_surface(ax, tris, colours, alpha=_opacity(transparency, 0.9))
    _finish_axes(ax, np.vstack([cube["xyz"], tris.reshape(-1, 3)]), **view)
    if created:
        plt.show()


def show_frames(xyz_text, interactive=True, nshow=5, interval=500, loop="forward", **view):
    """A multi-frame XYZ string: py3Dmol animation, or a static strip of nshow frames."""
    if interactive:
        v = py3Dmol.view(width=400, height=300)
        v.addModelsAsFrames(xyz_text, "xyz")
        v.setStyle({"stick": {"radius": 0.05}, "sphere": {"scale": 0.15}})     # thinner bonds (default is 0.25)
#        v.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
#        v.setViewStyle({"style": "outline", "color": "black", "width": 0.03}) # set outline
        v.setBackgroundColor("#f2f2f2")
        v.zoomTo()
        v.animate({"loop": loop, "interval": interval})
        v.show()
        return
    frames = _xyz_frames(xyz_text)
    picks = np.unique(np.linspace(0, len(frames) - 1, min(nshow, len(frames))).round().astype(int))
    fig, axs = plt.subplots(nrows=1, ncols=len(picks), figsize=(2.2 * len(picks), 2.4),
                            subplot_kw={"projection": "3d"}, squeeze=False)
    for ax, k in zip(axs[0], picks):
        _, symbols, xyz = frames[k]
        _draw_atoms(ax, symbols, xyz, scale=0.6)
        _finish_axes(ax, xyz, pad=0.6, **{"title": f"frame {k + 1} of {len(frames)}", **view})
    plt.show()


def show_mode(mol, vectors, interactive=True, amplitude=1.0, ax=None, **view):
    """A normal mode: py3Dmol vibration, or a static picture with displacement arrows.

    vectors: (natom, 3) displacement of each atom for this mode.
    """
    symbols, xyz = _molecule_atoms(mol)
    vectors = np.asarray(vectors).reshape(-1, 3)
    if interactive:
        lines = [str(len(symbols)), "mode"]
        for s, (x, y, z), (dx, dy, dz) in zip(symbols, xyz, vectors):
            lines.append(f"{s} {x:.6f} {y:.6f} {z:.6f} {dx:.6f} {dy:.6f} {dz:.6f}")
        v = py3Dmol.view(width=300, height=250)
        v.addModel("\n".join(lines), "xyz", {"vibrate": {"frames": 10, "amplitude": amplitude}})
        v.setStyle({"stick": {"radius": 0.05}, "sphere": {"scale": 0.15}})     # thinner bonds (default is 0.25)
#        v.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
#        v.setViewStyle({"style": "outline", "color": "black", "width": 0.03}) # set outline
        v.setBackgroundColor("#f2f2f2")
        v.zoomTo()
        v.animate({"loop": "backAndForth"})
        v.show()
        return
    ax, created = _new_3d_axes(ax, figsize=(3, 3))
    _draw_atoms(ax, symbols, xyz, scale=0.6)
    arrows = vectors / np.abs(vectors).max() * 0.6 * amplitude       # longest arrow 0.6 Å
    ax.quiver(*xyz.T, *arrows.T, color="black", linewidth=1.2, arrow_length_ratio=0.3)
    _finish_axes(ax, np.vstack([xyz, xyz + arrows]), pad=0.6, **view)
    if created:
        plt.show()


def show_vibration(wfn, mode=0, interactive=True, amplitude=1.0, **view):
    """Normal mode number `mode` (0 = lowest frequency) of a psi4.frequency wavefunction."""
    vib = wfn.frequency_analysis
    is_vib = np.array([t == "V" for t in vib["TRV"].data])          # skip translations/rotations
    vectors = vib["x"].data[:, is_vib][:, mode].reshape(-1, 3)
    show_mode(wfn.molecule(), vectors, interactive=interactive, amplitude=amplitude, **view)


def bond_length(mol, i, j):
    """Return the distance between atoms i and j in angstroms (0-based indices)."""
    i = i-1; j = j-1  # convert to 0-based indices
    xyz = mol.geometry().np * psi4.constants.bohr2angstroms
    return np.linalg.norm(xyz[i] - xyz[j])


def bond_angle(mol, i, j, k):
    """Return the i-j-k angle in degrees, with j as the central atom (0-based indices)."""
    i = i-1; j = j-1; k = k-1  # convert to 0-based indices
    xyz = mol.geometry().np
    a = xyz[i] - xyz[j]
    b = xyz[k] - xyz[j]
    cos_angle = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return np.degrees(np.arccos(cos_angle))


def dihedral(mol, i, j, k, l):
    """Return the i-j-k-l torsion angle in degrees, looking down the j-k bond (atoms counted from 1)."""
    xyz = mol.geometry().np[[i - 1, j - 1, k - 1, l - 1]]
    b1, b2, b3 = xyz[1] - xyz[0], xyz[2] - xyz[1], xyz[3] - xyz[2]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    x = np.dot(n1, n2)
    y = np.dot(np.cross(n1, n2), b2) / np.linalg.norm(b2)
    return np.degrees(np.arctan2(y, x))


def history_to_xyz(mol, history):
    """Build a multi-frame XYZ string from an optimization history."""
    symbols = [mol.symbol(i) for i in range(mol.natom())]
    frames = []
    for step, coords in enumerate(history["coordinates"]):
        xyz = np.asarray(coords) * psi4.constants.bohr2angstroms
        lines = [str(len(symbols)), f"step {step + 1}"]
        lines += [f"{s} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(symbols, xyz)]
        frames.append("\n".join(lines))
    return "\n".join(frames) + "\n"


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


def new_output(filename):
    """Send Psi4 output to filename (and its .log), first detaching the .log of any earlier job.

    psi4.set_output_file() alone leaves old .log files attached, so every later job
    would also write its log into them.
    """
    import logging
    logger = logging.getLogger("psi4")
    for handler in logger.handlers[:]:
        if isinstance(handler, logging.FileHandler):
            handler.close()
            logger.removeHandler(handler)
    psi4.set_output_file(str(filename), False)


def c1_copy(mol):
    """A new psi4 Molecule with the same atoms and charge, but no symmetry and no reorientation."""
    xyz = mol.geometry().np * BOHR2ANG
    lines = [f"{int(mol.molecular_charge())} {mol.multiplicity()}"]
    lines += [f"{mol.symbol(i)} {x:.10f} {y:.10f} {z:.10f}" for i, (x, y, z) in enumerate(xyz)]
    lines += ["symmetry c1", "no_reorient", "no_com"]
    return psi4.geometry("\n".join(lines))


def read_irc(filename):
    """Read an optking IRC trajectory file (coordinates in bohr) into a list of frames.

    Each frame is a dict: "Step" (int), "Energy" (hartree), "xyz" (natom x 3 array, Å)
    and "text" (the frame as XYZ text in Å, ready for show_frames).
    """
    lines = Path(filename).read_text().splitlines()
    frames, i = [], 0
    while i < len(lines) and lines[i].strip():
        natom = int(lines[i])
        title = lines[i + 1].split()                     # IRC trajectory step: 1. E: -232.8...
        atoms = [l.split() for l in lines[i + 2:i + 2 + natom]]
        symbols = [a[0] for a in atoms]
        xyz = np.array([a[1:4] for a in atoms], float) * BOHR2ANG
        text = f"{natom}\n{lines[i + 1].strip()}\n"
        text += "".join(f"{s} {x:10.6f} {y:10.6f} {z:10.6f}\n" for s, (x, y, z) in zip(symbols, xyz))
        frames.append({"Step": int(float(title[3])), "Energy": float(title[5]),
                       "text": text, "xyz": xyz})
        i += 2 + natom
    return frames
