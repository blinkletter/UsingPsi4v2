# The Helper Functions

Every notebook in this book starts by loading a small file of our own code, `psi4_helpers.py`, with the line `import psi4_helpers as ph`. The functions in it measure bonds and angles, plot an optimization, search output files and, most of all, **draw molecules, orbitals and surfaces**. Putting them in one file keeps the notebooks short and makes every picture look the same.

You do not need to understand this code to use the notebooks. But it is ordinary *Python*, and reading it is a good way to see how a cube file becomes an orbital picture. This page walks through the file section by section. The code shown here is read straight from `psi4_helpers.py`, so it is always the current version.

## What You Will Learn
- Which helper functions exist and which notebooks use them.
- How the `INTERACTIVE` switch chooses between a fixed picture and a 3D model you can rotate.
- What is inside a **cube file**, and how an **isosurface** is built from it.
- How to change the view, the transparency and the isovalue of a picture.

## Quick Reference

| Function | What it does | Used in |
|---|---|---|
| `ph.bond_length(mol, i, j)` | Distance between atoms `i` and `j` (Å) | *Structure*, *Optimization*, *Level of Effort*, *Bending Scan*, *Hydrogen Bonding* |
| `ph.bond_angle(mol, i, j, k)` | Angle `i`–`j`–`k` (degrees), `j` in the middle | *Structure*, *Optimization*, *Level of Effort*, *Hydrogen Bonding* |
| `ph.plot_history(history, title)` | Energy at each step of an optimization | *Optimization* |
| `ph.search_file(filename, text)` | Print the lines of a file that contain `text` | *Structure*, *Thermochemistry* |
| `ph.show_molecule(mol)` | Ball-and-stick model | *Structure*, *Hydrogen Bonding* |
| `ph.show_frames(xyz_text)` | Several structures: an animation, or a strip of pictures | *Optimization*, *Bending Scan* |
| `ph.show_mode(mol, vectors)` | A vibration: an animation, or arrows on the atoms | *Vibrations* |
| `ph.show_orbital(cube_path)` | One orbital from a cube file | *Hydrogen Bonding* |
| `ph.show_orbitals(cube_paths)` | Several orbitals side by side | *Molecular Orbitals* |
| `ph.show_density(cube_path)` | Electron density surface, optionally coloured by the ESP | *Surfaces* |

Atoms are numbered from 0, as everywhere in *Python*: in water built as `O`, `H`, `H`, the oxygen is atom 0.

## 1. Loading the File

The first cell of every notebook makes the helpers available. On your own computer, `psi4_helpers.py` sits one folder above the notebooks, so we add that folder to the places *Python* looks for modules. On *Colab* the file is downloaded from the book's *GitHub* repository. In outline (the real cell does a little more):

```python
if IN_COLAB:
    !wget -q -N https://raw.githubusercontent.com/blinkletter/UsingPsi4v2/main/psi4_helpers.py
else:
    sys.path.insert(0, "..")                    # psi4_helpers.py is one folder up

import psi4_helpers as ph

# 3D pictures: False = static images (as in the book); True = rotatable py3Dmol viewers
INTERACTIVE = False
```

The top of the file imports the libraries it needs and sets up two small tables: the colour and **covalent radius** of each element, used to draw atoms and decide which atoms are bonded, and the element symbol for each atomic number, used when reading cube files.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: from pathlib import Path
:end-before: def _new_3d_axes
```

*Psi4* works in **bohr** (atomic units of length) internally. `BOHR2ANG` converts to ångström, which is what chemists usually quote.

## 2. Measuring Structures

`bond_length()` takes the coordinates from the molecule, converts them to ångström and returns the length of the vector between two atoms. `bond_angle()` uses the **dot product** of the two bond vectors, $\cos\theta = \dfrac{\mathbf{a}\cdot\mathbf{b}}{|\mathbf{a}|\,|\mathbf{b}|}$. The angle does not depend on units, so no conversion is needed.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def bond_length(
:end-before: def plot_history(
```

## 3. Plotting an Optimization and Searching Output

`plot_history()` takes the `history` dictionary returned by `psi4.optimize(..., return_history=True)` and plots the energy of each step relative to the final one, in kJ/mol, in the house style of the book. `search_file()` prints every line of a text file that contains a word, with its line number. It is the quickest way to find a result in a long *Psi4* output file.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def plot_history(
```

## 4. Two Ways to Draw: the `INTERACTIVE` Switch

Every drawing helper has an `interactive` argument and the same overall shape:

```python
def show_something(..., interactive=True, ...):
    if interactive:
        # build a py3Dmol viewer and show it
        return
    # otherwise, draw a fixed picture with matplotlib
```

Why have two versions?

- A ***py3Dmol*** viewer runs in your web browser. You can rotate and zoom it, which is the best way to explore a molecule. But it only works in a live notebook, and to draw a surface it stores the **whole cube file** inside the notebook. With a dozen orbitals that adds up to tens of megabytes; the *Molecular Orbitals* notebook was once 72 MB.
- A ***matplotlib*** picture is an ordinary image of a few tens of kilobytes. It appears in the book exactly as it was drawn, but you cannot turn it.

The notebooks pass `interactive=INTERACTIVE` to every drawing helper, so one line in the first cell switches the whole notebook. `INTERACTIVE = False` gives the pictures you see in this book. Set it to `True` and rerun the notebook to get models you can rotate. (Called without the argument, the helpers draw the interactive version, as they did in earlier versions of this book.)

## 5. Drawing a Molecule

`show_molecule()` is the simplest drawing helper. The interactive version hands the molecule to *py3Dmol* as an $xyz$ string. The static version gets the element symbols and coordinates, then calls three small tools that all the static pictures share.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def show_molecule(
:end-before: def show_orbital(
```

The tools have names that start with an underscore (`_draw_atoms`), a *Python* convention for "used inside this file; you do not need to call it yourself".

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def _new_3d_axes(
:end-before: def _shaded(
```

- `_new_3d_axes()` makes a new 3D figure, unless we pass in an existing set of axes with `ax=`. That is how `show_orbitals()` puts several orbitals side by side in one figure.
- `_finish_axes()` makes the $x$, $y$ and $z$ scales equal (otherwise molecules come out squashed), places the camera and hides the axes. `elev` is the camera's height above the $xy$-plane and `azim` its angle around the $z$ axis, both in degrees; `zoom` makes the picture fill more of its box.
- `_draw_atoms()` draws a sphere for each atom and a line between two atoms when they are closer than 1.2 times the sum of their covalent radii. There is no bonding information in an $xyz$ file, so this distance test is how almost every molecule viewer decides where to put the bonds.

## 6. Cube Files and Isosurfaces

An orbital or an electron density is a value at **every point in space**. *Psi4* writes such a function to a **cube file** (`psi4.cubeprop()`): the values on a regular three-dimensional grid of points, typically 0.2 bohr apart. The file starts with a short header: two comment lines, the number of atoms and the corner (origin) of the grid, the number of points and the step along each axis, then one line per atom. After that come the values, many thousands of them.

`read_cube()` reads that layout into a dictionary, converting every distance to ångström.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def read_cube(
:end-before: def isosurface(
```

To draw a function of three dimensions we draw a **contour**: the surface on which the function has one chosen value, the **isovalue**. This surface is an **isosurface**, the 3D version of a contour line on a map. The standard way to find it is the **marching cubes** algorithm: walk through every little cube of eight neighbouring grid points, find where the value crosses the isovalue along its edges, and join those crossing points into triangles. We use the version in the *scikit-image* library.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def isosurface(
:end-before: def show_molecule(
```

`marching_cubes` returns the corners of the triangles (`verts`) in grid units, and which three corners make each triangle (`faces`). The line `cube["origin"] + verts @ cube["axes"]` converts grid units into positions in ångström: start at the corner of the grid and step along each axis.

## 7. Drawing Orbitals

An orbital has a **sign**. `show_orbital()` draws two isosurfaces: one where the orbital equals `+isoval` (blue) and one where it equals `-isoval` (red). The negative lobe is found by flipping the sign of every grid value and asking for `+isoval` again.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def show_orbital(
:end-before: def show_density(
```

`show_orbitals()` draws a list of cube files. In static mode it makes one figure with a panel for each orbital and hands each panel to `show_orbital()` through `ax=`.

A flat-coloured surface looks like a cut-out shape, so `_shaded()` adds simple lighting. For each triangle it finds the **normal** (the direction the triangle faces, from a cross product) and darkens triangles that face away from an imaginary light. `_opacity()` turns the `transparency` argument into the opacity that *matplotlib* and *py3Dmol* expect: transparency 0.75 means opacity 0.25. `_draw_surface()` then adds the triangles to the picture.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def _shaded(
:end-before: def _molecule_atoms(
```

```{note} The isovalue
:class: dropdown
The default `isoval=0.05` (in atomic units) draws the region where an orbital is fairly large. A **smaller** isovalue draws a bigger, puffier surface further out from the nuclei; a **larger** one shrinks the lobes towards the atoms. Comparing orbitals only makes sense at the same isovalue. *Psi4* writes a suggested value for each orbital in the second line of the cube file (the contour that encloses 85 % of the orbital density).
```

## 8. Drawing the Density and the Electrostatic Potential

`show_density()` draws the electron density at `isoval=0.001`, the usual molecular surface (see the *Surfaces* notebook). With `esp_cube_path` it also colours the surface by the **electrostatic potential**.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def show_density(
:end-before: def show_frames(
```

The ESP cube holds the potential at the same grid points as the density. To colour a triangle we need the potential at its centre, which usually falls between grid points. `map_coordinates()` from *SciPy* **interpolates**: it estimates the value from the surrounding grid points. Each value is then turned into a colour on a red–white–blue scale from `-esp_range` (red, negative) to `+esp_range` (blue, positive).

## 9. Animations, or Pictures of Motion

Two helpers draw something that moves.

`show_frames()` takes a **multi-frame** $xyz$ string: several structures, one after another, each starting with its number of atoms and a comment line. The interactive version plays them as an animation. The static version picks `nshow` frames, evenly spaced from first to last, and draws them in a row.

`show_mode()` takes a molecule and the **displacement vectors** of one normal mode (how far, and in which direction, each atom moves). *py3Dmol* can animate a vibration directly from these vectors. The static version draws them as arrows, scaled so the longest arrow is 0.6 Å.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def show_frames(
:end-before: def bond_length(
```

`show_frames()` uses one more small tool, `_xyz_frames()`, to split the multi-frame string into separate structures.

```{literalinclude} psi4_helpers.py
:lang: python
:start-at: def _molecule_atoms(
:end-before: def read_cube(
```

## 10. Changing a Picture

All the static drawing helpers accept these keyword arguments.

| Keyword | Default | Effect | Helpers |
|---|---|---|---|
| `elev` | 20 | Camera height above the $xy$-plane (degrees) | all |
| `azim` | 35 | Camera angle around the $z$ axis (degrees) | all |
| `zoom` | 1.4 | Larger values fill more of the panel | all |
| `isoval` | 0.05 (orbitals), 0.001 (density) | Value at which the surface is drawn | `show_orbital(s)`, `show_density` |
| `transparency` | set by each helper | 0 = solid, 1 = invisible; the notebooks use 0.75 (0.4 for the ESP map) | `show_orbital(s)`, `show_density` |
| `esp_range` | 0.05 | Potential (hartree/e) at the ends of the colour scale | `show_density` |
| `nshow` | 5 | Number of frames drawn | `show_frames` |
| `amplitude` | 1.0 | Arrow length (static) or size of the motion (animated) | `show_mode` |

*Psi4* places a small molecule in a standard orientation that follows its symmetry. Water lies in the $yz$-plane with its $C_2$ axis along $z$, so `elev=0, azim=0` looks straight at the molecular plane and `elev=0, azim=90` looks along the plane, edge-on. The edge-on view is the one that shows out-of-plane lone pairs.

```python
ph.show_orbitals(lmo_paths, mol=water, interactive=INTERACTIVE,
                 transparency=0.75, elev=0, azim=90)
```

```{tip} Try it
- In the *Molecular Orbitals* notebook, draw the Boys orbitals with `elev=0, azim=90` and then `azim=0`. In which view are the two lone pairs easiest to tell apart?
- Draw the $1b_1$ orbital with `isoval=0.02` and `isoval=0.1`. How much does the apparent size of the orbital depend on this choice?
- In the *Surfaces* notebook, change `esp_range` from 0.05 to 0.02. What happens to the colours?
```

## The Whole File

````{note} psi4_helpers.py
:class: dropdown
```{literalinclude} psi4_helpers.py
:lang: python
```
````

## Summary
- `psi4_helpers.py` holds the shared code; the notebooks load it with `import psi4_helpers as ph`.
- Every drawing helper can make a rotatable *py3Dmol* model or a small fixed *matplotlib* picture; `INTERACTIVE` in the first cell chooses.
- A cube file holds a function on a 3D grid. An isosurface, found by **marching cubes**, is the surface where that function has one chosen value.
- Orbitals are drawn as two isosurfaces, $\pm$`isoval`; the ESP is drawn by colouring the density surface.
- `elev`, `azim`, `zoom`, `isoval` and `transparency` change how a picture looks.
