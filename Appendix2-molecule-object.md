# The Molecule Object

When we build a molecule with `psi4.geometry()`, *Psi4* hands back a **Molecule object**:

```python
water = psi4.geometry("""
0 1
O
H 1 r
H 1 r 2 a
r = 0.96
a = 104.5
""", name="water")
```

`water` is more than the text we typed. It knows its atoms, masses, charge, symmetry and coordinates, and it carries **methods** (functions attached to the object) that report or change those things. Typing `dir(water)` lists them all: 156 names in *Psi4* 1.11. The 28 that begin with an underscore (`__init__`, `__repr__`, ...) are ordinary *Python* plumbing, and we can ignore them. This page explains the other 128, grouped by what they are for.

## What You Will Learn
- How to ask a Molecule object about its atoms, coordinates, symmetry and fragments.
- Which methods change a molecule, and how to make a safe copy first.
- How to save a structure as an XYZ file or a string.
- The traps: atoms numbered from 0, fragments from 1, and coordinates in bohr.

## How to Read This Page

- **Call methods with parentheses.** `water.natom()` gives `3`; `water.natom` without parentheses only names the method.
- **Atoms are numbered from 0.** In water built as `O`, `H`, `H`, the oxygen is atom 0 and the hydrogens are atoms 1 and 2. **Fragments** (the pieces between `--` lines) are numbered from **1**.
- **Coordinates are in bohr** (atomic units), even when the geometry was typed in ångströms. Multiply by `psi4.constants.bohr2angstroms` (0.529177...) to get ångströms.
- **Run `water.update_geometry()` first** if the molecule has not yet been used in a calculation. Until then, a Z-matrix has not been turned into coordinates and some methods return nothing useful.
- **Methods that start with `f`** (`fx`, `fsymbol`, `full_geometry`, ...) count **dummy atoms** too (the `X` atoms sometimes used in Z-matrices). For molecules without dummy atoms, the `f` and plain versions give the same answer.
- **Methods that start with `print_`** write to the *Psi4* output file, not to the notebook. Use `ph.search_file()` or open the `.out` file to see them.
- Some methods return *Psi4*'s own types: a `Matrix`, a `Vector` or a `Vector3`. Add `.np` to a `Matrix` or `Vector` to get a NumPy array. Read a `Vector3` one item at a time (`v[0]`, `v[1]`, `v[2]`) or with `print(v)`.

```{warning} Do not call list() on a Vector3
`list(water.xyz(1))` or `list(water.center_of_mass())` never finishes and the notebook cell hangs. Use `[v[0], v[1], v[2]]` or `np.array(water.geometry())[1]` instead.
```

All example output below is for water with $r$ = 0.96 Å and $\angle$HOH = 104.5°, as defined above.

## 1. Atoms: Counting and Identifying

| Method | What it gives | Water example |
|---|---|---|
| `natom()` | Number of real atoms | `3` |
| `nallatom()` | Number of atoms including dummies | `3` |
| `symbol(i)` | Element symbol of atom `i` | `water.symbol(0)` → `'O'` |
| `label(i)` | Label as typed in the input (`H1`, `C2`, ...) | `'O'` |
| `Z(i)` | Nuclear charge used in the calculation (0 for a ghost atom) | `8.0` |
| `true_atomic_number(i)` | Atomic number of the element, even for a ghost atom | `8` |
| `charge(i)` | Nuclear charge of atom `i` (same as `Z` unless an ECP changes it) | `8.0` |
| `mass(i)` | Mass in amu (most abundant isotope by default) | `15.9949` |
| `mass_number(i)` | Mass number *A*, or `-1` if *Psi4* doesn't know it | `16` |
| `basis_on_atom(i)` | Name of the basis set on atom `i` (only after a calculation has assigned one) | `'6-31G*'` |
| `fZ(i)`, `fcharge(i)`, `flabel(i)`, `fmass(i)`, `fsymbol(i)`, `ftrue_atomic_number(i)` | The same as above, but counting dummy atoms | |

A loop over the atoms is the usual way to use these:

```python
for i in range(water.natom()):
    print(i, water.symbol(i), water.mass(i))
```

## 2. Coordinates and Geometry

| Method | What it gives | Water example |
|---|---|---|
| `geometry()` | All coordinates as an `natom × 3` Matrix, in bohr | `water.geometry().np` |
| `full_geometry()` | The same, including dummy atoms | |
| `x(i)`, `y(i)`, `z(i)` | One coordinate of atom `i`, in bohr | `water.y(1)` → `-1.434` |
| `fx(i)`, `fy(i)`, `fz(i)` | The same, counting dummy atoms | |
| `xyz(i)` | All three coordinates of atom `i` as a `Vector3` | `print(water.xyz(1))` |
| `distance_matrix()` | Every atom–atom distance as a Matrix, in bohr | O–H 0.96 Å, H···H 1.518 Å after conversion |
| `atom_at_position(coord, tol)` | Index of the atom within `tol` bohr of `coord`, or `-1` | |
| `center_of_mass()` | Centre of mass as a `Vector3` (does not move the molecule) | `[0, 0, 0]` |
| `inertia_tensor()` | Moment of inertia tensor (amu·bohr²) as a Matrix | |
| `units()` | Units the geometry was *typed* in | `'Angstrom'` |
| `input_units_to_au()` | Factor that turns input units into bohr | `1.88973` |

Coordinates in ångströms as a NumPy array:

```python
xyz_angstrom = water.geometry().np * psi4.constants.bohr2angstroms
```

The bond lengths in ångströms (the helper `ph.bond_length()` does the same job):

```python
dist = water.distance_matrix().np * psi4.constants.bohr2angstroms
print(f"O–H = {dist[0, 1]:.3f} Å,  H···H = {dist[1, 2]:.3f} Å")
```

## 3. Properties of the Whole Molecule

| Method | What it gives | Water example |
|---|---|---|
| `name()` | The molecule's name (set with `name=` in `psi4.geometry`) | `'water'` |
| `molecular_charge()` | Overall charge | `0` |
| `multiplicity()` | Spin multiplicity, 2*S* + 1 | `1` |
| `comment()` | A free-text comment, empty unless set | `''` |
| `provenance()` | Which program and version created the object | a dictionary |
| `connectivity()` | Bonds, if any were supplied (normally empty; *Psi4* doesn't need them) | `[]` |
| `nuclear_repulsion_energy()` | Repulsion between the nuclei, in hartree | `9.168` |
| `nuclear_repulsion_energy_deriv1()` | Its gradient, `natom × 3` | |
| `nuclear_repulsion_energy_deriv2()` | Its second derivatives, `3natom × 3natom` | |
| `nuclear_dipole()` | Nuclear part of the dipole moment (a.u.); the electrons add the rest | `[0, 0, 0.978]` |
| `rotational_constants()` | Rotational constants *A*, *B*, *C* in cm⁻¹ | `[27.26, 14.52, 9.47]` |
| `rotor_type()` | Linear, symmetric top, asymmetric top, ... | `'RT_ASYMMETRIC_TOP'` |
| `rotational_symmetry_number()` | σ, used in the rotational entropy (*Thermochemistry*) | `2` |

The **nuclear repulsion energy** is a useful fingerprint: two geometries that give the same value to many decimal places are the same structure.

## 4. Symmetry

*Psi4* uses two point groups. The **full point group** is the true symmetry of the molecule. The **computational point group** is the largest **Abelian** subgroup ($D_{2h}$ or one of its subgroups), which is what the program actually uses to speed up the calculation. For water both are $C_{2v}$; for methane the full group is $T_d$ but the computational group is a subgroup such as $C_{2v}$.

| Method | What it gives | Water example |
|---|---|---|
| `get_full_point_group()` | Full point group | `'C2v'` |
| `get_full_point_group_with_n()` | The same with *n* left as a letter | `'Cnv'` |
| `full_pg_n()` | The *n* in that name | `2` |
| `schoenflies_symbol()` | Computational point group | `'c2v'` |
| `point_group()` | Computational point group as an object; `.symbol()` gives its name | `water.point_group().symbol()` |
| `irrep_labels()` | Labels of the irreducible representations, used for vibrations and orbitals | `['A1', 'A2', 'B1', 'B2']` |
| `symmetry_from_input()` | Point group typed with `symmetry` in the input, if any | `''` |
| `find_point_group(tol)` | Works out the computational point group | |
| `find_highest_point_group(tol)` | Highest possible computational point group | |
| `symmetrize(tol)` | Nudges the atoms so the structure has exact symmetry (within `tol` bohr) | |
| `set_point_group(pg)`, `reset_point_group(name)`, `form_symmetry_information(tol)` | Set the point group by hand. Expert use; the `symmetry` line in the geometry string is easier. | |

## 5. Z-Matrix Variables

A geometry written with variables such as `r` and `a` keeps those variables, and we can read or change them. This is how a scan (see *Bending Scan*) can change one angle at a time.

| Method | What it does | Water example |
|---|---|---|
| `has_zmatrix()` | `True` if any atom is given by a Z-matrix line | `True` |
| `is_variable(name)` | `True` if `name` is one of the variables | `water.is_variable("r")` → `True` |
| `get_variable(name)` | Current value of a variable (Å or degrees) | `water.get_variable("a")` → `104.5` |
| `set_variable(name, value)` | Changes a variable and rebuilds the geometry | `water.set_variable("a", 120.0)` |
| `update_geometry()` | Rebuilds the coordinates from the input and variables. Safe to call any time. | |
| `reinterpret_coordentry(flag)` | Whether `update_geometry()` should rebuild from the Z-matrix (`True`, the default) or keep the current Cartesians | |

## 6. Changing the Molecule

Most of these change the molecule **in place**. Make a copy with `clone()` first if you still need the original.

| Method | What it does |
|---|---|
| `clone()` | Returns an independent copy: `d2o = water.clone()` |
| `set_name(text)`, `set_comment(text)` | Set the name or comment |
| `set_molecular_charge(q)`, `set_multiplicity(m)` | Change charge or multiplicity. Better to set them on the first line of the geometry string. |
| `set_mass(i, mass)` | Change the mass of atom `i` in amu. Used for isotopes (*Isotope Effects*): `d2o.set_mass(1, 2.014101778)` |
| `set_nuclear_charge(i, Z)` | Change the nuclear charge of atom `i` (for ECPs; expert use) |
| `set_geometry(matrix)` | Replace all coordinates with an `natom × 3` Matrix in bohr |
| `set_full_geometry(matrix)` | The same, including dummy atoms |
| `translate(vector3)` | Move every atom by a `Vector3` in bohr |
| `move_to_com()` | Move the molecule so its centre of mass is at the origin |
| `fix_com(flag)`, `fix_orientation(flag)` | Stop *Psi4* moving or rotating the molecule. Same as `no_com` and `no_reorient` in the geometry string, which is the easier way. |
| `com_fixed()`, `orientation_fixed()` | Report whether those are switched on |
| `set_units(units)`, `set_input_units_to_au(factor)` | Change the input units. Expert use. |
| `add_atom(Z, x, y, z, symbol, mass, charge, label, A)` | Add an atom by hand (coordinates in bohr). Expert use; easier to write a new geometry string. |
| `set_connectivity(bonds)`, `set_provenance(dict)` | Store bonds or provenance information |

To give a molecule new coordinates from a NumPy array:

```python
xyz = np.array(water.geometry())        # bohr
xyz[1, 1] += 0.1                        # move H1 by 0.1 bohr along y
water.set_geometry(psi4.core.Matrix.from_array(xyz))
water.update_geometry()
```

## 7. Printing to the Output File

These write a report into the *Psi4* output file (`.out`). Nothing appears in the notebook.

| Method | What it writes |
|---|---|
| `print_out()` | Cartesian coordinates in the input units, with charge, multiplicity and masses |
| `print_out_in_angstrom()`, `print_out_in_bohr()` | The same in the stated units |
| `print_in_input_format()` | The geometry as it was typed (Cartesian or Z-matrix) |
| `print_cluster()` | Coordinates with `--` lines between fragments |
| `print_distances()` | Every interatomic distance (Å) |
| `print_bond_angles()` | Every bond angle (degrees) |
| `print_out_of_planes()` | Out-of-plane angles (molecules with four or more atoms) |
| `print_rotational_constants()` | Rotational constants in cm⁻¹ and MHz |

For water, `print_distances()` and `print_bond_angles()` add these lines to the output file:

```text
        Interatomic Distances (Angstroms)

        Distance 1 to 2 0.960
        Distance 1 to 3 0.960
        Distance 2 to 3 1.518

        Bond Angles (degrees)

        Angle 2-1-3:  104.500
        Angle 1-2-3:   37.750
        Angle 1-3-2:   37.750
```

Notice that the printout numbers atoms from **1**, unlike the methods themselves.

## 8. Saving and Converting

| Method | What it gives |
|---|---|
| `save_string_xyz_file()` | A standard XYZ file as a string (atom count, blank line, coordinates in Å). The book uses this to draw molecules. |
| `save_string_xyz()` | Charge and multiplicity line followed by coordinates in Å, ready to paste into `psi4.geometry()` |
| `save_xyz_file(filename, True)` | Writes an XYZ file to disk (`True` = ångströms) |
| `to_string(dtype=...)` | The structure in another format: `"xyz"`, `"psi4"`, `"orca"`, `"nwchem"`, `"gamess"`, ... Add `units="Angstrom"` for formats that would otherwise use bohr. |
| `create_psi4_string_from_molecule()` | A *Psi4* geometry string for the current state, Z-matrix variables included |
| `format_molecule_for_mol()` | An MDL molfile (`.mol`) block, which RDKit and many viewers can read |
| `to_arrays()` | NumPy arrays: coordinates (bohr), masses, element symbols, nuclear charges, uncertainties |
| `to_dict()` | Everything about the molecule as a *Python* dictionary |
| `to_schema(dtype=2)` | A dictionary in the QCSchema format used by other quantum chemistry tools |

`from_string`, `from_arrays`, `from_dict` and `from_schema` go the other way: they **build a new molecule**. They are called on the class rather than on an existing molecule:

```python
helium = psi4.core.Molecule.from_string("He 0 0 0")
```

For water, `save_string_xyz_file()` gives:

```text
3

 O    0.000000000000    0.000000000000   -0.065775570547
 H    0.000000000000   -0.759061990794    0.521953018286
 H    0.000000000000    0.759061990794    0.521953018286
```

The coordinates are not the ones we typed: *Psi4* has moved the centre of mass to the origin and lined up the molecule with its symmetry axes.

## 9. Fragments

A geometry string can be split into **fragments** with `--` lines, for example the two molecules of a water dimer (*Hydrogen Bonding*). Fragments are numbered from **1**. A fragment can be **Real** (a normal part of the calculation), **Ghost** (its basis functions are kept but its atoms have no nuclei or electrons, as in a counterpoise correction) or **Absent** (left out).

| Method | What it does | Water dimer example |
|---|---|---|
| `nfragments()` | Number of fragments | `2` |
| `get_fragments()` | Atom range of each fragment, as (first, last + 1), counting atoms from 0 | `[(0, 3), (3, 6)]` |
| `get_fragment_types()` | Real, Ghost or Absent for each fragment | `['Real', 'Real']` |
| `get_fragment_charges()`, `get_fragment_multiplicities()` | Charge and multiplicity of each fragment | `[0, 0]`, `[1, 1]` |
| `extract_subsets(real, ghost)` | A **new** molecule with the chosen fragments Real, Ghost or left out | `dimer.extract_subsets(1)` → first water alone; `dimer.extract_subsets(2, 1)` → second water with the first as ghost |
| `set_active_fragment(n)`, `set_active_fragments([...])` | Make fragments Real | |
| `set_ghost_fragment(n)`, `set_ghost_fragments([...])` | Make fragments Ghost | |
| `activate_all_fragments()`, `deactivate_all_fragments()` | Make every fragment Real, or none | |
| `BFS()` | Finds separate molecules by looking for bonded atoms (atoms counted from 0) | `([[0, 1, 2], [3, 4, 5]],)` |

For counterpoise corrections, `psi4.energy(..., bsse_type="cp")` does all of this for us; these methods are for when we want to do it by hand.

## 10. Basis Sets on Particular Atoms

Normally the `basis` option sets one basis set for every atom. These methods mix basis sets within one molecule, for example a larger basis on the atoms that matter. The `basis` block in an input is the usual way to do the same thing.

| Method | What it does |
|---|---|
| `set_basis_all_atoms(name, "BASIS")` | One basis set on every atom |
| `set_basis_by_symbol(symbol, name, "BASIS")` | A basis set on every atom of one element: `water.set_basis_by_symbol("O", "aug-cc-pvdz", "BASIS")` |
| `set_basis_by_label(label, name, "BASIS")` | A basis set on atoms with a given label (`H1`) |
| `set_basis_by_number(i, name, "BASIS")` | A basis set on one atom |
| `basis_on_atom(i)` | Reports the basis set on atom `i` (fails if none has been assigned yet) |

## 11. Specialist Tools

These are for advanced work and are not needed in this book.

| Method | What it does |
|---|---|
| `run_dftd3(func)`, `run_sdftd3(func)`, `run_dftd4(func)` | Grimme's dispersion correction (D3 or D4) for this geometry. *Psi4* calls these itself when we ask for a method such as `b3lyp-d3bj`. |
| `run_gcp(func)` | Grimme's geometric counterpoise correction for basis set superposition error |
| `B787(ref_mol)` | Lines this molecule up on top of `ref_mol` (shift, rotation and atom reordering) and reports the RMSD. Useful for checking whether two structures are the same. |
| `scramble(ref_mol)` | A test for `B787`: randomly moves and reorders a molecule, then checks that `B787` undoes it |

## Summary
- `psi4.geometry()` returns a Molecule object; `dir()` lists its 128 public methods.
- Ask questions with methods such as `natom()`, `symbol(i)`, `geometry()`, `get_full_point_group()` and `rotational_constants()`.
- Atoms count from 0, fragments from 1, and coordinates come back in bohr.
- `clone()` before changing a molecule with `set_mass()`, `set_variable()` or `set_geometry()`.
- `save_string_xyz_file()` and `to_string()` turn a molecule into text for other programs.

```{warning} AI Statement
:class: dropdown
This page was written by Claude (Anthropic) using Claude Code, at my request. Claude listed the methods with `dir()`, read their descriptions in *Psi4* 1.11, and ran each example on water and a water dimer to check the output shown. I selected, edited and checked all content.
```
