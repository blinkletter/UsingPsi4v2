# The Wavefunction Object


When we add `return_wfn=True` to a calculation, *Psi4* hands back a **Wavefunction object** as well as the energy:

```python
water = psi4.geometry("""
0 1
O
H 1 0.96
H 1 0.96 2 104.5
""", name="water")

psi4.optimize("hf/6-31g*", molecule=water)
E, wfn = psi4.energy("hf/6-31g*", molecule=water, return_wfn=True)
```

`E` is just a number. `wfn` is the whole finished calculation: the molecule, the basis set, the orbitals and their energies, the density, the PSI variables made by the job and, after a frequency job, the Hessian and the vibrational analysis. It keeps all of this until we overwrite it, so it is the safe place to go back to for results. Typing `dir(wfn)` lists its methods: 186 public names for an HF wavefunction in *Psi4* 1.11. **107** of them belong to every wavefunction, whatever method made it, and this page explains all of them, grouped by what they are for. The other 79 are the internal machinery of the SCF procedure; section 12 picks out the few we might want.

(whatyouwilllearnA5)=
## What You Will Learn
- What a wavefunction object holds, and why `return_wfn=True` is worth adding.
- How to get the orbital energies, orbital coefficients and the density from it.
- Where the gradient, the Hessian and the vibrational analysis are stored after a frequency job.
- The traps: results split by symmetry, wavefunctions with no orbitals, and two methods that crash the notebook.

## How to Read This Page

- **Call methods with parentheses.** `wfn.nmo()` gives `19`; `wfn.nmo` without parentheses only names the method.
- **Each calculation gives its own wavefunction.** A new job overwrites the global PSI variables (`psi4.variable(...)`), but not a `wfn` we kept from an earlier job. Give each wavefunction its own name (`wfn_hf`, `wfn_b3lyp`, ...) if you need to compare them.
- **Results are split by symmetry.** Water is $C_{2v}$, with four irreducible representations (irreps), `A1`, `A2`, `B1` and `B2`. Many methods return one block per irrep. Calling `.np` on such an object raises an error that suggests `.nph` instead, which gives a tuple with one NumPy array per irrep. Usually it is easier to ask for the result in the **AO basis**, which comes back as a single ordinary array. The `_subset` methods do this (sections 5 and 6).
- **"Alpha" and "beta"** are the two spin directions. In a closed-shell (RHF or RKS) calculation they are identical, so the `a` and `b` versions of each method (`epsilon_a`, `epsilon_b`) give the same answer. They differ only for open-shell molecules (UHF, ROHF).
- **Everything is in atomic units**: energies in hartree, distances in bohr. Multiply by `psi4.constants.hartree2ev` or `psi4.constants.hartree2kJmol` to convert energies.
- Many methods return *Psi4*'s own types: a `Matrix`, a `Vector` or a `Dimension`. Add `.np` to a `Matrix` or `Vector` to get a NumPy array (if it has only one block). Use `.to_tuple()` on a `Dimension`.
- **A method that has nothing to report returns `None`.** For example, `wfn.gradient()` after an energy job.

```{warning} Two methods crash the notebook
`wfn.atomic_point_charges()` and `wfn.mo_extents()` read results that `psi4.oeprop()` has to make first. Called before that, they do not raise an error: they kill the *Python* kernel, and the notebook has to be restarted. Run `psi4.oeprop(wfn, "MULLIKEN_CHARGES")` or `psi4.oeprop(wfn, "MO_EXTENTS")` before calling them (section 8).
```

All example output below is for HF/6-31G* water at its HF/6-31G* optimized geometry, as defined above.

## 1. Which Wavefunction Did We Get?

Every driver function accepts `return_wfn=True`, but they don't all return the same kind of object.

| Call | What `wfn` is | Has orbitals? |
|---|---|---|
| `psi4.energy("hf/...", return_wfn=True)` | An `RHF` object (`UHF`, `ROHF` or `CUHF` for open shells). DFT gives the same types. | Yes |
| `psi4.optimize(..., return_wfn=True)` | The wavefunction from the last step of the optimization | Yes |
| `psi4.gradient(..., return_wfn=True)` | The wavefunction with the gradient stored on it | Yes |
| `psi4.frequency("hf/...", return_wfn=True)` | HF has an analytic Hessian, so this is an `RHF` with the Hessian added | Yes |
| `psi4.frequency("b3lyp/...", return_wfn=True)` | Made by **finite differences** (many gradient jobs at displaced geometries). This is a bare `Wavefunction` that holds the energy, gradient, Hessian and frequencies of the reference geometry, but **no orbitals**. | **No** |
| Post-HF methods (`mp2`, `ccsd(t)`, ...) | The wavefunction of that method; the HF wavefunction it started from is kept inside (`reference_wavefunction()`) | Varies |

So if a notebook needs both frequencies and orbitals from a DFT method, it needs two calculations: a frequency job for the Hessian, and an energy job for the orbitals. The *Acids and Bases* notebook does exactly this.

`type(wfn).__name__` reports which kind we have: `'RHF'` for our water.

## 2. The Energy and Where It Came From

| Method | What it gives | Water example |
|---|---|---|
| `energy()` | The final energy of the calculation, in hartree | `-76.010721` |
| `name()` | Name of the wavefunction type. For SCF this is the reference (`'RHF'`) even for B3LYP; use `functional().name()` (section 12) for the functional. | `'RHF'` |
| `module()` | Which part of *Psi4* made it | `'scf'` |
| `molecule()` | The **Molecule object** used, at the geometry of the calculation. See the *Molecule Object* appendix. | `wfn.molecule().natom()` → `3` |
| `reference_wavefunction()` | The wavefunction a correlated method (MP2, CCSD, ...) started from, or `None` for an SCF wavefunction | `None` |
| `efzc()` | Energy of the frozen core orbitals, when `freeze_core` is on | `0.0` |
| `options()` | The options object used for the job. Expert use; `psi4.core.get_option()` is easier. | |
| `get_print()` | How much detail the job printed to the output file | `1` |

`wfn.molecule()` is the same object as the molecule we passed in, so after `psi4.optimize()` it already holds the optimized coordinates.

## 3. PSI Variables Stored on the Wavefunction

Each wavefunction keeps a copy of the PSI variables made by its own job. The *PSI Variables* appendix lists the variables themselves; these methods read and write them.

| Method | What it does | Water example |
|---|---|---|
| `variable(name)` | Returns one variable, a number or an array. Not case-sensitive. | `wfn.variable("SCF TOTAL ENERGY")` → `-76.010721` |
| `variables()` | A dictionary of every variable on this wavefunction | `sorted(wfn.variables())` |
| `has_variable(name)` | `True` if the variable has been set | `wfn.has_variable("MULLIKEN CHARGES")` → `False` until `oeprop` has run |
| `set_variable(name, value)` | Stores a number or array under a name of our choice | |
| `del_variable(name)` | Removes a variable | |
| `scalar_variable`, `scalar_variables`, `has_scalar_variable`, `set_scalar_variable`, `del_scalar_variable` | The same, numbers only. Older versions; the plain methods above do both. | |
| `array_variable`, `array_variables`, `has_array_variable`, `set_array_variable`, `del_array_variable` | The same, arrays only | |
| `potential_variable`, `potential_variables`, `has_potential_variable`, `set_potential_variable`, `del_potential_variable` | The same, for external potentials. Provisional; not needed in this book. | |

```{important} Thermochemistry is not stored on the wavefunction
The thermochemistry results of a frequency job (`ZPVE`, `ENTHALPY`, `GIBBS FREE ENERGY` and their corrections) go only into the **global** store. `wfn.variable("ENTHALPY")` raises a `KeyError`. Read them with `psi4.variable("ENTHALPY")` straight after the frequency job, before another calculation overwrites them, and save them in a *Python* variable. The *Thermochemistry* and *Acids and Bases* notebooks do this.
```

## 4. Gradient, Hessian and Frequencies

These are filled in by `psi4.gradient()`, `psi4.hessian()` and `psi4.frequency()`. After a plain energy job they return `None`.

| Method or attribute | What it gives | Water example |
|---|---|---|
| `gradient()` | Gradient of the energy, `natom × 3`, in hartree/bohr. Close to zero at an optimized geometry. | 3 × 3 Matrix |
| `hessian()` | The **Hessian**: second derivatives of the energy, `3natom × 3natom`, in hartree/bohr². | 9 × 9 Matrix |
| `frequencies()` | The vibrational frequencies in cm⁻¹ as a Vector, lowest first; translations and rotations left out. An imaginary frequency appears as a negative number. | `wfn.frequencies().np` → `[1826.0, 4069.8, 4188.5]` |
| `frequency_analysis` | **Not a method** (no parentheses): a dictionary with the full vibrational analysis, described below | `wfn.frequency_analysis["omega"]` |
| `set_gradient(m)`, `set_hessian(m)` | Store a gradient or Hessian on the wavefunction. `vibanal_wfn()` uses the stored Hessian. | |
| `compute_energy()`, `compute_gradient()`, `compute_hessian()` | Run the calculation again on this wavefunction. Used by the driver; we call `psi4.energy()` and friends instead. | |

The Hessian depends only on the potential energy surface. The frequencies depend on the Hessian **and** the masses, and the thermochemistry depends on the frequencies **and** the temperature. That is why `psi4.driver.vibanal_wfn(wfn)` can redo the vibrational analysis for new isotopes (*Isotope Effects*) or a new temperature (*Thermochemistry*) without repeating the expensive Hessian.

### The frequency_analysis dictionary

`wfn.frequency_analysis` (and the dictionary returned by `psi4.driver.vibanal_wfn(wfn)`) holds one entry per **mode**. There are 3N modes, translations and rotations included: 9 for water, the first 6 of which are translations and rotations. Each entry is a `Datum` with `.label`, `.units` and `.data`; the numbers are in `.data`. For the vectors, each **column** is one mode.

| Key | What it is | Units |
|---|---|---|
| `omega` | Frequency, as a **complex** number. Real modes have it in `.real`, imaginary modes in `.imag`. | cm⁻¹ |
| `TRV` | `'TR'` (translation or rotation) or `'V'` (vibration), as a **list** | |
| `gamma` | Symmetry label of each vibration (`'A1'`, `'B2'`, ...); `None` for translations and rotations | |
| `degeneracy` | How many modes share this frequency | |
| `x` | Normal mode as Cartesian displacements, normalized. The one to use for arrows or animations. | bohr |
| `w` | Normal mode as Cartesian displacements, not normalized | bohr |
| `q` | Normal mode in mass-weighted coordinates, normalized | bohr·u½ |
| `mu` | Reduced mass of the mode | u |
| `k` | Force constant of the mode (negative for an imaginary mode) | mdyn/Å |
| `theta_vib` | Characteristic vibrational temperature $h\nu/k_\mathrm{B}$. If it is far above *T*, the mode stays in its ground state. | K |
| `Xtp0`, `Qtp0`, `DQ0` | Classical turning point and RMS spread of the $v = 0$ level. Specialist; not used in this book. | bohr, bohr·u½ |
| `IR_intensity` | IR intensity. Only present when the dipole derivatives were calculated. | km/mol |

The vibrations alone, with their symmetry, reduced mass and force constant:

```python
vib = wfn.frequency_analysis
is_vib = np.array([t == "V" for t in vib["TRV"].data])     # TRV is a list, so build the mask this way
pd.DataFrame({"symmetry": np.array(vib["gamma"].data)[is_vib],
              "freq (cm-1)": vib["omega"].data[is_vib].real,
              "mu (u)": vib["mu"].data[is_vib],
              "k (mdyn/A)": vib["k"].data[is_vib],
              "theta_vib (K)": vib["theta_vib"].data[is_vib]}).round(2)
```

For HF/6-31G* water this gives the bend (A1, 1826 cm⁻¹, *k* = 2.1 mdyn/Å) and the two stretches (A1 at 4070 and B2 at 4189 cm⁻¹, *k* ≈ 10–11 mdyn/Å). All three have $\theta_\mathrm{vib}$ above 2600 K, which is why the vibrations add almost nothing but the ZPVE to the thermal energy at room temperature.

## 5. Electrons and Orbitals: Counting

Methods ending in `pi` ("per irrep") return a `Dimension` with one count for each irrep, in the order of `wfn.molecule().irrep_labels()`. Read it with `.to_tuple()`.

| Method | What it gives | Water example |
|---|---|---|
| `nalpha()`, `nbeta()` | Number of alpha and beta electrons | `5`, `5` |
| `nmo()` | Number of molecular orbitals | `19` |
| `nso()` | Number of basis functions (symmetry orbitals). Equal to `nmo()` unless some were dropped as near-duplicates. | `19` |
| `nirrep()` | Number of irreps in the computational point group | `4` |
| `nmopi()`, `nsopi()` | Orbitals per irrep | `(10, 1, 3, 5)` |
| `doccpi()` | Doubly occupied orbitals per irrep | `(3, 0, 1, 1)` |
| `soccpi()` | Singly occupied orbitals per irrep (open shells only) | `(0, 0, 0, 0)` |
| `nalphapi()`, `nbetapi()` | Alpha and beta electrons per irrep | `(3, 0, 1, 1)` |
| `nfrzc()` | Number of frozen core orbitals | `0` |
| `frzcpi()`, `frzvpi()` | Frozen core and frozen virtual orbitals per irrep | `(0, 0, 0, 0)` |
| `same_a_b_orbs()`, `same_a_b_dens()` | `True` if alpha and beta orbitals (or densities) are identical, as in RHF | `True` |
| `force_occpi(docc, socc)` | Forces an occupation pattern. Expert use; the `docc` and `socc` options are easier. | |

`doccpi()` read with the irrep labels says that water's five occupied orbitals are three of A1 symmetry, one B1 and one B2, matching the qualitative MO diagram ($1a_1^2\,2a_1^2\,1b_2^2\,3a_1^2\,1b_1^2$).

## 6. Orbital Energies and Orbital Coefficients

The `_subset` versions are the ones to use. They take two arguments:

- **basis**: `"AO"` gives one ordinary array covering all symmetries, sorted by energy. `"SO"` (or `"MO"` for the density) gives one block per irrep.
- **subset**: `"ALL"`, `"OCC"` (occupied), `"VIR"` (virtual, i.e. empty), or the frozen and active versions `"FROZEN_OCC"`, `"ACTIVE_OCC"`, `"ACTIVE_VIR"`, `"FROZEN_VIR"`.

| Method | What it gives | Water example |
|---|---|---|
| `epsilon_a_subset(basis, subset)`, `epsilon_b_subset(...)` | Orbital energies, in hartree, as a Vector | `wfn.epsilon_a_subset("AO", "OCC").np` → `[-20.558, -1.346, -0.714, -0.571, -0.498]` |
| `epsilon_a()`, `epsilon_b()` | All orbital energies, one block per irrep (use `.nph`) | |
| `Ca_subset(basis, subset)`, `Cb_subset(...)` | **Orbital coefficients**: one column per MO, one row per basis function | `wfn.Ca_subset("AO", "OCC").np.shape` → `(19, 5)` |
| `Ca()`, `Cb()` | All coefficients, one block per irrep | |
| `alpha_orbital_space(...)`, `beta_orbital_space(...)` | Orbital sets for correlated-method programs. Expert use. | |

The HOMO and LUMO energies in eV (the *Acids and Bases* and *Molecular Orbitals* notebooks use these):

```python
h2ev = psi4.constants.hartree2ev
homo = wfn.epsilon_a_subset("AO", "OCC").np[-1] * h2ev    # -13.56 eV
lumo = wfn.epsilon_a_subset("AO", "VIR").np[0] * h2ev     #  +5.80 eV
```

By Koopmans' theorem, $-\varepsilon_\mathrm{HOMO}$ estimates the first ionization energy: 13.6 eV here, against 12.6 eV measured for water.

To label each orbital with its symmetry, combine the per-irrep energies with the irrep labels:

```python
labels = wfn.molecule().irrep_labels()
eps = wfn.epsilon_a().nph                                   # one array per irrep
orbitals = sorted((e, labels[h]) for h in range(wfn.nirrep()) for e in eps[h])
# [(-20.558, 'A1'), (-1.346, 'A1'), (-0.714, 'B2'), (-0.571, 'A1'), (-0.498, 'B1'), (0.213, 'A1'), ...]
```

The HOMO of water is the B1 orbital: the oxygen *p* orbital perpendicular to the molecular plane.

## 7. Matrices: Density, Fock, Overlap and Core Hamiltonian

These are the matrices the SCF procedure is built from. We rarely need them in this book, but they are there for anyone who wants to see the machinery. All are square, size `nso × nso` (19 × 19 for water in 6-31G*).

| Method | What it gives |
|---|---|
| `Da_subset("AO")`, `Db_subset(...)` | **Density matrix** in the AO basis (one array). `"SO"` or `"MO"` give it per irrep. |
| `Da()`, `Db()` | Density matrix, one block per irrep |
| `Fa_subset("AO")`, `Fb_subset(...)` | **Fock matrix** in the AO basis. The orbitals and orbital energies are its eigenvectors and eigenvalues. |
| `Fa()`, `Fb()` | Fock matrix, one block per irrep |
| `S()` | **Overlap matrix** of the basis functions, one block per irrep |
| `H()` | **Core Hamiltonian**: kinetic energy plus attraction to the nuclei, one block per irrep |
| `lagrangian()` | Orbital Lagrangian, used for gradients of correlated methods. Expert use. |
| `get_density(name)` | A density by name. Marked experimental in *Psi4*. |
| `aotoso()` | The matrix that converts AO quantities into symmetry-blocked SO quantities |
| `no_occupations()` | Natural orbital occupations, for correlated methods that make them |

A quick check that the density holds the right number of electrons: $\mathrm{tr}(\mathbf{D}\mathbf{S})$ counts the alpha electrons.

```python
S_ao = psi4.core.MintsHelper(wfn.basisset()).ao_overlap().np
print(np.trace(wfn.Da_subset("AO").np @ S_ao))    # 5.0000
```

## 8. Properties after psi4.oeprop()

`psi4.oeprop(wfn, ...)` calculates one-electron properties from a finished wavefunction and stores them on it. These methods read the results. Remember the warning above: the first two crash the kernel if `oeprop` has not made their results.

| Method | What it gives | Made by | Water example |
|---|---|---|---|
| `atomic_point_charges()` | Atomic charges from the **most recent** charge calculation (Mulliken, Löwdin or MBIS), as a Vector | `oeprop(wfn, "MULLIKEN_CHARGES")` | `[-0.869, 0.434, 0.434]` |
| `mo_extents()` | Spatial extent ⟨x²⟩, ⟨y²⟩, ⟨z²⟩, ⟨r²⟩ of each MO, as a list of four Vectors | `oeprop(wfn, "MO_EXTENTS")` | |
| `esp_at_nuclei()` | Electrostatic potential at each nucleus, in hartree/e | `oeprop(wfn, "ESP_AT_NUCLEI")` | `[-22.314, -0.997, -0.997]` |
| `get_dipole_field_strength()` | Strength of an applied electric field, if one was set with the `perturb_h` options | | `[0.0, 0.0, 0.0]` |

It is safer to read charges by name, because `atomic_point_charges()` only holds the last kind calculated:

```python
psi4.oeprop(wfn, "MULLIKEN_CHARGES", "LOWDIN_CHARGES")
wfn.variable("MULLIKEN CHARGES")      # [-0.869, 0.434, 0.434]
```

The dipole moment needs no `oeprop` call; every SCF job stores it: `wfn.variable("SCF DIPOLE")` gives the vector in $e\,a_0$, here `[0, 0, 0.865]`, which is 2.20 D.

## 9. The Basis Set

| Method | What it gives | Water example |
|---|---|---|
| `basisset()` | The orbital basis set as a `BasisSet` object | `wfn.basisset().name()` → `'6-31G*'` |
| `get_basisset(key)` | An auxiliary basis, e.g. the density-fitting basis | `wfn.get_basisset("DF_BASIS_SCF").name()` |
| `set_basisset(key, basis)` | Stores an auxiliary basis on the wavefunction. Expert use. | |
| `sobasisset()` | The basis rearranged into symmetry orbitals | |
| `basis_projection(C, dim, old_basis, new_basis)` | Projects orbitals from one basis set into another. *Psi4* does this itself with the `basis_guess` option. | |
| `mintshelper()` | A `MintsHelper` object that calculates integrals (overlap, kinetic, dipole, ...) in this basis. Expert use. | |

The `BasisSet` object has its own methods. The useful ones are `nbf()` (number of basis functions, `19`), `nshell()` (number of shells, `10`) and `name()`.

## 10. Saving, Loading and Copying

| Method | What it does |
|---|---|
| `write_molden(filename)` | Writes the orbitals to a **Molden** file for viewers such as Avogadro, Jmol or IboView. The same as `psi4.molden(wfn, filename)`. |
| `write_nbo(filename)` | Writes the wavefunction in the input format of the NBO program |
| `to_file(filename)` | Saves the whole wavefunction to a NumPy `.npy` file. Lets a later notebook reuse a long calculation. |
| `psi4.core.Wavefunction.from_file(filename)` | Loads a wavefunction saved with `to_file()`. Called on the class, like `Molecule.from_string`. The result is a plain `Wavefunction`. |
| `build(molecule, basis)` | Builds an **empty** wavefunction for a molecule and basis, with no calculation done. Called on the class: `psi4.core.Wavefunction.build(water, "cc-pvdz")`. Used for guesses and tests. |
| `deep_copy(other)` | Copies all the data of another wavefunction into this one |
| `shallow_copy(other)` | Makes this wavefunction share the data of another one, so changing one changes both. Expert use. |
| `c1_deep_copy(basis)` | A copy with the symmetry switched off ($C_1$). Some add-on programs need this. |

Saving a wavefunction and getting it back:

```python
wfn.to_file(str(output_dir / "water_wfn"))                             # writes water_wfn.npy
wfn2 = psi4.core.Wavefunction.from_file(output_dir / "water_wfn.npy")
wfn2.energy()                                                          # -76.010721
```

`to_file()` needs the file name as text, hence the `str(...)`; given a `Path` it stops with an `AttributeError`. `from_file()` accepts either.

## 11. Changing the Wavefunction and Specialist Tools

We never need these in the book. They are for writing new methods or connecting *Psi4* to other programs.

| Method | What it does |
|---|---|
| `set_energy(E)` | Changes the stored energy (and the `CURRENT ENERGY` variable) |
| `set_name(text)`, `set_module(text)` | Change the name or module label |
| `set_lagrangian(m)` | Stores an orbital Lagrangian |
| `set_reference_wavefunction(wfn)` | Stores the reference wavefunction |
| `set_print(level)` | Changes how much the wavefunction prints |
| `get_scratch_filename(n)` | Path of a scratch file *Psi4* writes during the job |
| `PCM_enabled()` | `True` if the job used the PCM solvent model (*Solvation*) |
| `get_PCM()`, `set_PCM(pcm)` | The PCM solvent object |
| `external_pot()`, `set_external_potential(pot)` | Point charges or fields placed around the molecule (QM/MM-style calculations) |

## 12. Extras on SCF Wavefunctions

An HF or DFT wavefunction (`RHF`, `UHF`, `ROHF`, ...) has 79 more methods. Almost all are steps of the SCF procedure itself (`form_F`, `form_C`, `diis`, `damping_update`, `iteration_energies`, ...) that the program calls while it converges. A few are worth knowing:

| Method or attribute | What it gives | Water example |
|---|---|---|
| `functional()` | The density functional object. `.name()` gives `'HF'` for Hartree–Fock, `'B3LYP'` for B3LYP. | `wfn.functional().name()` → `'HF'` |
| `iteration_` | Number of SCF iterations needed (an attribute: no parentheses) | `8` |
| `get_energies(term)` | One part of the energy: `"Nuclear"`, `"One-Electron"`, `"Two-Electron"`, `"XC"` (DFT only) or `"Total Energy"` | `wfn.get_energies("Nuclear")` → `9.288` |
| `occupation_a()`, `occupation_b()` | Occupation number of each orbital, one block per irrep (`.nph`): 1 or 0 for HF, fractions for fractional-occupation calculations | |
| `compute_spin_contamination()` | Prints ⟨*S*²⟩ for an open-shell (UHF) wavefunction to the output file | |
| `stability_analysis()` | Checks whether a lower-energy solution exists. The `stability_analysis` option is the usual way (*Bond Dissociation*). | |
| `print_orbitals()`, `print_energies()` | Write the orbital energies or energy breakdown to the output file | |

For water, the energy parts add up to the total: 9.288 (nuclear repulsion) − 123.224 (one-electron) + 37.925 (two-electron) = −76.011 hartree.

## Summary
- Add `return_wfn=True` to keep the whole calculation in a Wavefunction object; it is not overwritten by later jobs.
- Get orbital energies and coefficients with `epsilon_a_subset("AO", ...)` and `Ca_subset("AO", ...)`; the `"AO"` versions avoid the symmetry blocks.
- After a frequency job, `hessian()`, `frequencies()` and the `frequency_analysis` dictionary hold the vibrational results, but the thermochemistry (`ENTHALPY`, `GIBBS FREE ENERGY`) is only in the global PSI variables.
- A DFT frequency job made by finite differences returns a wavefunction **with no orbitals**; run an energy job if you need them.
- Run `psi4.oeprop()` before `atomic_point_charges()` or `mo_extents()`, or the notebook crashes.

```{warning} AI Statement
:class: dropdown
This page was written by Claude (Anthropic) using Claude Code, at my request. Claude listed the methods with `dir()`, read their descriptions in *Psi4* 1.11 and the source of the vibrational analysis (`psi4/driver/qcdb/vib.py`), and ran each example on HF/6-31G* water to check the output shown. I selected, edited and checked all content.
```
