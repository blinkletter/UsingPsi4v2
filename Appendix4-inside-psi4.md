# What Is Inside psi4

Every notebook in this book starts with `import psi4`. That one line brings in the whole program. Typing `dir(psi4)` lists everything it contains: **180 names** in *Psi4* 1.11. A handful of them do almost all the work in this book. The rest are less common calculations, settings, error types, tools for the program's own tests, and the internal machinery that makes it run.

This page goes through all 180 names, grouped by what they are for. The first section covers the ones we use all the time; later sections can be skimmed.

(whatyouwilllearnA4)=
## What You Will Learn
- Which functions in `psi4` run calculations, and which ones set up the job.
- How to check what your installation can do: version, add-ons, available methods.
- How to catch a calculation that fails without stopping the whole notebook.
- What the many submodules are, so an unfamiliar name in an error message makes sense.

## How to Read This Page

- `dir(psi4)` lists **functions** (things to call, like `psi4.energy(...)`), **modules** (collections of more functions, like `psi4.core`), **classes** (blueprints for objects, like error types) and plain **values** (like `psi4.__version__`).
- Names are case-sensitive: `psi4.energy`, not `psi4.Energy`.
- To read the built-in help for any name, use `help(psi4.energy)`, or type `psi4.energy?` in a notebook cell.
- The examples were run in the `psi4p314` environment. Your numbers (memory, threads, add-ons) depend on your own installation.

## 1. The Everyday Functions

These are the names used throughout the book. Most notebooks need no others.

| Name | What it does | Used in |
|---|---|---|
| `geometry(text, name=...)` | Builds a **Molecule** from a geometry string. See the *Molecule Object* appendix. | Every notebook |
| `set_options({...})` | Sets calculation options: `basis`, `reference`, `scf_type`, `freeze_core`, ... | Every notebook |
| `set_memory("2 GB")` | Sets the memory *Psi4* may use. | Every notebook |
| `set_num_threads(4)` | Sets how many CPU cores *Psi4* may use. | Every notebook |
| `set_output_file("file.out", False)` | Sends the long text output to a file; `True` appends instead of overwriting. | Every notebook |
| `energy("hf/6-31g*", molecule=mol)` | **Single-point energy** at a fixed geometry. Add `return_wfn=True` to also get the wavefunction. | *Water*, *Bond Dissociation*, *Surfaces*, ... |
| `optimize(...)`, `opt(...)` | **Geometry optimization**: moves the atoms to the lowest-energy structure and returns its energy. `opt` is a short alias. | *Optimizing*, *Level of Effort*, ... |
| `frequency(...)`, `frequencies(...)`, `freq(...)` | **Vibrational frequencies** and thermochemistry. Three names for the same function. | *Vibrations*, *Thermochemistry*, ... |
| `variable(name)` | Reads a stored result (a **PSI variable**). See the *PSI Variables* appendix. | *Thermochemistry*, *Hydrogen Bonding*, ... |
| `set_variable(name, value)` | Stores your own value as a PSI variable. Rarely needed. | |
| `oeprop(wfn, "DIPOLE", "MULLIKEN_CHARGES", ...)` | One-electron **properties** (dipoles, atomic charges) from a finished wavefunction. | *Charges and Dipole*, *Solvation* |
| `cubeprop(wfn)` | Writes **cube files** of orbitals, densities or the electrostatic potential for drawing surfaces. | *Surfaces*, *Molecular Orbitals*, *Hydrogen Bonding* |
| `molden(wfn, "file.molden")` | Saves orbitals in Molden format for other viewers. | *Molecular Orbitals* |
| `pcm_helper(text)` | Passes the PCM solvent description to *Psi4*. | *Solvation* |
| `constants` | Physical constants and unit conversions: `psi4.constants.hartree2kJmol` (2625.50), `psi4.constants.bohr2angstroms` (0.529177), ... | Most notebooks |
| `core` | The compiled C++ heart of *Psi4*. Holds the `Molecule`, `Wavefunction` and `Matrix` types and tools such as `psi4.core.clean()`, `psi4.core.clean_options()` and `psi4.core.variables()`. | Most notebooks |
| `driver` | The *Python* layer that turns `energy()` and friends into the right sequence of C++ calls. We use `psi4.driver.vibanal_wfn` from it in *Thermochemistry*. | *Thermochemistry* |

## 2. Other Ways to Run a Calculation

| Name | What it does |
|---|---|
| `gradient(...)` | One **gradient** (force on each atom) at a fixed geometry, returned as an `natom × 3` Matrix. `optimize()` calls this at every step. |
| `hessian(...)` | The **Hessian** (second derivatives) at a fixed geometry, without the frequency analysis. |
| `vibanal_wfn(wfn)` | Repeats the frequency and thermochemistry analysis from a Hessian already stored in `wfn`, for example at a new temperature. |
| `properties(...)`, `prop(...)` | Runs a calculation and evaluates properties (dipole, polarizability, ...) in one call; returns the energy. |
| `cbs(...)` | **Complete-basis-set** schemes: extrapolates energies from a series of basis sets and adds corrections. The short form `psi4.energy("mp2/cc-pv[dt]z")` does the same thing more easily. |
| `tdscf(wfn, states=5)` | **Excited states** (TDDFT or TDHF) from a converged SCF wavefunction: excitation energies for a predicted UV–vis spectrum. |
| `optimize_geometric(...)` | Optimization with the external *geomeTRIC* optimizer. Normally reached with `psi4.optimize(..., engine="geometric")`. |
| `fsapt_analysis(...)` | **Functional-group SAPT**: splits a SAPT interaction energy among chemical groups. |
| `gdma(wfn)` | Distributed multipole analysis with Stone's GDMA program (an add-on; not installed in `psi4p314`). |
| `fchk(wfn, "file.fchk")` | Saves the wavefunction as a Gaussian formatted checkpoint file, readable by many viewers. |
| `basis_helper(text)` | Defines a custom or mixed basis set, for example a larger basis on one atom. |
| `activate(mol)` | Makes `mol` the default molecule for later calls. Passing `molecule=mol` is clearer. |
| `scf_helper(...)`, `scf_wavefunction_factory(...)` | Internal helpers that set up and run the SCF step. Used by the driver, not by us. |
| `set_module_options(...)` | Old way of setting options for one module. Deprecated; use `set_options` with `"module__option"` names. |
| `run_qcschema(...)` | Runs a job described as a QCSchema dictionary, the format used to connect quantum chemistry programs. |
| `process_input(text)` | Converts an old-style *Psi4* input file ("Psithon") into *Python*. Used when running `psi4 input.dat` from a terminal. |
| `mdi_run(...)`, `ipi_broker(...)` | Let another program (a molecular dynamics code, for example) drive *Psi4* through the MDI or i-PI interfaces. |

Check whether a gradient method exists before starting a long optimization:

```python
"mp2" in psi4.procedures["gradient"]       # True: analytic MP2 gradients exist
```

## 3. Your Installation and Its Settings

| Name | What it gives | Example |
|---|---|---|
| `__version__` | *Psi4* version | `'1.11'` |
| `version_formatter(fmt)` | The version in a chosen format | |
| `get_memory()` | Memory allowance in bytes | `2000000000` after `set_memory("2 GB")` |
| `get_num_threads()` | Number of threads | `4` after `set_num_threads(4)` |
| `addons(name)` | `True` if an optional add-on is installed; with no argument, a list (unavailable ones struck through) | `psi4.addons("pcmsolver")` → `True` |
| `executable` | Path to the `psi4` command-line program | |
| `envvar_scratch` | Scratch folder for large temporary files (from the `PSI_SCRATCH` environment variable) | |
| `get_input_directory()` | Folder of the input file when run from a terminal | |
| `print_header()` | Prints the *Psi4* banner and citation to the output file | |
| `citation_formatter(fmt)` | The paper to cite for *Psi4* | `psi4.citation_formatter("{doi}")` → `'10.1063/5.0006002'` |
| `copy_file_to_scratch(...)`, `copy_file_from_scratch(...)` | Move files in and out of the scratch folder by *Psi4*'s naming rules | |
| `libint2_configuration()`, `libint2_print_out()` | Details of the Libint2 integral library that *Psi4* was built with | |
| `test("smoke")` | Runs *Psi4*'s own test suite: `"smoke"` is a quick check that the installation works | |

In `psi4p314`, the installed add-ons include `pcmsolver` (PCM solvation), `ddx` (ddCOSMO/ddPCM solvation), `dftd3` (D3 dispersion) and `geometric` (alternative optimizer). `gdma`, `dftd4`, `gcp` and `cfour` are not installed.

## 4. Errors: When a Calculation Fails

When something goes wrong, *Psi4* **raises** one of these error types. Catching the right one with `try`/`except` lets a scan or loop carry on past a failed point instead of stopping the notebook.

| Name | When it is raised |
|---|---|
| `PsiException` | Parent of all *Psi4* errors; catching it catches them all. |
| `ConvergenceError` | Parent of the three convergence errors below. |
| `SCFConvergenceError` | The SCF did not converge within `maxiter` cycles. Common for radicals and stretched bonds. Holds the unconverged wavefunction as `.wfn`. |
| `OptimizationConvergenceError` | The geometry optimization ran out of steps. Holds the last wavefunction as `.wfn`, so the last geometry is not lost. |
| `TDSCFConvergenceError` | The excited-state calculation did not converge. |
| `MissingMethodError` | The method name is unknown, often a typing error (`"b3lpy"`). |
| `ManagedMethodError` | The method exists but not with these settings, for example a gradient for a method/`*_type` combination that has none. |
| `ValidationError` | An input value is wrong: a bad option value, an impossible charge or multiplicity. |
| `ParsingError` | A text input (geometry or option string) can't be read. |
| `UpgradeHelper` | Syntax that worked in an older *Psi4* version; the message says what to write now. |
| `PastureRequiredError` | The option needs a retired module from the separate *Psi4Pasture* collection. |
| `TestComparisonError` | A `compare_*` check (section 6) failed. |

The *Bond Dissociation* notebook uses this pattern:

```python
try:
    E = psi4.energy("scf", molecule=mol)
except psi4.SCFConvergenceError:
    E = np.nan            # record a gap and carry on with the next distance
```

## 5. Tables of What Psi4 Can Do

| Name | What it holds |
|---|---|
| `procedures` | A dictionary of every method name *Psi4* accepts, sorted by job type: `"energy"`, `"gradient"`, `"hessian"`, `"properties"`. In 1.11 it lists about 3100 energy methods (most are DFT functionals with and without dispersion). |
| `energy_only_methods` | Methods with energies but no gradients: SAPT variants, `cphf`, `tdhf`, `cis`, `efp`. Optimizations and frequencies are not possible with them. |
| `integrated_basis_methods` | Methods with their own built-in basis set, called without `/basis`: `hf-3c`, `pbeh-3c`, `b97-3c`, `r2scan-3c`, `wb97x-3c`, `gaussian-2`, `sns-mp2`. |
| `hooks` | Lists of your own functions to run before or after each `energy`, `optimize` or `frequency` call. Expert use. |

## 6. Checking and Printing Results

The `compare_*` functions are how *Psi4*'s test suite checks that answers are right. Each prints `PASSED` and returns `True`, or raises `TestComparisonError`. They can be used in teaching notebooks to check a student's result.

| Name | Compares |
|---|---|
| `compare_values(expected, computed, digits, label)` | Numbers or arrays to a given number of decimal places |
| `compare_integers(...)`, `compare_strings(...)`, `compare(...)` | Exact matches: integers, strings, booleans |
| `compare_arrays(...)`, `compare_matrices(...)`, `compare_vectors(...)` | NumPy arrays, `Matrix` and `Vector` objects |
| `compare_recursive(...)` | Nested dictionaries and lists |
| `compare_molrecs(...)` | Molecules converted with `to_dict()` |
| `compare_wavefunctions(...)` | Two wavefunctions |
| `compare_cubes(...)`, `compare_fchkfiles(...)`, `compare_moldenfiles(...)`, `compare_fcidumps(...)` | Output files of each type |

```python
psi4.compare_values(-76.0105, E, 4, "HF/6-31G* water")   # prints PASSED, returns True
```

Small text helpers:

| Name | What it does |
|---|---|
| `nppp(array)`, `nppp10(array)` | Format a NumPy array as tidy text (8 or 10 decimals, tiny values shown as 0); use with `print()` |
| `banner(text)` | Prints a heading into the output file |
| `message_box(text)` | Returns `text` inside a box of `!` characters |
| `levenshtein(a, b)` | How many letter changes turn one string into another: `levenshtein("b3lpy", "b3lyp")` → `2` |
| `find_approximate_string_matches(word, options, n)` | The "did you mean...?" behind error messages: `"b3lpy"` → `['b3lyp', 'b3lyps', 'b3lyp5']` |
| `sanitize_method(name)`, `docs_table_link(name, mode)` | Turn a method name into a link to its manual entry, for error messages |

## 7. Objects for Special Jobs

| Name | What it is |
|---|---|
| `EmpiricalDispersion` | Calculates a -D2/-D3/-D4 dispersion correction on its own, without a DFT calculation. |
| `QMMMbohr` | A set of external point charges (positions in bohr) that surround the molecule, for simple QM/MM models. |
| `QMMM` | Older version of `QMMMbohr`; use `QMMMbohr` instead. |
| `AtomicComputer` | One single-geometry task inside a larger job (for example one displacement of a finite-difference frequency). Internal to the driver. |

## 8. Integral Files (FCIDUMP)

FCIDUMP is a standard file format for passing integrals to other correlated-method programs.

| Name | What it does |
|---|---|
| `fcidump(wfn)` | Writes the integrals of a wavefunction to an FCIDUMP file |
| `fcidump_from_file(filename)` | Reads one back in |
| `energies_from_fcidump(data)` | Recomputes SCF and MP2 energies from those integrals as a check |

## 9. Plugins and Molecule Internals

| Name | What it does |
|---|---|
| `create_plugin(name, template)` | Makes a starting folder for writing your own *Psi4* add-on in C++ or *Python* |
| `sanitize_name(name)` | Tidies a plugin name |
| `pluginutil`, `endorsed_plugins` | Modules that create and load plugins |
| `dynamic_variable_bind(cls)` | Adds the *Python*-side methods (such as `run_dftd3`, `to_string`) to the C++ `Molecule` class when `psi4` is imported |
| `molecule_get_attr`, `molecule_set_attr` | Let you read and set Z-matrix variables as attributes, e.g. `water.r = 1.0` for a variable called `r` |

## 10. Psi4's Own Submodules

These modules hold the code behind everything above. We rarely call them directly, but their names appear in error messages and tracebacks.

**Shared data and utilities**

| Module | What it contains |
|---|---|
| `p4util` | Utility functions for the driver and users; many top-level names above come from here |
| `qcdb` | Older molecule, basis set and database tools |
| `qcel`, `qcelemental` | The *QCElemental* package (two names for one module): periodic table data, physical constants, unit conversions, molecule parsing. `psi4.constants` comes from here. |
| `molutil` | Functions that build and change molecules, including `geometry()` |
| `psif`, `psifiles` | Numbers of *Psi4*'s scratch files (two names for one module) |
| `header`, `metadata` | The banner, citation and version information |
| `extras` | Output-file setup, add-on checks and the test runner |
| `inputparser` | Converts Psithon input files into *Python* |
| `util` | Terminal colour helpers |

**How a job is run**

| Module | What it contains |
|---|---|
| `driver` | The main user functions: `energy`, `optimize`, `frequency`, ... |
| `driver_cbs`, `driver_cbs_helper` | Complete-basis-set extrapolation |
| `driver_findif` | Finite differences: gradients from energies, Hessians from gradients, when no analytic derivative exists |
| `driver_nbody` | Many-body and counterpoise (`bsse_type`) jobs |
| `driver_util` | Checks which analytic derivatives a method has |
| `task_base`, `task_planner` | Break a job into separate calculations (`AtomicComputer` tasks) and plan them |
| `schema_wrapper` | Runs QCSchema jobs (`run_qcschema`) |
| `procrouting` | Routes each method name to the code that runs it; contains the next eight |
| `proc`, `proc_util` | The functions that run each method (`run_scf`, `run_dfmp2`, ...) and their helpers |
| `proc_table`, `proc_data` | The `procedures` table and data behind it |
| `scf_proc` | The SCF iterations |
| `diis` | DIIS, the method that speeds up SCF convergence |
| `dft` | Definitions of the 1395 built-in DFT functionals: `psi4.dft.functionals` |
| `empirical_dispersion` | The -D dispersion corrections |

**Particular kinds of calculation**

| Module | What it contains |
|---|---|
| `sapt`, `fsapt` | SAPT and functional-group SAPT |
| `solvent` | The `ddx` (ddCOSMO/ddPCM) and `efp` (effective fragment potential) solvation models |
| `qmmm` | External point charges (`QMMMbohr`) |
| `response` | SCF response: polarizabilities and excited states (`tdscf`) |
| `roa` | Raman optical activity |
| `findif_response_utils` | Finite-difference helpers for response properties |
| `mcscf` | Multiconfigurational SCF (CASSCF), for molecules a single HF configuration describes poorly |
| `libcubeprop` | The code behind `cubeprop()` |
| `diatomic` | `anharmonicity()`: fits a potential energy curve of a diatomic and returns spectroscopic constants ($\omega_e$, $\omega_e x_e$, $r_e$, ...). It suits the *Bond Dissociation* notebook. |
| `frac` | Fractional occupation of orbitals, used to study ionization and tune functionals |
| `gaussian_n` | The G2 composite method (`psi4.energy("gaussian-2")`) |
| `aliases` | Ready-made composite recipes such as `sherrill_gold_standard` and `allen_focal_point` |
| `wrapper_database` | Runs a method over a whole benchmark set of molecules (`database()`) |
| `wrapper_autofrag` | Splits a molecule into fragments by bonding (`auto_fragments()`) |
| `interface_cfour` | Runs the separate CFOUR program from *Psi4* |
| `mdi_engine` | The MDI interface behind `mdi_run()` |

## 11. Standard Python Names

`psi4` imports some ordinary *Python* modules for its own use, and they show up in `dir(psi4)` too. Import them yourself (`import numpy as np`) rather than borrowing them from `psi4`.

| Names | What they are |
|---|---|
| `np` | NumPy |
| `os`, `sys`, `pathlib`, `Path`, `shutil` | Files, folders and the operating system |
| `json`, `re`, `copy` | JSON files, text patterns (regular expressions), copying objects |
| `logging`, `logger`, `warnings`, `atexit` | Logging messages, warnings and clean-up when *Python* exits |
| `pp` | A "pretty printer" for nested lists and dictionaries |
| `Dict`, `Optional`, `Tuple`, `Union` | Type hints used in the code |
| `__version__`, `__name__`, `__file__`, `__path__`, `__package__`, `__doc__`, `__spec__`, `__loader__`, `__cached__`, `__builtins__` | Module information that every *Python* package has. `__file__` shows where *Psi4* is installed. |

## Summary
- `dir(psi4)` lists 180 names, but the book needs only about 15: `geometry`, `set_options`, `set_memory`, `set_num_threads`, `set_output_file`, `energy`, `optimize`, `frequency`, `variable`, `oeprop`, `cubeprop`, `molden`, `pcm_helper`, `constants` and `core`.
- `__version__`, `addons()` and `procedures` tell you what your installation can do.
- Catch `psi4.SCFConvergenceError` and the other error types to keep loops and scans running.
- Most of the remaining names are submodules that the main functions use behind the scenes.

```{warning} AI Statement
:class: dropdown
This page was written by Claude (Anthropic) using Claude Code, at my request. Claude listed `dir(psi4)` in *Psi4* 1.11, read the docstrings and source code for each name, and ran the examples shown. I selected, edited and checked all content.
```
