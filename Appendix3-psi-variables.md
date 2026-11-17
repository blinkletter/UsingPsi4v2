
# PSI Variables

A *Psi4* calculation produces far more numbers than the single energy it returns. The rest are kept as **PSI variables**: results stored under fixed names written in capitals, such as `"SCF TOTAL ENERGY"`, `"ZPVE"` or `"GIBBS FREE ENERGY"`. Several notebooks in this book read them, for example the *Thermochemistry*, *Charges and Dipole Moments* and *Hydrogen Bonding* notebooks.

This page explains the variables that the kinds of job used in this book create. To make the list, we ran twelve typical jobs on water (and on a water dimer for the last two sections), cleared the store before each one, and recorded what `psi4.core.variables()` held afterwards. The values in the tables come from those runs.

(whatyouwilllearnA3)=
## What You Will Learn
- How to read, list and clear PSI variables.
- How the names are built, so an unfamiliar name can be decoded.
- Which variables each type of job creates, what they mean, and which ones are worth using.

## Reading PSI Variables

```python
E, wfn = psi4.energy("hf/6-31g*", molecule=water, return_wfn=True)

psi4.variable("SCF TOTAL ENERGY")        # one variable, from the global store
wfn.variable("SCF TOTAL ENERGY")         # the same, from the wavefunction
psi4.core.variables()                    # every variable: {name: value}
psi4.core.scalar_variables()             # only single numbers
psi4.core.array_variables()              # only arrays (dipoles, gradients, Hessians)
psi4.core.clean_variables()              # empty the global store
```

Things to know:

- **There are two stores.** The **global** store (`psi4.variable`, `psi4.core.variables`) holds the results of the most recent calculations and is overwritten by each new job. The **wavefunction** (`wfn.variable`, `wfn.variables()`) keeps its own copy that stays with that calculation. Use `return_wfn=True` whenever you need results from more than one job.
- **Some variables live only on the wavefunction.** Atomic charges from `psi4.oeprop()` (`MULLIKEN CHARGES`, `LOWDIN CHARGES`, `MBIS CHARGES`, ...) and the detailed SAPT terms (`E ELST10`, `E DISP20`, ...) do **not** appear in `psi4.core.variables()`. Read them with `wfn.variable()`.
- **Energies are in hartree** ($E_\mathrm{h}$). Multiply by `psi4.constants.hartree2kJmol` (2625.5) for kJ/mol. Dipoles are in $e\,a_0$; multiply by `psi4.constants.dipmom_au2debye` (2.5417) for debye.
- **Arrays** come back as *Psi4* Matrix objects. Use `np.asarray(...)` or `.np` to turn them into NumPy arrays.
- **The full catalogue** of every name *Psi4* can create is the manual's [PSI Variables glossary](https://psicode.org/psi4manual/master/glossary_psivariables.html). Many hundreds exist; this page covers the ones our jobs actually make.

## How the Names Are Built

Most names follow a few patterns. Once you know them, a new name is easy to read.

| Pattern | Meaning | Example |
|---|---|---|
| `<METHOD> TOTAL ENERGY` | Total electronic energy at that level | `MP2 TOTAL ENERGY` |
| `<METHOD> CORRELATION ENERGY` | Correlation part only: total minus the SCF energy | `CCSD CORRELATION ENERGY` |
| `<METHOD> TOTAL GRADIENT` / `TOTAL HESSIAN` | First / second derivatives of that energy | `HF TOTAL GRADIENT` |
| `<METHOD> DIPOLE` | Dipole moment at that level | `SCF DIPOLE` |
| `CURRENT ...` | Whatever the **last** job produced, whatever the method | `CURRENT ENERGY` |
| `SCF ...` and `HF ...` | `SCF` covers both HF and DFT; `HF` is set only for Hartree–Fock | `SCF TOTAL ENERGY` |
| `... ITERATIONS` | How many cycles an iterative step needed | `SCF ITERATIONS` |
| `<title> DIPOLE` | Properties from `psi4.oeprop(..., title="GAS")` are named after the `title` | `GAS DIPOLE`, `PCM DIPOLE` |

The `CURRENT` variables are the convenient ones to use in a loop that changes method, because the name stays the same.

## 1. Every SCF Job (HF and DFT)

**Created by:** any `energy`, `optimize` or `frequency` job, because every method starts with an SCF step. RHF and UHF jobs create exactly the same names. Values are for HF/6-31G* water.

| Variable | Water value | What it means |
|---|---|---|
| `CURRENT ENERGY` | −76.0105 | Energy of the most recent job, at the highest level it ran. The same number the function returns. |
| `CURRENT REFERENCE ENERGY` | −76.0105 | The SCF energy underneath the current result. For MP2 or CCSD(T) this is the HF energy the correlation was added to. |
| `SCF TOTAL ENERGY` | −76.0105 | Final SCF energy (HF or DFT), including any solvation or dispersion terms. |
| `HF TOTAL ENERGY` | −76.0105 | Final Hartree–Fock energy (not set for DFT). |
| `SCF ITERATION ENERGY` | −76.0105 | Energy of the last SCF cycle; the same as `SCF TOTAL ENERGY` once converged. |
| `SCF ITERATIONS` | 8 | Number of SCF cycles. A large number warns of a hard case (radicals, stretched bonds). |
| `SCF TOTAL ENERGIES` | array (9) | The energy at every SCF cycle: handy for plotting convergence. |
| `NUCLEAR REPULSION ENERGY` | 9.194 | Repulsion between the nuclei only. Depends only on the geometry; a quick check that two structures match. |
| `ONE-ELECTRON ENERGY` | −123.05 | Kinetic energy of the electrons plus their attraction to the nuclei. |
| `TWO-ELECTRON ENERGY` | 37.85 | Electron–electron repulsion (Coulomb minus exchange). The total is nuclear + one-electron + two-electron. |
| `HF KINETIC ENERGY` | 75.82 | Kinetic energy of the electrons. |
| `HF POTENTIAL ENERGY` | −151.83 | All potential energy terms together. |
| `HF VIRIAL RATIO` | 2.0025 | −*V*/*T*. Exactly 2 for an exact wavefunction at an equilibrium geometry (the virial theorem); close to 2 is a good sign. |
| `SCF DIPOLE` | array (3) | Dipole moment vector at the SCF level, in $e\,a_0$. |
| `CURRENT DIPOLE` | array (3) | Dipole from the most recent job. |
| `PCM POLARIZATION ENERGY` | 0 | Solvent contribution from PCM; zero unless PCM is on (section 8). |
| `DD SOLVATION ENERGY` | 0 | Solvent contribution from the ddCOSMO/ddPCM models (`ddx` option); zero unless used. |
| `PE ENERGY` | 0 | Contribution from polarizable embedding (QM/MM); zero unless used. |

## 2. DFT Jobs

**Created by:** a DFT method such as `b3lyp`. Added to the SCF variables above. Values are for B3LYP/6-31G* water.

| Variable | Water value | What it means |
|---|---|---|
| `DFT TOTAL ENERGY` | −76.4087 | DFT energy, including any dispersion correction. |
| `DFT FUNCTIONAL TOTAL ENERGY` | −76.4087 | DFT energy from the functional alone, **without** dispersion. |
| `DFT XC ENERGY` | −7.564 | The exchange–correlation energy: the part that distinguishes one functional from another. |
| `DFT VV10 ENERGY` | 0 | Non-local VV10 dispersion term; only nonzero for functionals such as ωB97X-V. |
| `GRID ELECTRONS TOTAL` | 10 | Electrons counted by integrating the density on the DFT grid. Should match the real count; a big difference means the grid is too coarse. |
| `GRID ELECTRONS ALPHA`, `GRID ELECTRONS BETA` | 5, 5 | The same for each spin. |
| `XC GRID RADIAL POINTS`, `XC GRID SPHERICAL POINTS`, `XC GRID TOTAL POINTS` | 75, 302, 65 914 | Size of the integration grid. Set with the `dft_radial_points` and `dft_spherical_points` options. |

## 3. Dispersion Corrections (-D3, -D4)

**Created by:** a method with a dispersion suffix such as `b3lyp-d3bj`.

| Variable | Water value | What it means |
|---|---|---|
| `DISPERSION CORRECTION ENERGY` | −0.00057 | Grimme's empirical correction for London dispersion, added to the DFT energy. |
| `B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY` | −0.00057 | The same number, named after the method. |

The correction is tiny for one water molecule but matters for stacking, hydrogen-bonded complexes and conformations of larger molecules.

## 4. MP2 Jobs

**Created by:** `mp2`. Added to the SCF variables. Values are for MP2/cc-pVDZ water with `freeze_core`.

| Variable | Water value | What it means |
|---|---|---|
| `MP2 TOTAL ENERGY` | −76.2284 | MP2 energy: HF plus the MP2 correlation energy. |
| `MP2 CORRELATION ENERGY` | −0.2016 | The correlation energy alone: the energy HF misses because electrons avoid each other. |
| `CURRENT CORRELATION ENERGY` | −0.2016 | Correlation energy of the most recent correlated job. |
| `MP2 DOUBLES ENERGY` | −0.2016 | Part from pairs of electrons excited together; for closed-shell molecules it is all of it. |
| `MP2 SINGLES ENERGY` | ≈ 0 | Part from single excitations; zero for RHF references, nonzero for ROHF. |
| `MP2 SAME-SPIN CORRELATION ENERGY` | −0.0508 | Correlation between electrons of the same spin. |
| `MP2 OPPOSITE-SPIN CORRELATION ENERGY` | −0.1508 | Correlation between electrons of opposite spin, the larger part. |
| `SCS-MP2 TOTAL ENERGY`, `SCS-MP2 CORRELATION ENERGY` | −76.2247, −0.1979 | **Spin-component-scaled** MP2: the two spin parts reweighted (6/5 and 1/3), which usually improves on MP2 at no extra cost. |
| `CUSTOM SCS-MP2 TOTAL ENERGY`, `CUSTOM SCS-MP2 CORRELATION ENERGY` | −76.2247, −0.1979 | The same with weights set by the `mp2_os_scale` and `mp2_ss_scale` options; equal to SCS-MP2 by default. |

## 5. Coupled-Cluster Jobs (CCSD and CCSD(T))

**Created by:** `ccsd` or `ccsd(t)`. A CCSD(T) job also creates all the MP2 variables, since MP2 is its starting guess. Values are for CCSD(T)/cc-pVDZ water with `freeze_core`.

| Variable | Water value | What it means |
|---|---|---|
| `CCSD(T) TOTAL ENERGY` | −76.2410 | The "gold standard" energy used for benchmarks in this book. |
| `CCSD(T) CORRELATION ENERGY` | −0.2142 | CCSD(T) correlation energy. |
| `(T) CORRECTION ENERGY` | −0.0030 | The perturbative triples correction: CCSD(T) minus CCSD. |
| `CCSD TOTAL ENERGY`, `CCSD CORRELATION ENERGY` | −76.2380, −0.2112 | The same without the triples. |
| `CCSD SINGLES ENERGY`, `CCSD DOUBLES ENERGY` | 0, −0.2112 | Contributions of single and double excitations. |
| `CCSD SAME-SPIN CORRELATION ENERGY`, `CCSD OPPOSITE-SPIN CORRELATION ENERGY` | −0.0453, −0.1659 | Spin components, as for MP2. |
| `CCSD ITERATIONS` | 11 | Number of CCSD cycles. |
| `CC T1 DIAGNOSTIC` | 0.006 | A warning light. Values above about 0.02 suggest a single HF reference is a poor description (stretched bonds, some radicals) and CCSD(T) may be unreliable. |
| `CC D1 DIAGNOSTIC`, `CC NEW D1 DIAGNOSTIC`, `CC D2 DIAGNOSTIC` | 0.011, 0.011, 0.123 | Related warning lights, read in the same way: larger means less trustworthy. |

The diagnostics are worth a look in the *Bond Dissociation* notebook, where the stretched O–H bond pushes them up.

## 6. Geometry Optimizations

**Created by:** `optimize`. Gradients are calculated at every step, so they appear along with the energy variables.

| Variable | Water shape / value | What it means |
|---|---|---|
| `OPTIMIZATION ITERATIONS` | 3 | Number of optimization steps taken. |
| `CURRENT GRADIENT` | 3 × 3 array | Force on each atom (*x*, *y*, *z*) at the last step, in $E_\mathrm{h}/a_0$. Close to zero at a converged structure. |
| `SCF TOTAL GRADIENT`, `HF TOTAL GRADIENT` | 3 × 3 arrays | The same, named by method. |

## 7. Frequency Jobs and Thermochemistry

**Created by:** `frequency` (also `hessian`, and `vibanal_wfn`, which re-runs the analysis at a new temperature). Thermochemistry uses the temperature and pressure options `T` (298.15 K) and `P` (101 325 Pa). Values are for HF/6-31G* water.

| Variable | Water value | What it means | Used in |
|---|---|---|---|
| `ZPVE` | 0.02297 | **Zero-point vibrational energy**: energy the molecule keeps even at 0 K. | *Thermochemistry*, *Isotope Effects*, *Hydrogen Bonding* |
| `ZERO K ENTHALPY` | −75.9877 | Electronic energy + ZPVE. | *Thermochemistry* |
| `THERMAL ENERGY CORRECTION` | 0.02581 | ZPVE plus the translational, rotational and vibrational energy at *T*. | *Thermochemistry* |
| `THERMAL ENERGY` | −75.9849 | Electronic energy + thermal energy correction (*U*). | |
| `ENTHALPY CORRECTION` | 0.02675 | Thermal energy correction + *RT*. | *Thermochemistry* |
| `ENTHALPY` | −75.9840 | **Enthalpy** *H* = electronic energy + enthalpy correction. | *Thermochemistry*, *Acids and Bases* |
| `GIBBS FREE ENERGY CORRECTION` | 0.00538 | Enthalpy correction − *TS*. | *Thermochemistry* |
| `GIBBS FREE ENERGY` | −76.0053 | **Gibbs energy** *G* = electronic energy + Gibbs correction. | *Thermochemistry*, *Acids and Bases* |
| `CURRENT HESSIAN` | 9 × 9 array | Second derivatives of the energy (force constants) in $E_\mathrm{h}/a_0^2$. The frequencies come from this. | |
| `SCF TOTAL HESSIAN`, `HF TOTAL HESSIAN` | 9 × 9 arrays | The same, named by method. | |
| `CURRENT DIPOLE GRADIENT` | 9 × 3 array | How the dipole changes as each atom moves. The IR intensities come from this. | |
| `SCF DIPOLE GRADIENT` | 9 × 3 array | The same, named by method. | |

Entropy has no variable of its own. Calculate it from $S = (H - G)/T$.

## 8. Properties with `psi4.oeprop()`

**Created by:** `psi4.oeprop(wfn, "DIPOLE", "QUADRUPOLE", "MULLIKEN_CHARGES", ..., title="HF")`. The `title` becomes the first word of the name. Without a title, *Psi4* uses the method name.

| Variable | Where | What it means | Used in |
|---|---|---|---|
| `HF DIPOLE` (`<title> DIPOLE`) | global and wfn | Dipole moment vector in $e\,a_0$. The length of the vector is the dipole moment. | *Charges and Dipole*, *Solvation* (`GAS DIPOLE`, `PCM DIPOLE`) |
| `HF QUADRUPOLE` (`<title> QUADRUPOLE`) | global and wfn | Quadrupole moment (6 unique components). | |
| `MULLIKEN CHARGES` | **wfn only** | Partial charge on each atom by Mulliken's method. Simple but very basis-set dependent. | *Charges and Dipole* |
| `LOWDIN CHARGES` | **wfn only** | Löwdin charges: similar idea, somewhat less basis dependent. | *Charges and Dipole* |
| `LOWDIN SPINS` | **wfn only** | Unpaired spin on each atom (zero for closed-shell molecules). | |
| `MBIS CHARGES` | **wfn only** | Charges from dividing up the electron density (Minimal Basis Iterative Stockholder). The most stable of the three. | *Charges and Dipole* |
| `MBIS DIPOLES`, `MBIS QUADRUPOLES`, `MBIS OCTUPOLES` | **wfn only** | Atomic multipoles from the same partitioning. | |
| `MBIS VALENCE CHARGES`, `MBIS VALENCE WIDTHS` | **wfn only** | Details of the MBIS model of each atom. | |
| `MBIS RADIAL MOMENTS <R^2>`, `MBIS RADIAL MOMENTS <R^3>`, `MBIS RADIAL MOMENTS <R^4>` | **wfn only** | Measures of the size of each atom's share of the density. | |

## 9. Solvation with PCM

**Created by:** an SCF job with `pcm: True`. No new names appear, but some values change. For B3LYP/6-31G* water in PCM water:

| Variable | Value | What it means |
|---|---|---|
| `PCM POLARIZATION ENERGY` | −0.00998 | Stabilization from the polarized solvent (about −26 kJ/mol). |
| `SCF TOTAL ENERGY` | −76.4178 | **Includes** the PCM term. This is the energy to use for solvated species. |
| `DFT TOTAL ENERGY` | −76.4078 | Does **not** include the PCM term. |

## 10. Counterpoise and Many-Body Jobs

**Created by:** `psi4.energy("mp2", molecule=dimer, bsse_type=["cp", "nocp"])`, where the molecule has two fragments separated by `--`. Values are for the MP2/cc-pVDZ water dimer. **CP** is counterpoise-corrected (each monomer calculated with the partner's basis functions present as ghost atoms); **NOCP** is uncorrected.

| Variable | Dimer value | What it means | Used in |
|---|---|---|---|
| `CP-CORRECTED INTERACTION ENERGY` | −0.00626 (−16.4 kJ/mol) | **Binding energy** corrected for basis set superposition error. | *Hydrogen Bonding* |
| `NOCP-CORRECTED INTERACTION ENERGY` | −0.01176 (−30.9 kJ/mol) | Uncorrected binding energy, too negative because each monomer borrows its partner's basis functions. | *Hydrogen Bonding* |
| `CP-CORRECTED TOTAL ENERGY`, `NOCP-CORRECTED TOTAL ENERGY` | −152.463, −152.469 | Total energy of the dimer by each scheme. | |
| `CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY`, `NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY` | −0.00626, −0.01176 | Part of the energy due to pairs of fragments. For a dimer, this is the interaction energy. | |
| `CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY`, `CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY` | 0, −0.00626 | Running totals by number of bodies; matter for trimers and larger clusters. | |
| `CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY`, `CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY` | −152.457, −152.463 | The same for total energies. | |
| `NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY`, `NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY`, `NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY`, `NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY` | 0, −0.01176, −152.457, −152.469 | The same running totals without the counterpoise correction. | |
| `N-BODY (1)@(1) TOTAL ENERGY` | −76.2284 | Fragment 1 in its own basis set. Read `(fragments)@(basis)`. | |
| `N-BODY (1)@(1, 2) TOTAL ENERGY` | −76.2291 | Fragment 1 in the basis of both fragments (partner as ghost). Lower, because of the borrowed basis functions. | |
| `N-BODY (2)@(2) TOTAL ENERGY`, `N-BODY (2)@(1, 2) TOTAL ENERGY` | −76.2284, −76.2333 | The same for fragment 2. | |
| `N-BODY (1, 2)@(1, 2) TOTAL ENERGY` | −152.469 | The whole dimer. | |
| `NBODY NUMBER` | 5 | Number of separate calculations run. | |
| `CALCINFO NATOM`, `CALCINFO NFR`, `CALCINFO NMC` | 6, 2, 1 | Number of atoms, fragments and model chemistries in the job. | |

After this job `CURRENT ENERGY` holds the interaction energy of the first `bsse_type`, not a total energy.

## 11. SAPT Jobs

**Created by:** `psi4.energy("sapt0", molecule=dimer)`. **SAPT** (symmetry-adapted perturbation theory) calculates the interaction energy of two fragments directly and splits it into physical parts. Values are for SAPT0/jun-cc-pVDZ on the water dimer.

The four main terms are the ones to use:

| Variable | Dimer value | What it means | Used in |
|---|---|---|---|
| `SAPT ELST ENERGY` | −0.01406 (−36.9 kJ/mol) | **Electrostatics**: attraction between the two unperturbed charge distributions. | *Hydrogen Bonding* |
| `SAPT EXCH ENERGY` | +0.01134 (+29.8 kJ/mol) | **Exchange**: Pauli repulsion as the electron clouds overlap. Always positive. | *Hydrogen Bonding* |
| `SAPT IND ENERGY` | −0.00335 (−8.8 kJ/mol) | **Induction**: each molecule polarizing the other, including charge transfer. | *Hydrogen Bonding* |
| `SAPT DISP ENERGY` | −0.00203 (−5.3 kJ/mol) | **Dispersion**: correlated motion of electrons (London forces). | *Hydrogen Bonding* |
| `SAPT TOTAL ENERGY` | −0.00809 (−21.2 kJ/mol) | Sum of the four: the SAPT interaction energy. Also in `CURRENT ENERGY`. | *Hydrogen Bonding* |

The same five numbers appear again under the level's own name, and in two variants:

| Variables | What they are |
|---|---|
| `SAPT0 ELST ENERGY`, `SAPT0 EXCH ENERGY`, `SAPT0 IND ENERGY`, `SAPT0 DISP ENERGY`, `SAPT0 TOTAL ENERGY` | Identical to the `SAPT ...` values, named for the SAPT0 level. Higher levels (SAPT2, SAPT2+) have their own names. |
| `SCS-SAPT0 ELST ENERGY`, `SCS-SAPT0 EXCH ENERGY`, `SCS-SAPT0 IND ENERGY`, `SCS-SAPT0 DISP ENERGY`, `SCS-SAPT0 TOTAL ENERGY` | Spin-component-scaled SAPT0: dispersion reweighted as in SCS-MP2. |
| `SSAPT0 ELST ENERGY`, `SSAPT0 EXCH ENERGY`, `SSAPT0 IND ENERGY`, `SSAPT0 DISP ENERGY`, `SSAPT0 TOTAL ENERGY` | Scaled SAPT0 (sSAPT0): exchange terms scaled to correct a known error of SAPT0. Recommended with the jun-cc-pVDZ basis set. |
| `SAPT2 DISP ENERGY` | Also set by a SAPT0 job; equal to the SAPT0 dispersion here. |

The rest are the individual terms that build up the four main ones. They are for specialists:

| Variables | What they are |
|---|---|
| `SAPT ELST10,R ENERGY` | First-order electrostatics (all of `SAPT ELST` at SAPT0). |
| `SAPT EXCH10 ENERGY`, `SAPT EXCH10(S^2) ENERGY` | First-order exchange, exact and in the common $S^2$ approximation. |
| `SAPT IND20,R ENERGY`, `SAPT EXCH-IND20,R ENERGY` | Second-order induction and its exchange partner. |
| `SAPT HF TOTAL ENERGY`, `SAPT HF(2) ENERGY`, `SAPT HF(2) ALPHA=0.0 ENERGY` | δHF: higher-order induction estimated from a dimer HF calculation, added to the induction term. |
| `SAPT CT ENERGY` | Charge-transfer part of induction. |
| `SAPT DISP20 ENERGY`, `SAPT EXCH-DISP20 ENERGY` | Second-order dispersion and its exchange partner. |
| `SAPT SAME-SPIN DISP20 ENERGY`, `SAPT SAME-SPIN EXCH-DISP20 ENERGY` | Same-spin parts of those, used for SCS-SAPT0. |
| `SAPT EXCHSCAL`, `SAPT EXCHSCAL1`, `SAPT EXCHSCAL3`, `SAPT ALPHA` | Exchange-scaling factors used for sSAPT0 (no units). |
| `SAPT ENERGY` | Unused here (0). |

The wavefunction also holds the same terms with short names (`E ELST10`, `E EXCH10`, `E IND20`, `E DISP20`, `E SAPT0`, `E SCS-SAPT0`, ...), which `psi4.core.variables()` doesn't list.

## Summary
- PSI variables store every result of a calculation under a fixed name; read them with `psi4.variable()` or `wfn.variable()`.
- `psi4.core.variables()` lists what the jobs run so far have made. Each new job overwrites the global store, so keep the `wfn` for results you need later.
- Names follow patterns: `<METHOD> TOTAL ENERGY`, `<METHOD> CORRELATION ENERGY`, `CURRENT ...`.
- Atomic charges and detailed SAPT terms are stored only on the wavefunction.
- Energies are in hartree and dipoles in $e\,a_0$.

```{warning} AI Statement
:class: dropdown
This page was written by Claude (Anthropic) using Claude Code, at my request. Claude ran twelve typical *Psi4* 1.11 jobs on water and a water dimer, recorded the PSI variables each created, and checked their meanings against the PSI variables glossary in the *Psi4* manual. I selected, edited and checked all content.
```
