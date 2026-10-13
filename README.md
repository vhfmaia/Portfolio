# Advanced Abaqus User Subroutines Portfolio

[![Status](https://img.shields.io/badge/Status-Portfolio%20Verified%20%26%20In%20Work-blue.svg)](#)
[![FEA: Abaqus/Standard](https://img.shields.io/badge/FEA-Abaqus%2FStandard-005691.svg)](https://www.3ds.com/)
[![Language](https://img.shields.io/badge/Language-Fortran%202008%20%2F%20F90-734f96.svg)](https://fortran-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

<div align="justify">
A repository of production-grade user subroutines for <b>Abaqus/Standard</b> and <b>3DEXPERIENCE SIMULIA</b>, written in <b>Fortran</b> for high-performance structural, fatigue, and manufacturing simulations.
</div>

---

## Portfolio Projects & Subroutine Implementations

| ID | Subroutine | Physics / Type | Engineering Application | Status |
| :---: | :--- | :--- | :--- | :---: |
| **01** | [`DLOAD`](01_dload_gear_fatigue/README.md) | Moving Surface Load / Quasi-Static | Drive + coast rolling load and fatigue framework for a welded helical crown gear and differential carrier | ✅ Verified |
| **02** | [`DFLUX`](02_dflux_laser_welding/README.md) | Moving Conical Flux / Thermo-Mechanical | Sequentially coupled laser welding with phase transformations and residual stress | ✅ Verified |
| **03** | [`UINTER`](03_uinter_welding_interface/README.md) | Surface Interaction / BTR Hot Cracking | 3-State contact interface (Raw → Welded → Cracked) with BTR criteria | ⚠️ Exploratory (R&D Backlog) |
| **04** | [`UMAT`](04_umat_cohesive_welding/README.md) | Cohesive Zone Modeling / Fracture Energy *G*<sub>c</sub> | Mesh-independent 3-state thermo-mechanical weld interface (`COH3D8`) | ⚠️ In Work (Testing Phase) |
| **05** | [`USDFLD`](05_usdfld_shaft_straightening/README.md) | Radial Field Variable / Nonlinear Kinematic Hardening | Elastoplastic straightening & springback of case-hardened transmission shaft (DOI: [10.1007/s00170-014-5688-8](https://doi.org/10.1007/s00170-014-5688-8)) | ✅ Verified |

---

### Project Showcase: 01. Helical Crown Gear & Carrier Assembly (`DLOAD`)

<p align="center">
  <a href="01_dload_gear_fatigue/README.md">
    <img src="99_images_and_miscelania/gear_torque_simulation.gif" alt="Abaqus DLOAD Welded Crown Gear and Carrier Assembly Simulation" width="100%" />
  </a>
</p>

<div align="justify">
Finite element simulation (1.49 million equations) of an automotive helical crown gear (20MC5) laser-welded with Inconel filler wire to its differential carrier (EN-GJS-600-10), with a 100° incomplete-penetration weld flaw. The tooth contact load is applied analytically by a <code>DLOAD</code> subroutine and rolled across the full 360° revolution, in drive and in coast (regeneration), without meshing the pinion.
</div>

#### Key Features & Verification:
- **Drive and Coast Load Cases:** active flank and torque selected per analysis step, covering manual drivelines and hybrid / electric drivelines with torque reversal.
- **Constant Transmitted Torque:** the pressure field integrates to $+270.03\text{ N}\cdot\text{m}$ (drive) and $-270.11\text{ N}\cdot\text{m}$ (coast) over the flank surfaces, with ripple below $0.15\%$ through the revolution.
- **Finite Element Load Path:** reaction torque at the bearing reference node $-270.13\text{ N}\cdot\text{m}$ against $-270.0$ applied ($0.05\%$), averaged over the full revolution; transverse and axial bearing forces within $0.1\%$ of gear force theory.
- **Welded Assembly:** bonded weld seam over 260°, frictional contact across the unbonded flaw, kinematic couplings at the bearing seats.
- **Fatigue Method:** stress-range evaluation per critical point with S-N curves and Palmgren-Miner damage summation for an OEM load spectrum.

<div align="justify">
Detailed formulation, kinematic equations, FEA verification, and post-processing scripts are available in <a href="01_dload_gear_fatigue/README.md">01_dload_gear_fatigue/</a>.
</div>

---

### Project Showcase: 02. Circumferential Laser Welding (`DFLUX`)

<p align="center">
  <a href="02_dflux_laser_welding/README.md">
    <img src="99_images_and_miscelania/conical_heat_source_improved_plot.svg" alt="Abaqus DFLUX 3D Conical Gaussian Heat Source and Welding Power Schedule" width="100%" />
  </a>
</p>

<div align="justify">
Sequentially coupled thermo-mechanical analysis of circumferential laser welding (<i>R</i> = 15 mm, <i>P</i><sub>abs</sub> = 1080 W, <i>v</i> = 1200 mm/min) with moving conical Gaussian heat source, metallurgical phase transformations (<code>ABQ_PHASE_TRANS</code>), hybrid elements (<code>C3D8H</code>), continuous annealing at 1530 °C, and residual stress prediction for 16MnCr5 case-hardening steel.
</div>

#### Key Thermal, Metallurgical & Mechanical Highlights:
- **3D Conical Gaussian Heat Source:** Analytically normalized peak flux *Q*<sub>0</sub> = 4.6935 × 10⁵ mW/mm³ with 10° linear power ramp-up, 360° steady weld, and 10° crater-filling ramp-down overlap.
- **Unified Material Architecture:** Single `Material_16MnCr5.inp` shared by both thermal and mechanical simulations, ensuring seamless data consistency.
- **Metallurgical Kinetics & Solidification (BTR):** Tracks liquidus/solidus melt pool kinetics and solid-state phase decomposition (Ferrite-Pearlite, Austenite, Bainite, Martensite via JMA and Koistinen-Marburger kinetics).
- **Quadratic Unsymmetric Thermal Convergence:** Direct unsymmetric solver formulation (`UNSYMM=YES`) resolving non-symmetric metallurgical tangent matrices, localized reference flux controls (`*CONTROLS, PARAMETERS=FIELD`) eliminating cold-element dilution, and uninhibited Newton-Raphson steps without line search under-relaxation.
- **Sequentially Coupled Structural Analysis:** Mitigates high-temperature volumetric locking via hybrid formulation elements (`C3D8H`), hydrostatic gravity stabilization, and finite-strain kinematics (`NLGEOM=YES`).

<div align="justify">
Detailed formulation, verification routines, simulation input decks, and automated post-processing pipeline are available in <a href="02_dflux_laser_welding/README.md">02_dflux_laser_welding/</a>.
</div>

---

### Project Showcase: 03 & 04. Advanced Interface Modeling Evolution

<div align="justify">
The portfolio demonstrates a complete engineering R&D cycle comparing two distinct methodologies to capture the 3 physical states of a welded interface (<b>Raw → Welded → Cracked</b>):
</div>

<p align="center">
  <img src="99_images_and_miscelania/uinter_vs_czm_comparison.svg" alt="Architectural Evolution from Contact Interaction UINTER to Cohesive Zone Modeling UMAT" width="100%" />
</p>

<div align="justify">
<ul>
  <li><b><a href="03_uinter_welding_interface/README.md">03. Exploratory Interface (UINTER)</a>:</b> Formulates the 3-state state machine on contact surfaces. Benchmarking identified critical mesh-dependency due to displacement jumps vs strain rates, prompting its classification as an exploratory R&D backlog module.</li>
  <li><b><a href="04_umat_cohesive_welding/README.md">04. Cohesive Zone Interface (UMAT)</a>:</b> Under active development and testing to evaluate mesh-independent fracture energy <i>G</i><sub>c</sub> regularization. Uses an exact 3 × 3 analytical <code>DDSDDE</code> Jacobian to investigate convergence without contact chattering or SDI cutbacks (strictly a research prototype under evaluation, not for production use).</li>
</ul>
</div>

---

## Directory Layout

```text
.
├── 01_dload_gear_fatigue/
│   ├── dload_crown_helical.f                  # User Subroutine (Modern Fortran DLOAD + UEXTERNALDB)
│   ├── Gear_Shaft_assembly_simulation.inp     # Master Input Deck: drive revolution
│   ├── Gear_Shaft_assembly_drive_coast.inp    # Master Input Deck: drive + coast revolutions
│   ├── Materials_Assembly.inp                 # 20MC5 and GJS600-10 material cards
│   ├── Gear_Shaft_geometry.step               # CAD geometry: crown gear + differential carrier
│   ├── compile_vtk_to_gif.py                  # Multi-viewport visualization compiler
│   ├── odb_to_vtk.py                          # Abaqus Python ODB to VTK extractor
│   ├── run_pipeline.sh                        # Headless Linux / HPC batch execution script
│   └── README.md                              # Comprehensive Technical Report
├── 02_dflux_laser_welding/
│   ├── Disk_heatsource_TH.inp                 # Master Thermal Analysis Deck (DC3D8)
│   ├── Disk_heatsource_ME.inp                 # Master Mechanical Analysis Deck (C3D8H, Sequentially Coupled)
│   ├── Geometry_TH.inp                        # Thermal Mesh Deck (DC3D8 Solid Bricks)
│   ├── Geometry_ME.inp                        # Mechanical Hybrid Mesh Deck (C3D8H, NSET_BASE)
│   ├── Material_16MnCr5.inp                   # Unified Thermo-Elasto-Plastic & Phase Transformation Deck
│   ├── dflux_disk_conical_gaussian.f          # User Subroutine (Modern Fortran TDC Conical Model)
│   ├── compile_vtk_to_gif.py                  # Multi-Viewport Animated Report Compiler
│   ├── odb_to_vtk.py                          # Abaqus Python ODB to VTK Extractor & ZIP Packager
│   ├── run_pipeline.sh                        # Headless Linux / HPC Batch Execution Pipeline
│   └── README.md                              # Comprehensive Technical Report
├── 03_uinter_welding_interface/
│   ├── uinter_btr_hot_cracking.f              # User Subroutine (UINTER State Machine in Fortran)
│   ├── dummy_btr_contact_patch.inp            # 2-Element Contact Patch Benchmark Deck
│   └── README.md                              # Technical Report & R&D Backlog Documentation
├── 04_umat_cohesive_welding/
│   ├── umat_cohesive_welding_btr.f            # Prototype User Material Subroutine (Modern Fortran)
│   ├── cohesive_weld_btr_verification.inp     # 3-Element Cohesive Zone Benchmark Deck (COH3D8)
│   └── README.md                              # Comprehensive Technical Report & Formulation
├── 99_images_and_miscelania/
│   ├── gear_torque_simulation.gif             # Full 360-degree three-viewport animation (6 FPS)  
│   ├── conical_heat_source_improved_plot.svg  # 4-Panel 3D Conical Heat Source & Trajectory Plot (SVG)
│   ├── conical_heat_source_schematic.svg      # Analytical Heat Source Geometry Schematic (SVG)
│   ├── sequentially_coupled_architecture.svg  # Multi-Physics Data Flow & Pipeline Architecture (SVG)
│   ├── uinter_3state_state_machine.svg        # UINTER 3-State Contact State Machine & Won Criteria (SVG)
│   ├── cohesive_3state_state_machine.svg      # Cohesive Zone UMAT Architecture & Softening Law (SVG)
│   ├── uinter_vs_czm_comparison.svg           # Interface Modeling Evolution & Bottlenecks (SVG)
│   └── conical_gaussian_3d.png                # Conical Heat Source Surface Elevation (PNG)
├── LICENSE
└── README.md
```

---

## Disclaimer

<div align="justify">
This repository is published for professional portfolio and academic demonstration purposes. All CAD dimensions, process parameters, and metallurgical kinetics presented herein are synthetic, open-literature benchmarks designed to illustrate advanced FEA and Fortran programming methodologies, and do not represent any confidential data.
</div>

---

## Development Approach & Use of AI

<div align="justify">
The engineering content of this portfolio is the author's own: the choice of problems, the physical and numerical modelling strategy, the finite element models (geometry, meshes, materials, contact, welds and boundary conditions), the simulation campaigns on Abaqus / 3DEXPERIENCE and the judgement of whether a result is physically right. AI coding assistants were used as a productivity tool, under the author's direction, to speed up the writing and adaptation of Fortran subroutines, pre- and post-processing scripts, numerical cross-checks and documentation. Generated code is not taken on trust: projects marked as verified are checked against an analytical solution, an independent numerical integration or the solver's own reaction output, as documented in each project, and work still in progress is labelled as such in the table above.
</div>

---

## Author & Contact

**Victor Maia**  
- **Email:** [vhfm08@gmail.com](mailto:vhfm08@gmail.com)  
- **GitHub:** [@vhfmaia](https://github.com/vhfmaia)  
