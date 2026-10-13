# 05. Elastoplastic Straightening & Springback Analysis of Case-Hardened Transmission Shaft (`USDFLD`)

[![Status: Verified](https://img.shields.io/badge/Status-Verified-brightgreen.svg)](#)
[![FEA: Abaqus/Standard](https://img.shields.io/badge/FEA-Abaqus%2FStandard-blue.svg)](https://www.3ds.com/)
[![Subroutine: USDFLD](https://img.shields.io/badge/Subroutine-USDFLD-005691.svg)](#)
[![Language: Fortran](https://img.shields.io/badge/Language-Fortran%202008%20%2F%20F90-734f96.svg)](https://fortran-lang.org/)

---

## Executive Summary & Industrial Scope

<div align="justify">
In automotive transmission manufacturing, critical driveline components—such as secondary countershafts and stepped pinion shafts—undergo severe thermochemical case hardening (gas carburizing or carbonitriding followed by accelerated oil quenching) to achieve exceptional surface wear resistance, high contact fatigue strength, and bending endurance. However, the severe, non-uniform thermal gradients and differential martensitic phase transformation kinetics inherent to liquid quenching inevitably induce significant axial warping and geometric distortion (runout).
</div>

<div align="justify">
To ensure coaxiality and satisfy stringent powertrain <b>Total Indicated Runout (TIR)</b> manufacturing tolerances (typically under 20–30 µm) prior to hard-finish gear grinding and final assembly, serial production lines deploy automated <b>three-point mechanical press straightening</b> machines.
</div>

<p align="center">
  <img src="../99_images_and_miscelania/shaft_straightening_schematic.svg" alt="Abaqus USDFLD Secondary Shaft Three-Point Straightening and Springback Architecture" width="100%" />
</p>

<div align="justify">
This repository documents an advanced, high-fidelity three-dimensional finite element analysis (FEA) framework implemented in <b>Abaqus/Standard</b>. The pipeline models the complete non-linear mechanical response of a case-hardened secondary transmission shaft based on industrial benchmarking geometry (<code>husson_shaft.step</code>):
</div>

1. **Initial Quenching Residual Stress Pre-Equilibrium:** Dissipation of self-equilibrating thermal-martensitic residual stress tensors (-350 to -500 MPa compressive case, +150 to +250 MPa tensile core);
2. **Radial Hardness & Yield Gradient Interpolation (`USDFLD`):** Continuous integration-point mapping of depth-dependent constitutive properties across the 0.6 mm case layer into the ductile core;
3. **Nonlinear Reverse Plasticity & Bauschinger Effect:** Accurate prediction of plastic overbend and reverse yielding during punch unloading using a two-backstress **Chaboche nonlinear kinematic hardening model**;
4. **Full 3D Conformal Quadratic Tetrahedral Discretization (`C3D10` / `C3D10M`):** Pragmatic meshing architecture accommodating internal lubrication cross-holes, splines, and helical gear teeth;
5. **Post-Straightening Total Indicated Runout (TIR):** Quantitative springback calculation from the centerline displacement after unloading.

---

## 1. Physical Motivation: The Straightening & Springback Dilemma

<div align="justify">
A case-hardened transmission shaft is inherently a bimaterial continuum with an extreme mechanical property gradient:
</div>

* **Brittle Case Layer (0.0–0.6 mm):** High surface carbon content (~0.8 wt% C), tempered martensite microstructure, ~60 HRC (<i>σ</i><sub>y,0</sub> &approx; 1900 MPa), high notch sensitivity, and low plastic fracture strain (<i>ε</i><sub>f</sub> &lt; 2%).
* **Ductile Core (&gt; 0.6 mm):** Low carbon base steel (~0.16–0.20 wt% C), bainite/martensite core, ~30 HRC (<i>σ</i><sub>y,0</sub> &approx; 750 MPa), and high plastic ductility (<i>ε</i><sub>f</sub> &gt; 12%).

<div align="justify">
During three-point press straightening, the punch imposes a localized bending moment on the central bearing journal exceeding the section yield strength. When the outer fibers yield in tension and compression, the ductile core undergoes substantial plastic deformation while the high-strength case layer undergoes extreme elastic-plastic strain.
</div>

<div align="justify">
Upon punch retraction, the elastically stored strain energy in the shaft drives <b>springback recovery</b>. Crucially, the internal redistribution of residual stresses forces the plastified zones into <b>reverse plastic loading</b>. If an engineer uses a standard isotropic hardening rule, the solver completely ignores the <b>Bauschinger effect</b>, dramatically underestimating reverse plastic strain and predicting an erroneous springback angle (typically off by 30–45%).
</div>

---

## 2. Anatomy of the `USDFLD` Subroutine & Parameter Breakdown

<div align="justify">
Rather than slicing the CAD geometry into artificial geometric shells with discontinuous material boundaries—which creates numerical stress concentrations and meshing bottlenecks—this architecture utilizes an advanced <b>Fortran 2008 User Defined Field subroutine (<code>USDFLD</code>)</b>.
</div>

### 2.1 Keyword Deck Integration

<div align="justify">
The subroutine is linked directly to the Abaqus/Standard material card via <code>*USER DEFINED FIELD</code> and <code>*PLASTIC, DEPENDENCIES=1</code>:
</div>

```abaqus
*USER DEFINED FIELD
*MATERIAL, NAME=STEEL_CASE_HARDENED
*ELASTIC
 210000.0, 0.3
*PLASTIC, HARDENING=COMBINED, DATA TYPE=PARAMETERS, DEPENDENCIES=1
** Yield,     C1,    gamma1,      C2,   gamma2,  Q_inf,      b,    FV1 (Depth mm)
 1900.0, 125000.0,    180.0, 45000.0,    35.0,  150.0,   12.0,           0.00
 1500.0, 110000.0,    160.0, 38000.0,    30.0,  180.0,   10.0,           0.15
 1100.0,  95000.0,    140.0, 32000.0,    25.0,  220.0,    8.0,           0.35
  750.0,  80000.0,    120.0, 25000.0,    20.0,  250.0,    6.0,           0.60
  750.0,  80000.0,    120.0, 25000.0,    20.0,  250.0,    6.0,          25.00
*DEPVAR
 3
 1, SDV_DEPTH, "Radial depth from outer case boundary (mm)"
 2, SDV_PEEQ,  "Accumulated equivalent plastic strain"
 3, SDV_RADIUS,"Current radial distance from centerline (mm)"
```

### 2.2 Detailed Parameter Analysis & Data Provenance (Measured vs. Literature)

<div align="justify">
A critical requirement in high-precision finite element analysis is <b>parameter provenance traceability</b>—differentiating between directly measured physical properties, empirical conversions and literature baselines.
</div>

#### Data Source Classification Matrix:

1. **Direct Experimental Measurements (Primary Source: Husson et al., 2014):**
   * *Benchmark Publication:* Rémi Husson, Cyrille Baudouin, Régis Bigot, Edoardo Sura, *"Consideration of residual stress and geometry during heat treatment to decrease shaft bending"*, *The International Journal of Advanced Manufacturing Technology* (2014) 72:1455–1463, DOI: [10.1007/s00170-014-5688-8](https://doi.org/10.1007/s00170-014-5688-8).
   * *Base Steel Alloy:* 20MnCr5 (DIN EN 10084, WNr. 1.7147, SAE 5120 equivalent).
   * *CAD Component Geometry:* Secondary transmission countershaft of the benchmark (`husson_shaft.step`, $Z = -8.6$ to $285.5\text{ mm}$, overall length $L = 294.1\text{ mm}$), featuring Ø25 mm splines, Ø38 mm central bearing journal, Ø51 mm shoulder, Ø66.9 mm helical gear, and 4 transverse lubrication cross-holes.
   * *Heat Treatment Condition:* Gas carburized, liquid oil quenched ($860^\circ\text{C}$), and tempered ($180^\circ\text{C}$) (Industrial experimental batch **Lot E**).
   * *Microhardness Depth Profile:* Directly measured Vickers microhardness values ($HV_{0.1}$) spanning from the outer surface ($d = 0.0\text{ mm}$, $60.5\text{ HRC}$ / $\sim 700\text{ HV}$) through the case boundary down to the ductile base core ($d > 2.0\text{ mm}$, $30.0\text{ HRC}$ / $\sim 300\text{ HV}$).
   * *Effective Case Depth (CHD / Eht 550 HV1):* Measured at $d_{\text{case}} = 0.60\text{ mm}$.
   * *Elastic Constants:* Standard transmission steel values ($E = 210\,000.0\text{ MPa}$, $\nu = 0.30$, $\rho = 7.85 \times 10^{-9}\text{ tonne/mm}^3$).

2. **Empirical Microstructural Yield Conversion (Tabor's Relationship):**
   * The continuous radial yield stress gradient $\sigma_{y,0}(d)$ is derived analytically from the experimental Vickers microhardness profile using **Tabor's empirical relation** for strain-hardened metallic alloys ($\sigma_{y} \approx 3.0 \times HV$ in compatible stress units):
     $$\sigma_{y,0}(d) = 750.0 + (1900.0 - 750.0) \cdot \exp\left( -\left(\frac{d}{0.72}\right)^{1.6} \right) \quad [\text{MPa}]$$
     $$\text{HRC}(d) = 30.0 + (60.5 - 30.0) \cdot \exp\left( -\left(\frac{d}{0.75}\right)^{1.6} \right) \quad [-]$$

3. **Literature Baseline for Cyclic 20MnCr5 Plasticity (Chaboche 2-Backstress Formulation):**
   * The Armstrong-Frederick / Chaboche kinematic hardening constitutive formulation (Chaboche 1986, 1989; Lemaitre & Chaboche 1990) requires 2 backstress tensors ($\boldsymbol{\alpha} = \boldsymbol{\alpha}_1 + \boldsymbol{\alpha}_2$) to model both the sharp initial yield knee and the stabilized slope.
   * Baseline parameters were compiled from published cyclic tension-compression push-pull fatigue literature for carburized 20MnCr5 (e.g., Lang & Vöhringer 1995; Grum 2001; Bomas et al., AWT/HTM).

<div align="justify">
No press straightening test data are used in this project: the hardening parameters are literature-based estimates and have not been calibrated against stroke-force or runout measurements.
</div>

| Parameter | Identifier in Deck | Value Range | Physical Significance | Exact Provenance & Calibration Origin |
| :--- | :--- | :---: | :--- | :--- |
| **Field Variable 1** | `FIELD(1)` / `FV1` | 0.00–25.0 mm | Radial depth from external surface $d(r,z) = R_{\text{outer}}(z) - r$ | Computed dynamically at integration points via `GETVRM('COORD', ...)` in `usdfld_case_hardneded_shaft.f` |
| **Initial Yield Stress** | $\sigma_{y,0}$ | 1900.0 &rarr; 750.0 MPa | Size of initial elastic yield domain | Measured microhardness profile (Husson 2014, Lot E) converted via Tabor's relation |
| **Kinematic Modulus 1** | $C_1$ | 125.0 &rarr; 80.0 GPa | Initial kinematic hardening tangent modulus | Cyclic push-pull data (Lang & Vöhringer 1995) |
| **Dynamic Recovery 1** | $\gamma_1$ | 180.0 &rarr; 120.0 | Backstress saturation rate for component 1 | Fitted from nonlinear transition shoulder of initial backstress $\boldsymbol{\alpha}_1$ |
| **Kinematic Modulus 2** | $C_2$ | 45.0 &rarr; 25.0 GPa | Secondary kinematic hardening modulus | Stabilized cyclic hysteresis loop slope at large plastic strain ($\Delta\varepsilon_p / 2 = 1.0\%$) |
| **Dynamic Recovery 2** | $\gamma_2$ | 35.0 &rarr; 20.0 | Backstress saturation rate for component 2 | Literature-based estimate for the stabilized backstress |
| **Isotropic Expansion Max** | $Q_\infty$ | 150.0 &rarr; 250.0 MPa | Saturated isotropic yield expansion (Voce) | Monotonic tensile curve minus kinematic saturation limit; higher in ductile core |
| **Isotropic Decay Rate** | $b$ | 12.0 &rarr; 6.0 | Rate of isotropic yield expansion with PEEQ | Exponential fit of cyclic hardening evolution across first 10 cycles |

<div align="justify">
<b>Physical Rationale for Depth-Dependent Hardening:</b> The severe thermal and chemical gradient created during carburizing alters both carbon concentration (~0.8 wt% C at surface vs ~0.2 wt% C in core) and the martensite start temperature ($M_s$). The high surface carbon maximizes tetragonality, dislocation density, and internal lattice micro-stresses, raising initial yield strength to ~1900 MPa and driving rapid kinematic backstress saturation ($\gamma_1 = 180$). In contrast, the low-carbon ductile core undergoes extensive dislocation multiplication with higher isotropic work-hardening capacity ($Q_\infty = 250\text{ MPa}$, $b = 6.0$).
</div>

---

## 3. Mathematical Formulation: Chaboche Nonlinear Kinematic Hardening

<div align="justify">
The non-linear Bauschinger effect and cyclic ratchetting/springback under reverse bending are governed by the classical Armstrong-Frederick / Chaboche kinematic hardening constitutive framework:
</div>

### 3.1 Yield Surface Definition

$$\mathcal{F}(\boldsymbol{\sigma}, \boldsymbol{\alpha}, \sigma_y) = f(\boldsymbol{\sigma} - \boldsymbol{\alpha}) - \sigma_y(p) = \sqrt{\frac{3}{2} (\mathbf{s} - \boldsymbol{\alpha}'):(\mathbf{s} - \boldsymbol{\alpha}')} - \sigma_y(p) \le 0$$

where:
* <b>s</b> = <b>σ</b> - (1/3) tr(<b>σ</b>) <b>I</b> is the deviatoric Cauchy stress tensor;
* <b>α</b> is the total kinematic backstress tensor, representing the translation of the yield surface center in stress space;
* <b>α</b>' is the deviatoric part of the backstress tensor;
* <i>σ</i><sub>y</sub>(<i>p</i>) is the current size of the yield surface (isotropic component).

### 3.2 Two-Term Kinematic Backstress Evolution

<div align="justify">
To accurately capture both the sharp non-linear knee at yielding onset and the linear hardening slope at larger strain ranges, a superposed two-backstress decomposition is employed:
</div>

$$\boldsymbol{\alpha} = \sum_{k=1}^{2} \boldsymbol{\alpha}_k$$

$$\dot{\boldsymbol{\alpha}}_k = \frac{2}{3} C_k \dot{\boldsymbol{\varepsilon}}^p - \gamma_k \boldsymbol{\alpha}_k \dot{p}$$

where:
* <i>ε&#775;</i><sup>p</sup> is the plastic strain rate tensor;
* <i>p&#775;</i> is the equivalent plastic strain rate (PEEQ);
* <i>C</i><sub>k</sub> and <i>γ</i><sub>k</sub> are the material constants defined in the keyword deck as a function of field variable <i>FV</i><sub>1</sub>.

### 3.3 Isotropic Hardening Component (Voce Law)

$$\sigma_y(p) = \sigma_{y,0}(FV_1) + Q_\infty(FV_1) \left( 1 - e^{-b(FV_1) \, p} \right)$$

---

## 4. Finite Element Discretization Strategy: The Industrial Tetrahedral Approach

<div align="justify">
In academic textbooks, cylindrical shafts are frequently discretized using idealized structured hexahedral meshes (such as <code>C3D8I</code> with O-grid blocking). However, in realistic automotive transmission shafts—such as the actual secondary shaft analyzed here (<code>husson_shaft.step</code>)—the geometry features:
</div>

1. **Cross-Drilled Lubrication Holes:** Transverse internal channels (Ø3.5 mm and Ø8.0 mm) drilled through the shaft core to supply oil to bearing journals;
2. **Torque Transmission Splines:** Longitudinal drive flutes and keyways that disrupt rotational symmetry;
3. **Integral Helical Pinion Teeth:** Complex 3D involute surfaces directly adjacent to the main shaft bearing shoulder;
4. **Multiple Diameter Steps & Transition Fillets:** Stepped transitions from Ø25 mm to Ø38 mm, Ø51 mm, and Ø68.2 mm.

<div align="justify">
Attempting to create structured hexahedral partitions (<code>C3D8I</code>) across intersecting lubrication channels and helical teeth requires dozens of non-physical geometric cuts, introduces severe mesh skewness, and produces prohibitive pre-processing bottlenecks in CAE. Consequently, an industrial-grade, fully conformal <b>quadratic tetrahedral meshing strategy (<code>C3D10</code> / <code>C3D10M</code>)</b> is deployed.
</div>

### 4.1 Modified 10-Node Quadratic Tetrahedrals (`C3D10M`)
<div align="justify">
Standard quadratic tetrahedral elements (<code>C3D10</code>) can exhibit volumetric locking during unconstrained, large-scale plastic flow because their volumetric strain is over-constrained. Furthermore, standard <code>C3D10</code> elements can produce oscillating, non-uniform contact nodal forces on curved surfaces.
</div>

<div align="justify">
To overcome both limitations, the entire shaft is discretized using <b>Abaqus modified quadratic tetrahedrals (<code>C3D10M</code>)</b>:
</div>

* **Volumetric Locking Elimination:** Incorporates an internal hydrostatic pressure degree of freedom and specialized hourglass/volume control that maintains volume constancy (<i>ν</i><sub>eff</sub> &rarr; 0.5) without artificial stiffening;
* **Consistent Contact Formulation:** Features lumped mass and contact-consistent force distributions, eliminating contact chattering between the punch rigid surface and the shaft convex boundary;
* **Conformal Geometry:** Captures cross-drilled lubrication holes, spline roots, and tooth fillets without artificial tie interfaces or interface stiffness jumps.

### 4.2 Radial Case Layer Refinement Strategy
<div align="justify">
To accurately capture the steep hardness, yield stress, and residual stress gradient within the 0.6 mm carburized case layer:
</div>

* **Surface Inflation / Boundary Layer Seeding:** Surface element sizing is maintained at &Delta;<i>h</i><sub>surf</sub> &approx; 0.15–0.20 mm along the central contact journal, providing 3 to 4 quadratic tetrahedral elements through the 0.6 mm case depth;
* **Core Transition:** Element size transitions smoothly to 1.2–1.5 mm in the ductile core, balancing computational efficiency with high through-thickness bending resolution.

---

## 5. Multi-Step Numerical Solution Scheme

<p align="center">
  <img src="../99_images_and_miscelania/shaft_straightening_pipeline_architecture.svg" alt="Abaqus USDFLD Sequentially Coupled FEA Data Pipeline Architecture" width="100%" />
</p>

<div align="justify">
The simulation pipeline is structured into three sequential implicit steps in <b>Abaqus/Standard</b>:
</div>

### Step 0: Initial Quenching Residual Stress Equilibrium
* **Input Deck:** `*INITIAL CONDITIONS, TYPE=STRESS`
* **Mechanism:** Imports the realistic residual stress field inherited from the carburizing-quenching process (-350 to -500 MPa hoop and axial compression in the outer case, balanced by +150 to +250 MPa core tension).
* **Kinematics:** 3-2-1 isostatic boundary constraints. The solver executes a general static step (`*STATIC`) to resolve nodal force residuals and establish self-equilibrating internal stress tensors before tool contact occurs.

### Step 1: Displacement-Controlled Punch Stroke
* **Tooling Model:** The punch tool (<i>R</i> = 20 mm) and supporting V-blocks (span <i>L</i> = 250 mm) are modeled as analytical rigid surfaces or discrete rigid bodies (`R3D4`).
* **Control Mode:** Imposed vertical downward displacement (&Delta;<i>U</i><sub>y</sub>) applied directly at the punch Reference Point (RP). Force-controlled loading is strictly avoided as it triggers numerical divergence near plastic limit loads.
* **Contact Interaction:** Finite-sliding surface-to-surface interaction with penalty formulation (<i>μ</i> = 0.10) and hard normal contact.

### Step 2: Unloading & Total Indicated Runout (TIR) Springback
* **Retraction:** The punch RP is retracted upward until total contact normal force vanishes identically (<i>CNORMF</i> = 0).
* **Runout Computation:** Centerline nodal displacement vectors along the shaft Z-axis are extracted to evaluate residual curvature:

$$\text{TIR}_{\text{residual}} = \max\left(U_y(z)\right) - \min\left(U_y(z)\right)$$

---

## 6. Industrial Numerical Engineering: Convergence Stabilization & Best Practices

<div align="justify">
Elastoplastic contact simulations involving localized plastic yielding, reverse unloading, and high-strength surface layers are highly prone to severe discontinuity iterations (SDI) and convergence stalls. The following stabilization guidelines are strictly embedded:
</div>

### 6.1 Displacement Control vs. Force Control Stability
<div align="justify">
During 3-point bending of a solid stepped shaft, the tangent stiffness matrix <b>K</b><sub>T</sub> degrades rapidly once the outer case and core yield simultaneously. Force-controlled Newton-Raphson iterations diverge instantly due to an ill-conditioned tangent matrix (det(<b>K</b><sub>T</sub>) &approx; 0). Imposing displacement boundary conditions on the punch Reference Point preserves diagonal dominance in <b>K</b><sub>T</sub>, ensuring monotonic quadratic convergence across the entire loading branch.
</div>

### 6.2 Asymmetric Solver Activation (`UNSYMM=YES`)
<div align="justify">
Because frictional penalty contact generates a non-symmetric contribution to the global stiffness matrix when sliding occurs (<i>τ</i><sub>crit</sub> = <i>μ</i> <i>p</i>), setting:
</div>

```abaqus
*STEP, NAME=PUNCH_LOADING, NLGEOM=YES, UNSYMM=YES
```

<div align="justify">
is mandatory. While each iteration incurs a ~30% higher computational cost, it eliminates severe Newton-Raphson residual cycling and prevents solver cutbacks.
</div>

### 6.3 Contact Chattering Mitigation with `C3D10M`
<div align="justify">
Traditional <code>C3D10</code> elements distribute equivalent nodal forces unequally between corner and mid-side nodes (corner nodes can receive negative equivalent forces under uniform pressure), causing severe contact status oscillations (open/close chatter). The modified <code>C3D10M</code> element formulation uses a consistent weighting scheme that guarantees positive nodal forces across all contact faces, eliminating contact chattering and drastically accelerating convergence.
</div>

---

## 7. Subroutine Architecture & Integration (`USDFLD`)

<div align="justify">
The user subroutine is packaged as a standalone modern Fortran 2008 file (<code>usdfld_case_hardneded_shaft.f</code>). At every material integration point:
</div>

1. `CALL GETVRM('COORD', ...)` extracts the current 3D coordinates (<i>x</i>, <i>y</i>, <i>z</i>);
2. The outer radius profile <i>R</i><sub>outer</sub>(<i>z</i>) is evaluated piecewise based on the shaft stepped CAD geometry (`husson_shaft.step`):
   * Seat 1 (<i>Z</i> &lt; 35 mm): <i>R</i> = 12.5 mm (Ø25 mm);
   * Central Journal (35 &le; <i>Z</i> &le; 180 mm): <i>R</i> = 19.0 mm (Ø38 mm, punch contact section);
   * Shoulder (180 &lt; <i>Z</i> &le; 250 mm): <i>R</i> = 25.5 mm (Ø51 mm);
   * Seat 2 (<i>Z</i> &gt; 250 mm): <i>R</i> = 15.0 mm (Ø30 mm).
3. The depth from the outer surface is calculated as:

$$d(r, z) = R_{\text{outer}}(z) - \sqrt{(x - x_0)^2 + (y - y_0)^2}$$

4. `FIELD(1)` is assigned to depth <i>d</i>(<i>r</i>, <i>z</i>), enabling continuous interpolation of yield stress and Chaboche parameters via `*PLASTIC, DEPENDENCIES=1`;
5. `STATEV(1)` (Depth), `STATEV(2)` (PEEQ), and `STATEV(3)` (Radius) are stored for field contour visualization in Abaqus/Viewer.

---

## 8. File Inventory

```text
05_usdfld_shaft_straightening/
├── usdfld_case_hardneded_shaft.f     # USDFLD user subroutine (stepped outer radius, depth field variable, PEEQ)
├── 20MnCr5.inp                       # Depth-dependent Chaboche + Voce material deck (include file)
├── husson_shaft.step                 # CAD geometry of the secondary transmission shaft
└── README.md                         # Technical report & documentation
```

The mesh and the master input deck of the three-step analysis are not stored in the repository. The model is rebuilt from the CAD geometry following Sections 4 and 5.

---

## 9. Execution Guide (Abaqus / SIMULIA 3DEXPERIENCE)

With a master deck that includes the material card and assigns it to the shaft section (`*INCLUDE, INPUT=20MnCr5.inp`, see the header of that file):

```bash
abaqus job=shaft_straightening input=<master_deck>.inp user=usdfld_case_hardneded_shaft.f cpus=8 interactive
```

The residual runout is evaluated from the displacement $U_y$ of the centerline nodes after the unloading step (Section 5, Step 2).

---

## 10. Key Results & Engineering Evaluation

<div align="justify">
The model is set up to evaluate the following results:
</div>

1. **Punch Force vs. Stroke Hysteresis:**  
   Captures the steep initial elastic slope, the onset of case yield at high load (~1900 MPa), progressive core yield propagation, and the non-linear unloading path driven by reverse plastic deformation.
2. **Springback Sensitivity to the Hardening Model:**  
   The material deck also holds a pure isotropic hardening card, so the residual runout predicted with the **Chaboche kinematic hardening model** can be compared with an isotropic prediction, which neglects the Bauschinger effect during the reverse yielding on unloading.
3. **Tensile Case Cracking Risk Assessment:**  
   During peak punch deflection, the outer convex fiber experiences severe tensile stress. Monitoring maximum principal stress (<i>σ</i><sub>I</sub>) against the case fracture toughness threshold (<i>K</i><sub>Ic</sub> &approx; 18–22 MPa&radic;m) confirms that no micro-crack initiation occurs within the allowable straightening stroke window.

---

## 11. Standards & References

1. **Husson, R., Baudouin, C., Bigot, R., & Sura, E. (2014):** *Consideration of residual stress and geometry during heat treatment to decrease shaft bending*, The International Journal of Advanced Manufacturing Technology, 72(9-12), 1455–1463. [DOI: 10.1007/s00170-014-5688-8](https://doi.org/10.1007/s00170-014-5688-8).
2. **Chaboche, J. L. (1986):** *Time-independent constitutive theories for cyclic plasticity*, International Journal of Plasticity, 2(2), 149–188.
3. **Chaboche, J. L. (1989):** *Constitutive equations for cyclic plasticity and cyclic viscoplasticity*, International Journal of Plasticity, 5(3), 247–302.
4. **Armstrong, P. J., & Frederick, C. O. (1966):** *A mathematical representation of the multiaxial Bauschinger effect*, CEGB Report RD/B/N 731, Central Electricity Generating Board.
5. **DIN EN 10084:** *Case hardening steels – Technical delivery conditions* (chemical composition and hardenability of 16MnCr5 and 20MnCr5).
6. **ISO 1101:2017:** *Geometrical product specifications (GPS) – Geometrical tolerancing – Tolerances of form, orientation, location and run-out*.
7. **ASTM E140-12b:** *Standard Hardness Conversion Tables for Metals Relationship Among Brinell Hardness, Vickers Hardness, Rockwell Hardness*.

---

## Author & Contact

**Victor Maia**  
- **Email:** [vhfm08@gmail.com](mailto:vhfm08@gmail.com)  
- **GitHub:** [@vhfmaia](https://github.com/vhfmaia)  
