# 06. Hydrogen Transport for Weld Cold Cracking (`UMATHT`)

[![Status: In Work](https://img.shields.io/badge/Status-In%20Work%20(Abaqus%20run%20pending)-orange.svg)](#)
[![FEA: Abaqus/Standard](https://img.shields.io/badge/FEA-Abaqus%2FStandard-blue.svg)](https://www.3ds.com/)
[![Subroutine: UMATHT](https://img.shields.io/badge/Subroutine-UMATHT-005691.svg)](#)
[![Language: Fortran](https://img.shields.io/badge/Language-Fortran%202008%20%2F%20F90-734f96.svg)](https://fortran-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](../LICENSE)

---

## Executive Summary

> [!WARNING]
> **DEVELOPMENT STATUS: IN WORK**  
> The subroutine and its 1D verification deck are written; the formulation has been checked against the analytical solutions with an independent 1D finite element replica of the same constitutive update (section 4). The Abaqus run of the deck is pending, and the coupling to the welding model of project 02 is the next stage (section 5).

Hydrogen-assisted cold cracking, also called delayed cracking, is a typical failure of laser welds between case-hardening steels and cast irons. During cooling, the austenite of the heat-affected zone (HAZ) transforms to martensite. Hydrogen is far more soluble in austenite but diffuses there several orders of magnitude more slowly, so it is left supersaturated in a hard, highly stressed microstructure. It then diffuses towards the hydrostatic tension and gathers in traps, and a crack can appear hours to days after welding.

This project models the transport part of that chain with a <b>User Material for Heat Transfer (<code>UMATHT</code>)</b>. The nodal temperature degree of freedom carries the lattice hydrogen concentration. This is the heat transfer analogy of Díaz <i>et al.</i> (2016), and it keeps the same <code>DC3D8</code> meshes, solver and post-processing as the welding thermal analysis of project 02.

---

## 1. Physical Model

Mass balance of lattice ($C_L$) and trapped ($C_T$) hydrogen, with Oriani local equilibrium between the two populations:

$$\frac{\partial (C_L + C_T)}{\partial t} + \nabla \cdot \mathbf{J} = 0$$

$$\mathbf{J} = -D_L \nabla C_L + \frac{D_L C_L \bar{V}_H}{R T} \nabla \sigma_h$$

$$\frac{\theta_T}{1-\theta_T} = K_T \frac{\theta_L}{1-\theta_L}, \qquad K_T = \exp\left(\frac{W_B}{R T}\right), \qquad \theta_L = \frac{C_L}{N_L}, \quad C_T = N_T\, \theta_T$$

Where:
* $D_L$: lattice diffusivity, Arrhenius per phase and mixed with the austenite fraction $f_\gamma$: $D_L = (1-f_\gamma) D_{0,\alpha} e^{-Q_\alpha/RT} + f_\gamma D_{0,\gamma} e^{-Q_\gamma/RT}$;
* $\bar{V}_H$: partial molar volume of hydrogen; the second term of $\mathbf{J}$ drives hydrogen up the gradient of the hydrostatic stress $\sigma_h$ (Sofronis & McMeeking, 1989);
* $W_B$, $N_T$: trap binding energy and trap density; $N_T$ can include the dislocation traps created by plastic strain (Kumnick & Johnson, 1980): $\log_{10} N_T\,[\text{sites/m}^3] = 23.26 - 2.33\, e^{-5.5\,\varepsilon_p}$.

---

## 2. Heat Transfer Analogy in the `UMATHT`

| Heat transfer (Abaqus) | Hydrogen transport | Implementation |
| :--- | :--- | :--- |
| Temperature `NT11` | Lattice concentration $C_L$ | Primary unknown |
| Internal energy `U` (density = 1) | Total hydrogen $C_L + C_T$ | Increment from the exact trapped amounts at the start and the end of the increment |
| `DUDT` | $1 + \partial C_T / \partial C_L = D_L / D_{\text{eff}}$ | $\dfrac{\partial C_T}{\partial C_L} = \dfrac{N_T}{N_L}\dfrac{K_T}{(1-\theta_L+K_T\theta_L)^2}$ |
| `FLUX` ($-k\nabla T$) | $\mathbf{J}$ | Fick + stress drift |
| `DFDG` | $-D_L\,\mathbf{I}$ | Exact |
| `DFDT` | $D_L \bar{V}_H \nabla\sigma_h / (R T)$ | Exact |

Because the increment of <code>U</code> comes from the exact Oriani equilibrium rather than from <code>DUDT·ΔC</code>, the hydrogen balance holds at any trap occupancy. Traps created by plastic strain, or emptied by heating ($K_T$ falls with $T$), act as sinks and sources of lattice hydrogen without any extra term.

### Inputs

| `PROPS` | Symbol | Unit | | Field variable | Quantity | Unit |
| :---: | :--- | :--- | --- | :---: | :--- | :--- |
| 1 | $D_{0,\alpha}$ (ferrite / martensite) | mm²/s | | 1 | Temperature | °C |
| 2 | $Q_\alpha$ | J/mol | | 2 | Austenite fraction $f_\gamma$ | – |
| 3 | $D_{0,\gamma}$ (austenite) | mm²/s | | 3–5 | $\nabla\sigma_h$ (x, y, z) | MPa/mm |
| 4 | $Q_\gamma$ | J/mol | | 6 | Equivalent plastic strain | – |
| 5 | $N_L$ | mol/m³ | | | | |
| 6 | $\bar{V}_H$ | mm³/mol | | | | |
| 7 | $W_B$ | J/mol | | | | |
| 8 | $N_T$ (constant part) | mol/m³ | | | | |
| 9 | Kumnick-Johnson traps (1 = on) | – | | | | |

State variables (`*DEPVAR 4`): $C_L$, $C_T$, $D_L$, $N_T$. Concentrations can be in any one unit (mol/m³ here; 1 wppm ≈ 7.8 mol/m³ in steel), as long as $N_L$ and $N_T$ use the same one.

---

## 3. Reference Values for a Case-Hardening Steel

| Parameter | Value | Source |
| :--- | :--- | :--- |
| $D_{0,\alpha}$, $Q_\alpha$ (bcc iron) | 7.23 × 10⁻² mm²/s, 5.69 kJ/mol | Kiuchi & McLellan (1983); $D_L \approx 7 \times 10^{-3}$ mm²/s at 20 °C |
| $N_L$ | 5.1 × 10²⁹ sites/m³ = 8.47 × 10⁵ mol/m³ | Sofronis & McMeeking (1989) |
| $\bar{V}_H$ | 2.0 × 10⁻⁶ m³/mol = 2000 mm³/mol | Sofronis & McMeeking (1989) |
| $W_B$ (dislocations) | 60 kJ/mol | Sofronis & McMeeking (1989) |
| $D_{0,\gamma}$, $Q_\gamma$ (austenite) | 1.0 mm²/s, 50 kJ/mol | **Order-of-magnitude estimate by the author** ($D_\gamma \approx 10^{-9}$ mm²/s at 20 °C), to be replaced by data for the grade |

Martensite has a lower effective diffusivity than bcc iron because of its lath boundaries and dislocation density; in this model that enters through the trap terms, not through $D_L$.

---

## 4. Verification Deck and Expected Results

<code>hydrogen_transport_verification.inp</code> holds three independent bars of 40 <code>DC3D8</code> elements with a 1 × 1 mm section, run in one transient step of 1000 s with fixed increments of 0.1 s. The expected values below are the analytical solutions; the "1D replica" column is an independent finite element implementation of the same <code>UMATHT</code> update (lumped capacity, 2-point Gauss flux, backward Euler) and shows the discretisation error that the Abaqus run should reproduce.

### 4.1 Bar A: Permeation Transient (Lattice Only)

$C(0) = 1$, $C(L) = 0$, $L = 1$ mm, $D_L = 0.01$ mm²/s. Exit flux (Devanathan-Stachurski permeation):

$$\frac{J(t)}{J_{ss}} = 1 + 2\sum_{n=1}^{\infty} (-1)^n \exp\left(-\frac{n^2 \pi^2 D t}{L^2}\right), \qquad J_{ss} = \frac{D C_0}{L} = 0.01, \qquad t_{\text{lag}} = \frac{L^2}{6D} = 16.67\text{ s}$$

| *t* [s] | 5 | 10 | 16.67 | 25 | 50 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Analytical $J/J_{ss}$ | 0.0340 | 0.2929 | 0.6167 | 0.8305 | 0.9856 |
| 1D replica | 0.0353 | 0.2873 | 0.612 | 0.8268 | 0.9851 |

### 4.2 Bar B: Permeation with Oriani Traps

$W_B$ = 10 kJ/mol ($K_T$ = 60.507 at 20 °C) and $N_T$ = 125 967 mol/m³, so that $N_T K_T / N_L = 9$. At low occupancy ($\theta_L K_T < 10^{-4}$) the effective diffusivity is $D_{\text{eff}} = D_L / (1 + N_T K_T / N_L) = D_L / 10$: same steady flux as bar A, ten times slower transient.

| *t* [s] | 50 | 100 | 166.7 | 250 | 500 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Analytical ($D_{\text{eff}}$) | 0.0340 | 0.2929 | 0.6167 | 0.8305 | 0.9856 |
| 1D replica | 0.0349 | 0.2931 | 0.6165 | 0.8302 | 0.9856 |

### 4.3 Bar C: Stress-Driven Equilibrium

Closed bar ($\mathbf{J} \cdot \mathbf{n} = 0$ at both ends), $L = 2$ mm, initial $C = 1$, uniform $\nabla\sigma_h$ = 500 MPa/mm (field variable 3). At equilibrium the flux vanishes everywhere and the hydrogen redistributes towards the tensile end, conserving its total amount:

$$C(x) = \bar{C}\, \frac{k L\, e^{k x}}{e^{k L} - 1}, \qquad k = \frac{\bar{V}_H\, \partial\sigma_h/\partial x}{R T} = 0.41028\text{ mm}^{-1}$$

| *x* [mm] | 0 | 0.5 | 1.0 | 1.5 | 2.0 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Analytical *C* | 0.64521 | 0.79212 | 0.97249 | 1.19392 | 1.46576 |
| 1D replica (*t* = 1000 s) | 0.64518 | 0.79209 | 0.97245 | 1.19388 | 1.46573 |

The total amount stays at 2.000000 (the integral of *C* over the length).

---

## 5. Coupling to the Welding Simulation (Roadmap)

The intended analysis chain reuses the thermal and metallurgical results of the welded crown gear and differential carrier of projects 01 and 02:

1. **Thermal + metallurgical** (project 02, `DFLUX`): temperature history and austenite / martensite fractions.
2. **Mechanical** (sequential): residual stresses with transformation strains → $\sigma_h$ and the equivalent plastic strain.
3. **Hydrogen transport** (this `UMATHT`, same mesh): diffusion during cooling, then a long room-temperature step (hours to days) for the delayed accumulation in the HAZ. The initial hydrogen goes into the fusion zone (diffusible hydrogen of the process, e.g. a few ml/100 g).
4. **Fracture**: the cohesive law of [project 04](../04_umat_cohesive_welding/README.md), extended with a hydrogen-dependent strength, $\sigma_c(\theta) = \sigma_c(0)(1 - 1.0467\,\theta + 0.1687\,\theta^2)$ (Serebrinsky, Carter & Ortiz, 2004), placed along the HAZ boundary.

Open points before step 3 can run on the weld model:
- The nodal temperatures of the thermal ODB can be read as field variable 1 (same mesh); the austenite fraction, $\nabla\sigma_h$ and $\varepsilon_p$ are integration-point results of other analyses and need a mapping (`UFIELD`, or a coupled temperature-displacement analysis that computes $\nabla\sigma_h$ in the same run).
- **Solubility jump between phases:** austenite dissolves much more hydrogen than ferrite / martensite. With a single $N_L$ and $C_L$ continuous, this jump (and the supersaturation it causes on transformation) is not represented yet; the usual remedy is to solve for a normalised concentration $\phi = C_L / S(T, f_\gamma)$, as in the Abaqus `*MASS DIFFUSION` procedure.

---

## 6. File Inventory

```text
06_umatht_hydrogen_diffusion/
├── umatht_hydrogen_transport.f            # UMATHT: Oriani trapping, stress drift, phase / temperature diffusivity
├── hydrogen_transport_verification.inp    # 3-bar 1D verification deck (permeation, traps, stress drift)
└── README.md                              # This document
```

## 7. How to Run the Verification Deck

```bash
abaqus job=hydrogen_transport_verification user=umatht_hydrogen_transport.f interactive
```

Check against section 4: the summed <code>RFL11</code> of <code>A_OUT</code> and <code>B_OUT</code> in the <code>.dat</code> file (exit flux; its magnitude is <i>J</i> for the 1 mm² section, steady value 0.01), and <code>NT11</code> of <code>C_PROFILE</code> at the end of the step.

---

## 8. References

1. Oriani, R.A. (1970). The diffusion and trapping of hydrogen in steel. *Acta Metallurgica*, 18, 147–157.
2. Sofronis, P., McMeeking, R.M. (1989). Numerical analysis of hydrogen transport near a blunting crack tip. *Journal of the Mechanics and Physics of Solids*, 37, 317–350.
3. Kumnick, A.J., Johnson, H.H. (1980). Deep trapping states for hydrogen in deformed iron. *Acta Metallurgica*, 28, 33–39.
4. Kiuchi, K., McLellan, R.B. (1983). The solubility and diffusivity of hydrogen in well-annealed and deformed iron. *Acta Metallurgica*, 31, 961–984.
5. Díaz, A., Alegre, J.M., Cuesta, I.I. (2016). Coupled hydrogen diffusion simulation using a heat transfer analogy. *International Journal of Mechanical Sciences*, 115–116, 360–369.
6. Serebrinsky, S., Carter, E.A., Ortiz, M. (2004). A quantum-mechanically informed continuum model of hydrogen embrittlement. *Journal of the Mechanics and Physics of Solids*, 52, 2403–2430.
