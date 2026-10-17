# Metallurgical Phase Transformation Cards: One-Element Verification

### How the material cards of the laser weld model (project 02) behave in the Abaqus metallurgical phase transformation framework, one property at a time

The crown gear (20MnCr5), the differential carrier (EN-GJS-600-10) and the filler wire (Ni 6082) of [project 02](../02_dflux_laser_welding/) are modelled with the metallurgical phase transformation framework of Abaqus/Standard (`ABQ_PHASE_TRANS` parameter and property tables, Abaqus Materials Guide, "Modeling Metallurgical Phase Transformations"). This folder checks each part of those cards on one or two elements, where the response can be compared with the card data and with hand calculation, before it is trusted in a model of 2.8 million equations.

Status: the first set of test decks is written; results, the property plots against temperature, the model TTT / CCT diagrams and the sources of every table will be added as the tests are run.

---

## Tests

| Deck | Element | What it checks |
| :--- | :--- | :--- |
| `tests/PT01_cct_20MnCr5.inp`, `tests/PT01_cct_GJS600_10.inp` | DC3D8, temperature imposed | phase fractions after austenitising at 1000 C and cooling at 0.3, 1, 3, 10, 30, 100, 300 and 1000 C/s (one cycle per rate): the model CCT diagram |
| `tests/PT02_retention_rule.inp` | 3 x DC3D8, synthetic card | the retained-fraction table of the martensitic transformation: 80 % ferrite formed at 600 C, then a quench, with and without the rule (0, 0), (0.25, 0.25), (0, 1); tells whether its starting proportion is the absolute austenite fraction |
| `tests/PT03_mixture_elastic.inp` | C3D8, 20MnCr5 | Young's modulus of 50 % austenite + 50 % martensite (quench to 316 C, Koistinen-Marburger) against the rule of mixtures expected from the interpolation in the volume fraction ratios |
| `tests/PT04_latent_heat.inp` | 2 x DC3D8, synthetic card | latent heat of fusion released linearly against `*LATENT HEAT, SMOOTH` (adiabatic element, constant body heat flux) |
| `tests/PT05_dilatometry_20MnCr5.inp`, `..._GJS600_10.inp`, `..._Inconel82.inp` | C3D8, free | simulated dilatometer: total and thermal strain against temperature on heating and on cooling at 10 and 300 C/s, with the transformation strains |

Temperatures are imposed on every node (heat transfer: `*BOUNDARY` on degree of freedom 11; stress: `*TEMPERATURE`), with fixed increments of about 2 C. The synthetic cards of PT02 and PT04 are written for the test (constant properties); the other decks include the cards of project 02, which are copied next to the deck before running:

```bash
# from 02_dflux_laser_welding/THERMAL_RLS: Material_phase_trans_table_types.inp, Material_20MnCr5.inp,
# Material_GJS600_10.inp, Material_Inconel82_RLS.inp
abaqus job=PT01_cct_20MnCr5 input=PT01_cct_20MnCr5.inp interactive
# results: *EL PRINT tables in PT01_cct_20MnCr5.dat (centroid values every increment); the ODB keeps the same
# history output as a backup
```

Units: mm, s, tonne, N, MPa, mJ, C. Abaqus/Standard 2025x.

## Reference

Abaqus Materials Guide, 2025 edition (Dassault Systemes): Modeling Metallurgical Phase Transformations; Parameter and Property Table Type Reference; Latent Heat; Specifying Material Data as Functions of Temperature and Independent Field Variables. The guide and its verification inputs are cited, not reproduced; the decks and cards in this folder are written for this project.
