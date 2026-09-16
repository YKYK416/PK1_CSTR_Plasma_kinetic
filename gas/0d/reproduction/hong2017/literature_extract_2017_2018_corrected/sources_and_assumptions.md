# Sources, corrections, and implementation choices

## Directly transcribed model content

| Content | Paper location | Implementation |
| --- | --- | --- |
| Species | 2018 corrigendum, corrected Table 1 | `SPECIES` block: N2(v1-v8), H2(v1-v3), electronic states, ions, NHx, and four surface species. |
| Electron-impact excitation, ionization, vibrational excitation/de-excitation, N2/H2 dissociation, attachment | 2017 Table 2: W1-W10, W22-W30, W130 | BOLSIG+ declarations/links; fixed analytic fits where printed. |
| V-T and V-V kinetics | 2017 section 2.1.2, equations (1)-(24), Tables 3-4 | Explicit embedded rate calculations and forward/reverse reactions. |
| Neutral, ion, three-body, wall and negative-ion reactions | 2017 Table 2: W31-W140 | Explicit reaction lines and coefficients. |
| Surface chemistry | 2017 section 2.2, equations (30)-(38), Tables 5-6 | Equation (33) for adsorption/E-R, equation (37) for L-H, and equation (38) for N2 dissociative adsorption. |

## Mandatory 2018 corrigendum updates applied

- W4 uses the Carrasco et al. polynomial, not BOLSIG+.
- W81 is excluded because it duplicates W17.
- W97 is `N2(a'1) + H2 -> N2 + 2H`.
- W98 includes H2(v1-v3), H2(B3), H2(B1), H2(C3), H2(A3), and Rydberg H2; its barrier disappears above `0.3 Ev > 16600 K`.
- W29 and W74 are excluded from the corrected model; W6 and W7 remain excluded, while W4 and W8 are retained.
- Corrected H2 molecular constants are used: omega_e=4401 cm^-1, chi_e=0.02758, DeltaE=107.7 K and E10=5988.2 K.
- Corrected L-H activation energies are S15=0.3 eV and S16=0.2 eV.

## Explicit assumptions required for an executable input

1. **Surface choice:** Table 5 contains values for Al2O3 / ND / metal. This input intentionally selects the metal column. To reproduce the other surfaces, replace only the stated sticking probabilities and L-H diffusion barrier (Al2O3 0.5 eV; ND 0.3 eV; metal 0.2 eV).
2. **Vibrational energies for W98:** the papers state that precise Laporta energies were used but do not print them. Values in the file are calculated from the corrected Table 4 H2 constants: v=1 0.515556 eV, v=2 1.001014 eV, v=3 1.456374 eV. They are transparent approximations, not a substitute for the external Laporta table.
3. **Incident energy in equation (38):** the table gives a function of normal incident energy Ez but the paper does not give the coding convention. The reconstruction uses Ez=k_B*Tg in eV (thermal normal-energy scale); this is isolated in the `EZ` assignment.
4. **BOLSIG+ data:** BOLSIG+ reaction labels require a compatible cross-section database. The supplied database came from the pre-existing local ZDPlasKin runtime, so it is recorded as an external dependency rather than as literature-extracted data.
5. **Reactor constants:** V/A=0.007 cm, Lambda=0.1 mm, ST=1e15 cm^-2 and D(300 K)=0.79 cm^2/s follow sections 2.1.5 and 2.2.2. Surface species are converted to equivalent volume density, ST/(V/A), exactly as described in section 2.2.1.

## Deliberate exclusions

The file excludes non-paper surface entropy/BEP closures that appeared in the older local `build_full/kinet.inp`. It also removes two neutral reactions present in that file but absent from Table 2. This keeps the mechanism tied to Hong 2017/2018 rather than blending later modifications.
