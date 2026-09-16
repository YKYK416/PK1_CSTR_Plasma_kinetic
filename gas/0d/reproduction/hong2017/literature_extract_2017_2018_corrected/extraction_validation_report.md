# Extraction validation report

## Completed checks

- ZDPlasKin 2.0a preprocessor successfully read the input and wrote `zdplaskin_m.F90`: 48 loaded species, 504 reaction records, 41 BOLSIG+ references and 404 embedded Fortran lines.
- Corrected species list includes all entries in corrigendum Table 1, including H2 Rydberg.
- All corrigendum deltas have been searched in the input: W81 absent; corrected W97 present; W98 covers vibrational, electronic and Rydberg H2; W29/W74 absent; S15/S16 corrected energies present.
- The gas-phase portion retains the explicit vibration-rate implementation used for W11-W18 and reaction groups W19-W140.
- Surface reactions implement the 18 Table 5 reaction families under the metal-surface choice, with site-density conversion made explicit.

## Remaining non-uniqueness

The preprocessor emitted pre-existing long embedded-Fortran-line warnings and two duplicate-reaction warnings in the BOLSIG excitation block, but no syntax error. The two articles do not contain the original authors' literal code, exact BOLSIG+ cross-section file/version, Laporta energy table, initial densities, or a unique numerical prescription for the incident normal energy in the N2 sticking fit. These inputs cannot be claimed as extracted solely from the PDFs and are documented in `sources_and_assumptions.md`.

## Recommended runtime check

Run the ZDPlasKin preprocessor in a disposable working directory with `kinet_hong_2017_2018_corrected.inp` renamed to `kinet.inp` and `bolsigdb.dat` alongside it. Then verify every BOLSIG+ process resolves against the intended database and set the initial `Surf` population to ST/(V/A) before a physical simulation.
