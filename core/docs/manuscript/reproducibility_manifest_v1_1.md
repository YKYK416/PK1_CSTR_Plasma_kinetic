# Reproducibility manifest — manuscript version 1.1

## Purpose and scope

This manifest identifies the local inputs and primary analysis artefacts supporting
the continuous-wave CSTR results reported in Figures 1–18 of the main manuscript. It is an internal
pre-archive manifest, not a substitute for a citable public release. The reported
manuscript results comprise only calculations that satisfy the stated terminal-window
P0 numerical-convergence criterion.

## Mechanism provenance

The calculation uses the locally documented, literature-traceable reconstruction
of the corrected Hong N2/H2 plasma-ammonia mechanism. It must not be represented
as the unique original-author archival input deck. Reconstruction assumptions and
source tracing are documented in the sources-and-assumptions and
extraction-validation reports under the corrected reconstruction folder.

| Artefact | Relative path | SHA-256 |
|---|---|---|
| Reconstructed kinetic input | ../Reproduction/2017Hong/literature_extract_2017_2018_corrected/kinet.inp | D9B0DA025E1EEEEF42FA58963DADCC923366D89AD0621A8AB4516DA966794F38 |
| Electron-impact database | ../Reproduction/2017Hong/literature_extract_2017_2018_corrected/bolsigdb.dat | DB210584EEE6EF78A8CD420F5385F80AB359F8A9547A0E05B8ACBF59EE56182E |

## Primary accepted continuous-wave map

| Item | Definition |
|---|---|
| Reactor boundary | 0D continuously fed and exhausted CSTR; 300 K; residence time 10 ms |
| Electron boundary | Electron density imposed at 1.17 × 10^8 cm^-3 |
| Map | x_N2 = 0.1–0.9 in 0.1 steps; E/N = 20–240 Td in 20 Td steps |
| Accepted cases | 108 continuous-wave CSTR cases |
| Acceptance gate | State and NH3 residuals <= 10^-3 across three consecutive terminal numerical-analysis windows; the stored window-energy diagnostic is retained for driver consistency and is not an energy-balance test |
| Explicit exclusions | Self-consistent electron density/power, pulse-frequency and duty-cycle claims, surface chemistry, experimental calibration |

## Data, figures, and scripts

| Artefact | Relative path | SHA-256 |
|---|---|---|
| Dense direct-NH map data | ../analysis/figures_20260829/Figure_3_dense_EN_N2_source_data.csv | 699B7015D11B45525A80613547F22E5ED2DE0EB8366CA90F539B8C603A8D4C54 |
| Local reaction-control analysis | ../工具脚本/hong_reaction_control.py | 6AE743BF0E9D4E54BCB5BE6372C01E3E7C94A6017E8DA34DBF053B23AB2A013C |
| Residence-time robustness script | ../工具脚本/hong_tau_robustness.py | 0E0394771F2806C06AB0744B10973A6DDE869DB60D592438661977F2AACB3786 |
| Residence-time case table | ../analysis/p8_tau_robustness_20260831_r3/P8_tau_cases.csv | 3C4AE7C99E674D22F25FDB6E4C437AC00308FD9B4F70D5EBC23681A6A792FD21 |
| Figure 17 | ../analysis/p8_tau_robustness_20260831_r3/figures/Figure_17_residence_time_robustness.png | 418D8941585F313C5CA21070AF37DE0CEBB46522E8ED7B5F421746CB1F354117 |
| H2* scenario-envelope script | ../工具脚本/hong_h2star_envelope.py | 8A25D8B92FA2CC3170B6445CB2059CCB1FA67B468AC027BC47BD18D2978FD235 |
| H2* scenario case table | ../analysis/p9_h2star_scenario_envelope_20260831_r3/P9_scenario_cases.csv | ED929AE66BB6F7432FA7269475D52C5B5BE7135F9D9D6676D47B31B5A73FA9C9 |
| Figure 18 | ../analysis/p9_h2star_scenario_envelope_20260831_r3/figures/Figure_18_H2star_scenario_envelope.png | 695AB98E2406FB44003FE460026C3D3058DDA641CC3896D26BB34E3DC5F39BE0 |
| Word exporter | ../工具脚本/export_manuscript_word.py | 47116C53FC1ED19E1F59FD7B67D4915B103048ADF66EAEFC04F4FB20D7EBC3B4 |
| Main figures 1–16 and graphical abstract | ../analysis/figures_20260829/ | Embedded in manuscript |

## Analysis interpretation safeguards

- Direct-NH ranks use terminal-window integrated forward ROP, not an instantaneous
  rate snapshot.
- The named H2* family is a mechanism-level aggregation; it is not a new
  spectroscopic assignment.
- Local control coefficients use symmetric 1/1.10 and 1.10 perturbations and are
  local diagnostics, not a global uncertainty quantification.
- Diagnostic family disabling supports causal relevance at the stated reference
  condition; it does not establish an invariant global rate-limiting step.
- Pulse calculations that fail P0 are excluded from the reported result set.
- Figure 17 contains 20 P0-accepted CW-CSTR cases spanning 0.5–30 ms at four
  contrasting states; it is a targeted transport screen, not a full response surface.
- Figure 18 contains 14 P0-accepted independent 1/3–3× scenarios for three
  H2*-relevant reaction groups at two reference states; it is a deterministic
  sensitivity envelope, not statistical uncertainty quantification.

## Release checklist before submission

1. Freeze a release containing input files, exact code revision, environment
   specification, source data, generated figures, and this manifest.
2. Upload the release to a public repository that assigns a permanent DOI.
3. Replace the manuscript data-availability placeholder with the DOI, version,
   access date, and license.
4. Re-run the figure and Word-generation workflow from the frozen release and
   record the verification output.
