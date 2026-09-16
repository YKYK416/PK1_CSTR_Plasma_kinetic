# Pre-submission readiness record — manuscript v1.2

**Date:** 1 September 2026  
**Scope:** Editorial and scientific-readiness audit of `manuscript_submission_en_v1.md` and its editable Word export. This is not a journal-specific compliance certificate.

## Current manuscript claim boundary

The manuscript supports a model-specific finding: within the reconstructed corrected Hong gas-phase N2/H2 CSTR mechanism, the named electronically excited-H2* family is the largest integrated direct-NH source in the accepted continuous-wave map. It also supports a separate source-traceability and topology audit of an Fe(110) DFT–microkinetic branch. It does **not** report a gas–surface flux comparison, a self-consistent plasma prediction, pulse results, energy efficiency, or experimentally calibrated catalyst performance.

## Completed internal checks

| Check | Status | Evidence |
|---|---|---|
| Claim–evidence alignment | Passed | Direct-NH hierarchy, causal disabling, transport scenarios, deterministic rate scenarios, and prospective surface work are separated in the Abstract, Results, Discussion, captions, and Table 3. |
| Numerical-scope wording | Passed | P0 is defined as a terminal-window convergence gate. Its fixed-CW window-energy diagnostic is explicitly not an energy-balance or efficiency test. |
| Cross-model integrity | Passed | The Fe(110) branch is reproduced at its author boundary and treated as a structurally distinct model. The unknown field-transfer factor \(\beta\) is not set to unity or fitted to NH3. |
| Figure inventory | Passed | One graphical abstract and Figures 1–25 are embedded and captioned; all Markdown figure paths resolve locally. |
| Citation closure | Passed | In-text numerical citations and the 16-entry reference list are closed. The Shao–Mesbah bibliographic entry was cross-checked against DOI 10.1021/jacsau.3c00654. |
| Editable Word export | Passed | The v1.2 DOCX has valid XML parts, 26 embedded media items, 26 inline shapes, 3 tables, and the expected new text. A Word-generated review PDF contains 35 pages and the expected Figure 25, data-availability, and reference text. |

## Mandatory items before external submission

1. Select the exact journal and article type; then verify the current official author instructions, manuscript format, data policy, graphical-abstract policy, declarations, and initial-submission requirements.
2. Replace all `TBD` author, affiliation, correspondence, funding, and competing-interest statements with author-approved information.
3. Release the frozen source data, reconstructed inputs, scripts, generated outputs, and provenance ledgers in a citable repository; insert the archive DOI, version, licence, and access statement.
4. Have all authors approve the mechanistic wording, especially the distinction between the Hong gas result and the Fe(110) source audit.
5. Perform an independent scientific/code review of the source reconstruction and the reported numerical outputs before submission.

## Highest-impact work for a stronger-journal submission

| Priority | Required evidence | Why it matters |
|---|---|---|
| 1 | A boundary-harmonized hybrid calculation or experiment with independently calibrated \(F_s=\beta E_{\mathrm{bulk}}\), catalyst inventory, area-to-volume ratio, site density, and inlet/outlet operator | It would convert the present gas–surface comparison guardrail into a real competing-pathway test. |
| 2 | Experimental validation at the three pre-specified discrimination states, using the blocked factorial protocol in Figure 25 | It would test whether the model’s upstream direct-NH gate and downstream turnover separation survive reactor reality. |
| 3 | A numerically accepted, self-consistent pulsed model with phase-resolved electron kinetics | It would address the principal operating-mode limitation without back-filling claims from the current rejected pulse branch. |
| 4 | Formal correlated uncertainty analysis of key electron-impact and NH-entry coefficients | It would replace the present transparent deterministic scenarios with probabilistic robustness statements. |

## Submission decision

The manuscript is internally coherent and ready for author review as a detailed model-and-mechanism manuscript. It is **not yet ready for external submission** because author/declaration information, a public archival release, a selected venue, and independent cross-boundary or experimental validation remain outstanding. Those are real evidence and compliance gaps, not formatting issues.
