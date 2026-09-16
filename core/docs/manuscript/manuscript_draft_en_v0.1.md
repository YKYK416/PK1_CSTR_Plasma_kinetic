---
status: Internal detailed English draft; not submitted and not journal-formatted
article_type: Computational plasma-kinetics mechanism study
target_style: Plasma Sources Science and Technology / Journal of Physics D
citation_style: Numbered citations represented by BibTeX keys during drafting
---

# Electronically Excited H2 Sustains the Dominant NH Formation Route across a Dense N2–E/N Map in a Plasma CSTR

**[Author name(s) to be inserted]**

**[Affiliation(s) to be inserted]**

## Abstract

Plasma-assisted ammonia synthesis contains coupled electron, vibrational, electronic-excitation, radical, ion, and surface-chemistry processes. Consequently, an NH3 concentration alone cannot identify the precursor route that produces NH, a central intermediate in nitrogen hydrogenation. Here, we implement the Hong N2/H2 nonequilibrium plasma mechanism, including its published corrigendum, in a zero-dimensional continuous stirred-tank reactor (CSTR) framework. The CSTR boundary is coupled to a strict period-steady-state gate: the maximum relative change of significant species, the cycle-to-cycle NH3 residual, and the cycle-to-cycle energy residual must each remain at or below 10^-3 for three consecutive cycles. NH source fractions are then obtained by integrating reaction-of-progress data only over the final three accepted cycles. We evaluate 108 independent continuous-wave cases spanning N2 inlet fractions of 0.1–0.9 and reduced electric fields of 20–240 Td at 300 K, a 10 ms residence time, and a fixed electron density of 1.17×10^8 cm^-3. The four explicitly represented electronically excited H2 channels remain the largest NH source in every accepted case, contributing 84.892–99.664% of the integrated NH production. Thus, this model and boundary set do not predict a crossover of the dominant NH formation route. Instead, high-field and nitrogen-rich conditions redistribute secondary flux through Rydberg-H2*, N(2P)+H2, and three-body association channels. The steady NH3 concentration exhibits a distinct intermediate-field window: its maximum occurs at 140 Td for N2 fractions from 0.1 to 0.8 and at 120 Td for N2=0.9. Disabling the named electronic-H2* channels lowers the baseline NH3 concentration by 83.31%, whereas diagnostic Rydberg-loss perturbations do not reverse their dominance. A symmetric local kinetic-control panel identifies the downstream rate controls: NH+H2+M to NH3+M dominates the local productivity response at the optimum (S=1.035), whereas H+NH2+M to NH3+M (S=0.974) and H+ + NH3 to NH3+ + H (S=-0.856) oppose one another at high field. These results separate NH-source topology from net NH3 accumulation and provide a numerically auditable steady-state baseline for future pulse-afterglow studies.

**Keywords:** plasma-assisted ammonia synthesis; NH radical; reaction-of-progress analysis; electronically excited hydrogen; CSTR; zero-dimensional kinetics; reduced electric field.

## 1. Introduction

Ammonia synthesis in nonthermal plasmas is governed by a reaction environment that differs fundamentally from a thermal reactor. Electron-impact processes create atoms, ions, vibrationally excited molecules, and electronically excited molecules, while their subsequent gas-phase and surface interactions determine whether nitrogen-containing intermediates are hydrogenated, recombined, or destroyed. Mechanistic interpretation is therefore not equivalent to reporting an outlet NH3 concentration. A chemically meaningful analysis must identify the precursor routes that create NH and determine whether their importance changes when the feed composition or electrical excitation is varied. Reviews of plasma-driven ammonia synthesis have emphasized that radical-driven, plasma-enhanced catalytic, and surface-enhanced pathways can coexist and that their relative relevance depends on the discharge and catalytic boundary conditions [Rouwenhorst2020].

Detailed plasma-kinetics models provide a route to state-resolved mechanism testing. Hong *et al.* developed a kinetic model for atmospheric-pressure N2–H2 plasma catalysis that includes electron kinetics, vibrationally excited N2 and H2, electronically excited states, and prescribed surface reactions [Hong2017]. The subsequent corrigendum corrected several reaction and rate-data entries and is part of the mechanism provenance used here [Hong2018Corrigendum]. In that modelling framework, NH is an early precursor to ammonia formation. Later zero-dimensional DBD modelling showed that vibrational kinetics can affect ammonia predictions when the temporal and spatial character of the discharge is represented through microdischarge and uniform-plasma components [VantVeer2020], while kinetic analysis in plasma-enhanced catalysis has identified conditions in which vibrationally excited N2 can alter surface activation [Rouwenhorst2019]. These studies motivate detailed pathway analysis, but they do not by themselves establish whether a change in N2 fraction or reduced field should cause a crossover of the *dominant* NH-producing reaction family in a common, validated steady-state setting.

That distinction is important for two reasons. First, a pathway fraction measured during a transient or in a closed batch reactor can be influenced by the accumulation of slow products and ions rather than by a reproducible reactor state. Second, a dominant source of NH need not be the condition that maximizes the final NH3 inventory. Net NH3 accumulation also depends on subsequent NHx hydrogenation, destruction reactions, residence time, and the imposed electron boundary. A study that equates the largest instantaneous NH source with the largest NH3 yield risks conflating chemical topology with reactor performance.

This work tests a falsifiable question: within the Hong mechanism under a controlled CSTR/continuous-wave (CW) boundary, does the dominant NH formation route switch across a dense N2–E/N plane? Rather than presuming a pathway switch, we first establish a periodic steady-state gate and then integrate individual reaction-of-progress (ROP) records only after the gate has been passed. The analysis covers 108 cases with N2 inlet fractions from 0.1 to 0.9 and E/N from 20 to 240 Td. We then use reaction disabling and a parameterized Rydberg-loss perturbation to distinguish a high pathway fraction from a causally important pathway.

The study makes three contributions. (i) It introduces a CSTR boundary and a transparent P0 convergence criterion for pathway comparisons in this local implementation of the Hong mechanism. (ii) It supplies a dense, reproducible N2–E/N map of terminal steady-state NH-source fractions and demonstrates that the named electronically excited H2 channels remain dominant throughout the accepted grid. (iii) It separates this source-topology result from the NH3 product window, which peaks at intermediate reduced field and is not predicted by the electronic-H2* NH-source fraction alone. The conclusions are deliberately restricted to the stated fixed-electron-density, 0D CSTR/CW boundary. Unvalidated pulse-afterglow results are not used.

## 2. Methods

### 2.1. Mechanism provenance and controlled model boundary

The calculation uses the local ZDPlasKin implementation of the Hong N2–H2 nonequilibrium plasma mechanism and incorporates the published 2018 corrigendum [Hong2017, Hong2018Corrigendum]. The generated kinetic module contains 48 state species and 469 elementary reactions. Its state set includes ground-state, vibrationally excited, electronically excited, radical, ionic, and prescribed surface species. The pathway conclusions in this paper concern **gas-phase NH source reactions**; prescribed surface species remain part of the reaction environment but are not assigned the CSTR feed/outflow term described below.

The total neutral number density is specified as

`n_tot = 2.446 × 10^19 (300 / Tg) cm^-3`,

where `Tg` is in kelvin. Unless explicitly stated otherwise, the controlled baseline uses `Tg = 300 K`, a continuous applied field, a residence time `tau = 1.0 × 10^-2 s`, and a fixed electron density `ne = 1.17 × 10^8 cm^-3`. The inlet mixture is N2/H2, with the N2 fraction varied as specified in Section 2.4. A fixed electron density is used here as a controlled electron-supply boundary. It is not a charge-self-consistent discharge prediction and must not be interpreted as a measurement of the electron density in a specific device.

### 2.2. CSTR formulation

The original closed zero-dimensional treatment was extended with a CSTR term inside the solver right-hand side. For each heavy gas-phase species *i*,

`dn_i/dt = omega_i + (n_feed,i - n_i)/tau`,

where `omega_i` is the chemical source term, `n_feed,i` is the feed number density, `n_i` is the instantaneous number density, and `tau` is the residence time. The term is applied to heavy gas species only. The electron, bookkeeping species M and S, SURF-labelled species, and surface species are excluded. This exclusion prevents the imposed material-flow term from acting as an artificial electron or surface-site source.

The CSTR model is an average residence-time boundary, not a spatially resolved plasma reactor model. It cannot represent filamentary microdischarge structure, sheath physics, axial gradients, local heating, or physical wall transport. Its role is narrower: it prevents a closed-batch inventory from being mistaken for a periodic flow-reactor state and makes pathway fractions comparable across operating points with the same prescribed residence time.

### 2.3. Period-steady-state acceptance criterion

Each case is propagated for up to 180 field cycles. A case is accepted only when all three period-to-period quantities are not larger than `1.0 × 10^-3` for three consecutive cycles:

1. the largest relative endpoint change among significant species (`rel_state_max`);
2. the NH3 cycle residual normalized by the NH3 inventory (`rel_dNH3_cycle`); and
3. the relative cycle-energy change (`rel_energy_cycle`).

The NH3 residual is evaluated using the larger of the start-of-cycle NH3 density, the end-of-cycle NH3 density, and a monitor floor. This normalization avoids declaring convergence only because a small absolute difference is compared with an arbitrarily large or small scale. Cases that do not pass this P0 gate are rejected before any pathway integration. The gate also requires a terminal stable streak of three cycles.

The importance of this distinction was assessed at the reference condition N2=0.3333, `Tg=300 K`, `E/N=120 Td`, 1 kHz CW, and fixed `ne`. Table 1 compares the closed-batch and CSTR outcomes. The closed calculation remains outside the P0 criterion after 60 cycles, whereas the CSTR reaches the acceptance condition at cycle 55. The comparison establishes the numerical boundary for the rest of the paper; it is not a claim that a 0D CSTR replaces a spatially resolved reactor model.

| Boundary | Maximum cycles | Status | Final NH3 (cm^-3) | `rel_state_max` | `rel_dNH3_cycle` |
|---|---:|---|---:|---:|---:|
| Closed batch | 60 | Not accepted | 8.292×10^14 | 1.73×10^-2 | 1.09×10^-2 |
| CSTR, `tau=10 ms` | 120 | Accepted at cycle 55 | 1.534×10^14 | 8.39×10^-4 | 6.58×10^-4 |

**Table 1.** Comparison demonstrating why all mechanism conclusions in this manuscript are based on accepted CSTR states rather than closed-batch endpoint values.

### 2.4. Dense N2–E/N scan

The principal scan contains nine N2 inlet fractions, `0.1, 0.2, ..., 0.9`, and twelve reduced fields, `20, 40, ..., 240 Td`, for 108 cases in total. Every point is executed as an isolated case with its own input record, output directory, P0 history, and reaction-rate output. The calculation is CW (`duty=1`), at 1 kHz, 300 K, the fixed electron-density boundary stated above, and `tau=10 ms`.

One hundred and six points passed on the first attempt. Two low-field points, 20 Td at N2=0.1 and N2=0.3, completed through the registered numerical retry path. Their final `rel_state_max` values were `6.40×10^-4` and `7.83×10^-4`, respectively, and both had `rel_dNH3_cycle < 6×10^-9`. They therefore satisfy the same P0 gate as the first-pass cases. Their retry status is retained in the released data table and marked by crosses in the dense pathway maps; it is not hidden by averaging or replacement.

### 2.5. Reaction-of-progress integration and NH-source classification

The rate output records time, time increment, cycle number, reaction label, and reaction rate. For each accepted point, only records belonging to the final three contiguous P0-steady cycles are retained. The right-rule integrated source flux of pathway *j* is

`Phi_j = sum_k r_j(t_k) Delta t_k`,

where *k* indexes the accepted records for reaction family *j*. The pathway share is `100 Phi_j / sum_j Phi_j`. The denominator is restricted to the primary NH-forming families listed below so that the reported values answer a specific question: how is **positive NH formation** distributed among the selected direct sources?

| Label | Reaction family represented in the ROP integration |
|---|---|
| Vib-H2 | `N + H2(V1/V2/V3) -> H + NH` |
| Named-H2* | `N + H2(B3SIG/B1SIG/C3PI/A3SIG) -> H + NH` |
| Rydberg-H2* | `N + H2(RYDBERG_SUM) -> H + NH` |
| N(2D)+H2 | `N(2D) + H2 -> H + NH` |
| N(2P)+H2 | `N(2P) + H2 -> H + NH` |
| Association | `H + N + (N2 or H2) -> NH + (N2 or H2)` |

The ROP values are pathway-resolved model outputs, not directly measured reaction rates. They are reported together with the last-cycle window in the source data so that the exclusion of startup transients can be checked independently.

### 2.6. Diagnostic perturbations and endpoint checks

The principal 120 Td, N2=0.3333 CSTR reference point was evaluated in four isolated mechanism branches: an unmodified baseline, a branch with Rydberg-H2* to NH disabled, a branch with all four named electronic-H2* to NH reactions disabled, and a branch with both excited-N pathways disabled. Each branch received a copied kinetic input, was preprocessed and compiled independently, and was accepted only after the same P0 test.

A separate diagnostic set adds a direct Rydberg wall-relaxation loss at 1, 10, 100, or 1000 times a named-H2* wall-relaxation reference. This is a parameterized lifetime-sensitivity test, not a measurement-based Rydberg quenching model. It tests whether the conclusion about named electronic H2* depends on a single assumed lifetime for the Rydberg aggregate.

Two endpoint checks were also retained. The first evaluates 800 K at 60, 120, and 240 Td for N2=0.3333. The second evaluates fixed electron-density factors of 0.1, 1, and 10 at 120 Td, 300 K, and N2=0.3333. These are explicitly described as endpoint robustness checks, not as a dense temperature or electron-density parameter study.

### 2.7. Local reaction-control panel

To distinguish direct-NH-source dominance from product-level control, a targeted local kinetic-control panel was evaluated at the global NH3-productivity maximum (140 Td, N2=0.1) and at the high-field, N2-rich loss-dominated state (240 Td, N2=0.9). The panel contains eight rate groups: named electronic-H2* to NH, Rydberg-H2* to NH, N(2D,2P)+H2 to NH, H+N+M to NH+M, H+NH2+M to NH3+M, NH+H2+M to NH3+M, H+ + NH3 to NH3+ + H, and NH3+ + NH3 to NH4+ + NH2. Each group was multiplied in a copied kinetic input by 1.10 and by 1/1.10; all other CSTR/P0 settings were unchanged. The central logarithmic sensitivity of a response *y* to the group multiplier *k* is `S = [ln(y_1.10) - ln(y_1/1.10)] / [ln(1.10) - ln(1/1.10)]`.

The responses are CSTR outlet productivity, the ROP-integrated NH3 destruction-to-formation ratio, and the named electronic-H2* fraction of direct NH formation. All 34 calculations (two baselines plus two perturbation directions for eight groups at two conditions) passed the same P0 gate. This is a numerical local-control test, not a measurement-based uncertainty quantification of the elementary rate coefficients.

### 2.8. Targeted kinetic-control atlas and turnover-regime definitions

The two-point panel was extended to a targeted atlas with N2 inlet fractions of 0.1, 0.3, 0.5, 0.7, and 0.9 and E/N values of 20, 60, 100, 140, 180, and 240 Td. The three pre-specified groups were NH+H2+M to NH3+M, H+NH2+M to NH3+M, and H+ + NH3 to NH3+ + H. Central differences at 1/1.10 and 1.10 therefore required 180 perturbed calculations. A copied metrics-only version of the pulse driver retained the same kinetic equations, DVODE integration, terminal species output, and P0 convergence output, but suppressed the unused time-resolved ROP rows. Its 140 Td, N2=0.1 sensitivity for NH+H2+M to NH3+M agreed with the full-ROP calculation to `4.44×10^-16` in absolute sensitivity.

For Figure 12, the full 108-point turnover data set is classified descriptively into formation-starved, productive-hydrogenation, transition-turnover, and ion-loss-cancellation regimes. Ion-loss cancellation is defined as `D_NH3/P_NH3 >= 0.9`; productive hydrogenation as productivity at least 10% of the map maximum with `D/P < 0.5`; and formation starvation as productivity below 10% of the maximum with `D/P < 0.1`. Remaining accepted points are transition turnover. These thresholds are transparent descriptive boundaries for visual synthesis, not universal physical phase boundaries.

## 3. Results

### 3.1. The P0 gate establishes a common steady-state basis for all pathway comparisons

Figure 1 makes the numerical scope of the study explicit. The CSTR/CW calculation reaches a periodic state at the reference condition, whereas the matched closed-system regression continues to accumulate NH3 and remains outside the P0 threshold after the same finite time. The CSTR state and NH3 residuals both cross the `10^-3` acceptance threshold at the terminal stable window. Across the dense scan, every one of the 108 accepted points lies below the threshold; the two registered retry cases remain visibly identified rather than silently replaced. Consequently, all ROP integrations used below refer to equal three-cycle, period-steady windows rather than unequal transient endpoints.

**Figure 1.** Numerical and reactor-boundary validation for the CSTR/CW study. (A) 0D plasma-CSTR formulation and terminal P0 window. (B) cycle-end NH3 inventory for the CSTR and closed-system regression. (C) period-to-period convergence at the reference condition (120 Td, N2=0.3333). (D) P0 residual audit for the 108-point scan; crosses identify the two accepted retry cases.

### 3.2. The dense map rejects a dominant-pathway crossover within the accepted CSTR space

Figure 3 presents the full N2–E/N map. The main result is unambiguous within the defined model boundary: the Named-H2* family is the largest integrated NH source in every accepted case. Its share ranges from 84.892% to 99.664%. The lowest value occurs at 240 Td and N2=0.9, where the named electronic-H2* contribution nevertheless remains larger than every competing direct NH family. Thus, the present scan does not support wording such as “the dominant NH pathway switches from electronic excitation to vibration” or “a high-field nitrogen-rich crossover occurs.”

The absence of a dominant crossover does not imply that the network is invariant. Secondary sources respond continuously to field and composition. Rydberg-H2* reaches a maximum contribution of 4.353% at 240 Td and N2=0.1. N(2P)+H2 reaches 8.922% at 20 Td and N2=0.9. The association family reaches 7.906% at 240 Td and N2=0.9. These increases produce a measurable redistribution of the residual NH source, especially at nitrogen-rich conditions, but none overtakes the Named-H2* contribution. The largest low-vibrational-H2 fraction is only `4.45×10^-14%`, which excludes the low H2(v1–v3) family as a competitive direct NH source in this specific representation and scan.

The conclusion is strengthened by the numerical audit. All 108 points pass the P0 acceptance rule. The maximum `rel_state_max` in the full data set is `8.40×10^-4`, the maximum NH3 residual is `8.33×10^-4`, the cycle-energy residual is zero within the output precision, and every accepted case ends with a stable streak of three cycles. Consequently, Figure 3 is a map of accepted steady windows rather than a map of unequal transient durations.

**Figure 3.** Dense pathway map for the 108 accepted CSTR/CW cases: (A) steady NH3 density, (B) Named-H2* contribution to direct NH production, (C) Rydberg-H2* contribution, and (D) N(2P)+H2 contribution. Crosses indicate the two accepted retry cases. Source fractions are ROP integrals over the final three P0-steady cycles.

### 3.3. Reaction disabling establishes a strong causal role for named electronic H2*

High ROP share alone does not prove that a pathway materially controls the product inventory, because a reaction can be part of a highly recycled network. The isolated branch calculations address this issue at the 120 Td, N2=0.3333 reference point. In the unmodified branch, the terminal NH source is 99.020% Named-H2*, 0.641% Rydberg-H2*, 0.292% N(2P)+H2, and only a negligible fraction from low-vibrational H2.

Removing the four named electronic-H2* to NH reactions lowers the final NH3 density from `1.534×10^14` to `2.559×10^13 cm^-3`, an 83.31% reduction. Disabling the Rydberg-H2* source lowers NH3 by 0.60%, while disabling the N(2D) and N(2P) sources together lowers NH3 by 0.32%. The magnitude separation between these perturbations supports the interpretation that the Named-H2* pathways are not merely correlated with NH3 at the reference condition; they are the dominant causal inlet to the NH-containing portion of this model network.

| Diagnostic branch | Final NH3 (cm^-3) | Change relative to baseline | Interpretation |
|---|---:|---:|---|
| Baseline | 1.534×10^14 | — | Unmodified reference mechanism |
| Disable Rydberg-H2* to NH | 1.524×10^14 | -0.60% | Secondary pathway at the reference point |
| Disable Named-H2* to NH | 2.559×10^13 | -83.31% | Strong causal contribution to the NH3 network |
| Disable N(2D)/N(2P) to NH | 1.529×10^14 | -0.32% | Minor direct contribution at the reference point |

**Table 2.** Isolated diagnostic perturbations at 120 Td, N2=0.3333, 300 K, fixed electron density, and `tau=10 ms`. A disabled branch is a diagnostic counterfactual and not a replacement reaction mechanism.

The Rydberg-loss sensitivity further distinguishes the central result from an uncertainty in the aggregate Rydberg state. Increasing the added Rydberg direct-loss multiplier from 1 to 1000 decreases the Rydberg-H2* NH fraction from 0.628% to 0.029%. Over the same set, the Named-H2* fraction remains between 99.03% and 99.63%, and final NH3 changes only slightly from `1.532×10^14` to `1.525×10^14 cm^-3`. The test does not identify the physical Rydberg lifetime. Instead, it shows that a broad parameterized change to the direct Rydberg loss does not overturn the Named-H2* conclusion at this reference point.

**Figure 2.** Mechanism evidence at selected CSTR points. Panels A and B summarize steady NH3 and Named-H2* pathway fractions in the initial 3×3 subset; panel C shows the causal reaction-disabling calculation; panel D shows the Rydberg-loss sensitivity. The dense scan supersedes the sparse subset for the global N2–E/N conclusion.

### 3.4. The NH3 response has an intermediate-field window

The pathway map and product response must be interpreted separately. The final NH3 density across the dense grid ranges from `5.597×10^-4` to `2.557×10^14 cm^-3`. Figure 4A shows that NH3 is very low at 20–40 Td, rises sharply at moderate field, remains high over an intermediate-field band, and falls strongly at the highest fields. The global maximum, `2.557×10^14 cm^-3`, occurs at N2=0.1 and 140 Td.

The location of the maximum is systematic. For N2 fractions from 0.1 through 0.8, the largest accepted NH3 density occurs at 140 Td. At N2=0.9, the maximum moves to 120 Td and falls to `4.805×10^12 cm^-3`. The monotonic decrease in the peak NH3 density with increasing N2 fraction is shown in Figure 4B. At 240 Td, the final NH3 inventory is only `7.1×10^-5` to `1.7×10^-3` of the peak value at the same feed composition. This result establishes an intermediate-field product window under the imposed boundary; it does not constitute a calibrated energy-efficiency prediction for a physical reactor.

Figure 4C makes the topology–product distinction visible. Many high-NH3 and low-NH3 points have a Named-H2* NH fraction near 95–100%. Conversely, the low-field, nitrogen-rich points combine high or moderately high Named-H2* fractions with extremely low final NH3. Therefore, even a persistent dominant NH source does not determine the final product inventory. Subsequent NHx reactions, loss processes, and the CSTR residence-time balance remain essential to the product-level response.

**Figure 4.** NH3 field window and topology–product separation. (A) NH3 response to E/N at each inlet N2 fraction. (B) peak NH3 concentration and the field at which it occurs for each N2 fraction. (C) final NH3 versus the Named-H2* direct NH-source share; colour represents E/N and marker area represents N2 fraction. These are deterministic model outputs, so replicate error bars are not applicable.

### 3.5. Endpoint checks do not reveal a low-vibrational-H2 crossover

The dense grid is performed at 300 K and a fixed electron-density boundary. Two endpoint checks test whether the central conclusion is immediately reversed by moderate changes outside that grid. At 800 K and N2=0.3333, the Named-H2* share is 99.876%, 99.292%, and 96.442% at 60, 120, and 240 Td, respectively. The largest low-vibrational-H2 share at these endpoints is 0.00373%, which is still far below the dominant family.

At 300 K, 120 Td, and N2=0.3333, increasing the fixed electron density from 0.1 to 10 times its reference value changes the final NH3 density but retains the Named-H2* share between 98.756% and 99.020%. Figure 5 reports both endpoint tests. Their purpose is deliberately limited: they establish that the dense-grid observation is not trivially reversed at the tested temperature and electron-supply endpoints. They do not support a claim of global robustness across temperature, pulse frequency, duty cycle, or charge-self-consistent electron density.

**Figure 5.** Endpoint robustness checks. (A) NH-source fractions at 800 K as a function of E/N. (B) fixed-electron-density sensitivity at 120 Td and N2=0.3333. The figure is an endpoint check, not a full multidimensional scan.

### 3.6. Flux-weighted hypergraphs retain reaction arity across mechanistic regimes

Figure 6 resolves three representative operating states as species--reaction bipartite hypergraphs: a low-field state (20 Td, N2=0.1), the inventory/productivity maximum (140 Td, N2=0.1), and a high-field N2-rich state (240 Td, N2=0.9). Each coloured square is an explicit reaction hypernode, and each edge preserves a required co-reactant; its width is proportional to the terminal-three-cycle ROP. The drawings therefore avoid the common but misleading reduction of a multireactant elementary step to a simple pairwise arrow. The named electronic-H2* entry branch is prominent in every state, while the high-field graph visibly develops larger association and secondary-excitation edges before the downstream NHx network.

**Figure 6.** Flux-weighted reaction hypergraphs at low field, the NH3 maximum, and high field/N2-rich conditions. Species nodes, reaction hypernodes, co-reactants, and downstream NHx formation are explicitly retained. Edge widths are terminal-three-cycle ROP fluxes, not fitted pathway weights.

### 3.7. Graph topology supports, but does not substitute for, the ROP pathway ranking

Figure 7 separates three distinct forms of evidence. The dense field of the ratio of named-H2* flux to the largest secondary source makes the absence of a pathway crossover quantitative. Reaction-resolved and family-resolved flux panels show the hierarchy at the three representative states, while the N-lineage projection provides a topology-only weighted-betweenness diagnostic. The full reaction hypergraphs contain 434, 468, and 468 records at the low-field, optimum, and high-field representative states, respectively; all have zero N- and H-atom imbalance in the parsed records.

The symmetry ledger is intentionally conservative. Ground-state H2 is assigned its complete `X 1Sigma_g+` term, whereas B3SIG/A3SIG, B1SIG, and C3PI are classified only by the spin and orbital labels recorded by the Hong mechanism because g/u and +/- information is not carried by those abbreviations. The Rydberg aggregate is term-unresolved. Accordingly, group theory is used here to state the evidence boundary of the labels, not to impose a collision selection rule or infer an elementary-rate ordering. Terminal ROP integration supplies the pathway ranking.

**Figure 7.** Graph-theoretical and symmetry-aware interpretation of the reaction network. (A) dense direct-NH dominance ratio; (B) source-family decomposition; (C) N-lineage weighted betweenness; (D) reaction-level hierarchy; (E) state-symmetry evidence ledger; and (F) interpretation and atom-balance audit. Topological metrics do not replace flux-based pathway ranking.

### 3.8. NH3 accumulation follows a formation--destruction competition

Figure 8 moves from direct NH entry to product turnover. At the productivity maximum, the terminal NH3 formation and destruction fluxes are `3.562×10^16` and `1.024×10^16 cm^-3 s^-1`, respectively (`D/P=0.287`). At 240 Td and N2=0.9, those large opposing terms become `1.459×10^15` and `1.458×10^15 cm^-3 s^-1` (`D/P=0.999436`), leaving an outlet productivity of only `8.308×10^11 cm^-3 s^-1`. The leading high-field loss is H+ + NH3 -> NH3+ + H, followed by NH3+ + NH3 -> NH4+ + NH2. Thus, strong direct NH entry can coexist with almost complete product-turnover cancellation.

**Figure 8.** Terminal NH3 source--sink turnover. (A--C) largest elementary formation and destruction terms at three representative states, annotated with total formation `P`, destruction `D`, and `D/P`. (D--F) accepted 108-point fields of NH3 formation, destruction, and their ratio. All fluxes are ROP integrals over the terminal three P0-steady cycles.

### 3.9. Material-balance performance descriptors locate the productive window without claiming energy efficiency

Figure 9 maps CSTR NH3 outlet productivity, N2 conversion, N-atom utilization, and the productivity--conversion frontier. The global productivity maximum coincides with the intermediate-field inventory maximum (140 Td, N2=0.1), while the maximum N2 conversion is `0.2527%` at 240 Td and N2=0.1. N-atom utilization to NH3 peaks at `7.567%` at 120 Td and N2=0.1. These distinct extrema show why conversion, product outlet, and utilization should not be treated as interchangeable objectives. The descriptors are based on the imposed CSTR material balance; the fixed-electron-density power record is not used to report a physical energy efficiency.

**Figure 9.** CSTR material-balance performance map. (A) NH3 outlet productivity, (B) N2 conversion, (C) N-atom utilization to NH3, and (D) productivity--conversion frontier. The star marks the global productivity maximum; marker area in panel D denotes N2 inlet fraction.

### 3.10. Product-level kinetic control switches downstream of the direct NH gate

The local-control panel tests a different question from the direct-NH-source map: which selected rate groups change the net CSTR NH3 outlet after the complete NHx network responds? At the global productivity maximum, the largest local productivity sensitivity is `S=1.035` for NH+H2+M to NH3+M, followed by `S=0.229` for H+NH2+M to NH3+M; the named electronic-H2* to NH group has `S=0.0313`. At 240 Td and N2=0.9, H+NH2+M to NH3+M (`S=0.974`) and H+ + NH3 to NH3+ + H (`S=-0.856`) are opposing formation and loss levers. The `D/P=0.999436` state has small ratio sensitivities because formation and destruction co-vary near cancellation, making productivity the more discriminating response.

**Figure 10.** Local reaction-control analysis from symmetric 1/1.10 and 1.10 rate multipliers. (A) CSTR NH3-productivity sensitivities at the global productivity maximum (140 Td, N2=0.1). (B) productivity sensitivities in the high-field, N2-rich loss-dominated state (240 Td, N2=0.9). (C) sensitivities of the NH3 destruction-to-formation ratio at the latter state. (D) sensitivities of the named electronic-H2* direct-NH-source fraction. Positive values increase the stated response; negative values suppress it. All 34 cases passed the terminal P0 acceptance gate.

### 3.11. A targeted atlas separates formation control from the high-field loss lever

All 180 atlas calculations pass the terminal P0 gate. Across the 30 sampled states, NH+H2+M to NH3+M has productivity sensitivities from `-0.0374` to `4.1385`, H+NH2+M to NH3+M spans `-0.0235` to `1.7704`, and proton-driven NH3 ionization spans `-0.8556` to `0.0275`. NH+H2+M to NH3+M is the largest selected absolute local control at 22 points, H+NH2+M to NH3+M at seven, and ionization at one weak-control point. This count does not negate high-field loss: at 240 Td the proton-driven ionization sensitivity is negative at every composition (`-0.713` to `-0.856`).

**Figure 11.** Targeted kinetic-control switching atlas for NH3 CSTR productivity. (A--C) central logarithmic sensitivities of the three pre-specified groups to productivity, evaluated at five N2 fractions and six E/N values. (D) group with the largest absolute local sensitivity; a minus sign denotes negative control. Every cell is an explicit 1/1.10 and 1.10 calculation, without interpolation; black frames mark independently calculated `D/P >= 0.9` states.

### 3.12. A two-layer regime map links direct NH entry to NH3 turnover

The 108 accepted points separate into 31 formation-starved, 36 productive-hydrogenation, 13 transition-turnover, and 28 ion-loss-cancellation cases under the stated descriptive thresholds. This product-level partition is distinct from direct NH entry: named electronic-H2* remains the largest direct NH source in the sampled atlas, spanning 84.9--99.6%, including ion-loss-cancellation states. Figure 12 therefore combines an upstream electronic-H2* NH-entry gate with a downstream NHx-hydrogenation/product-loss competition. The conditional design implication is to operate where NHx hydrogenation provides positive productivity control and to avoid the `D/P >= 0.9` cancellation region; it is not an energy-efficiency optimization.

**Figure 12.** Two-layer mechanistic regime map. (A) turnover regimes from all 108 accepted CSTR/CW points. (B) largest selected local productivity control in the targeted 30-point atlas. (C) causal synthesis connecting persistent named electronic-H2* direct NH entry with operating-condition-dependent NHx hydrogenation and ion-mediated NH3 loss. Regime thresholds are descriptive definitions stated in Methods.

## 4. Discussion

The most defensible interpretation of the dense map is a conditional no-crossover result. Within the state representation of the corrected Hong mechanism and the present CSTR/CW boundary, the four named electronic-H2* reactions form a robust NH-production gate. Changing N2 fraction and E/N changes the amount of flux entering secondary families, but the reallocation remains insufficient to displace the named electronic-H2* family. This statement is stronger than reporting that one pathway is large at one reference point because it is tested across 108 accepted steady states and is supported by a reaction-disabling calculation.

The finding should not be read as a universal mechanism for plasma-assisted ammonia synthesis. The literature distinguishes multiple plasma-driven and plasma-enhanced catalytic regimes, whose relative importance can change with surface chemistry, temperature, pressure, power deposition, and discharge mode [Rouwenhorst2020]. Models that represent microdischarge structure and vibrational kinetics can also produce sensitivities that are absent under a fixed CW boundary [VantVeer2020]. The present work instead provides a reproducible baseline: if later simulations with a validated time-dependent electron boundary exhibit a pathway crossover, the difference can be attributed to an explicit change in boundary physics rather than to an unrecognized lack of periodic convergence.

The intermediate-field NH3 maximum offers a second lesson. Electronically excited H2* can remain the major direct NH source even when final NH3 collapses at high field. A useful mechanistic analysis must therefore separate at least three questions: which reactions form NH, which reactions subsequently move NH through the NHx network, and which reactions determine net NH3 accumulation under flow. The local-control panel makes this distinction quantitative. At the product maximum, downstream three-body NH-to-NH3 conversion provides the leading infinitesimal productivity control. At high field, NH3 formation through H+NH2+M competes directly with proton-driven NH3 ionization. Thus the persistent direct-NH-source topology and the switch in downstream product control are compatible, rather than competing, mechanistic statements.

The targeted atlas adds an important qualification. A control category should not be inferred from one representative state alone. NH+H2+M to NH3+M is the largest selected local productivity control across most of the sampled space, while H+NH2+M to NH3+M increases in relative importance toward nitrogen-rich, higher-field conditions. Proton-driven NH3 ionization is not always the largest absolute selected sensitivity, but it changes sign and becomes strongly negative across the entire 240 Td row. This is precisely why the turnover ratio and the local productivity sensitivities are shown together: `D/P` identifies the cancellation regime, and the sensitivity identifies a kinetic lever capable of moving the system away from it.

Several limitations define the scope of the manuscript. First, fixed electron density is an imposed electron-supply condition, not a self-consistent solution of the discharge. Second, the 0D CSTR averages the reactor and does not resolve spatial nonuniformity, physical wall transport, or local plasma–surface coupling. Third, the Rydberg aggregate and its perturbation are kinetic-model constructs; the lifetime sensitivity is not an experimental determination of Rydberg quenching. Fourth, the study contains no independent experimental validation. Finally, low-frequency, long-afterglow pulse cases generated DVODE warnings under the present electron boundary. They fail the same numerical validity standard used here and are intentionally excluded. A pulse paper should begin only after a time-resolved, conservative electron afterglow boundary is implemented and checked against tolerance, step-size, and electron-decay-time sensitivities.

These restrictions also point to useful experiments. Optical or laser-based diagnostics capable of tracking NH/NHx and electronic-state proxies during a controlled CSTR-like residence time would test the predicted separation between NH-source topology and outlet NH3. In modelling, the next priority is not a larger unvalidated duty-cycle map. It is a verified transient electron model that can preserve numerical and chemical consistency throughout the afterglow. The P0/ROP workflow presented here can then be reused without changing its evidence standard.

## 5. Conclusions

A CSTR/P0 workflow was used to map direct NH production routes in a corrected Hong N2–H2 plasma-kinetics mechanism. The dense scan contains 108 accepted CW cases at 300 K, a 10 ms residence time, and fixed electron density, with N2 inlet fraction from 0.1 to 0.9 and E/N from 20 to 240 Td. Four named electronically excited-H2 reactions are the largest NH source in every case, contributing 84.892–99.664% of the integrated direct NH production. The model therefore does not predict a crossover of the dominant NH source in this parameter space. Instead, high-field and nitrogen-rich conditions increase secondary Rydberg-H2*, N(2P)+H2, and association contributions while leaving Named-H2* dominant.

The product response follows a different pattern. Final NH3 is maximized at an intermediate reduced field, at 140 Td for N2=0.1–0.8 and at 120 Td for N2=0.9. The contrast between persistent Named-H2* dominance and a strong high-field decrease in NH3 demonstrates that direct NH-source topology alone cannot predict final product accumulation. The 30-point local-control atlas shows that NH+H2+M to NH3+M is the leading selected production lever at 22 sampled states, H+NH2+M to NH3+M at seven, while proton-driven NH3 ionization provides a strong negative high-field lever (`S=-0.713` to `-0.856` at 240 Td). Reaction disabling confirms the causal importance of the named electronic-H2* pathways at the reference point, and the Rydberg-loss diagnostic does not reverse the conclusion. The resulting data set is a numerically auditable steady-state baseline. It should not be extrapolated to pulse-afterglow switching until a validated time-dependent electron boundary is available.

## Data and code availability

The definitive Figure 1--12 vector and raster outputs are in `../analysis/figures_20260829/`; superseded duplicates are segregated in `legacy_pre_unified_redraw_20260830/` and are not manuscript figures. The 108-point terminal steady-state pathway table is available at `../analysis/figures_20260829/Figure_3_dense_EN_N2_source_data.csv`. The Figure 2 causal-perturbation source table is `../analysis/figures_20260829/Figure_2_CSTR_mechanism_map_source_data.csv`; the field-window table and peak conditions are `Figure_4_productivity_window_source_data.csv` and `Figure_4_peak_conditions.csv`. Figure 6 reaction fluxes and Figure 7 graph metrics are `Figure_6_reaction_flux_table.csv` and `Figure_7_hypergraph_metrics.csv`; Figure 8--9 turnover records are `Figure_8_9_NHx_CSTR_summary.csv` and `Figure_8_NHx_reaction_turnover.csv`. The 34-case local-control result and sensitivity tables are `../analysis/p6_local_reaction_control_20260830/P6_control_case_results.csv` and `../analysis/p6_local_reaction_control_20260830/P6_local_sensitivities.csv`. The 180-case control atlas and its 90-row sensitivity table are `../analysis/p7_control_switching_atlas_20260830/P7_control_cases.csv` and `../analysis/p7_control_switching_atlas_20260830/P7_control_sensitivities.csv`. Reproducible plotting scripts include `../工具脚本/plot_hong_dense_cstr_map.py`, `../工具脚本/plot_hong_graph_pathways.py`, `../工具脚本/plot_hong_nhx_turnover.py`, `../工具脚本/plot_hong_reaction_control.py`, and `../工具脚本/plot_hong_control_atlas_and_regimes.py`. The analysis excludes all pulse cases that failed the stated numerical-validity criterion.

## Declarations

### Ethics statement

Not applicable. This study contains no human participants, animals, clinical data, or personal data.

### Author contributions

To be completed by the authors using the CRediT taxonomy before submission. At minimum, the final record should identify contributors to conceptualization, methodology, software, validation, formal analysis, visualization, writing, supervision, and funding acquisition where applicable.

### Funding

To be completed by the authors. No funding source is inferred in this draft.

### Conflict of interest

To be declared by the authors before submission. No statement is inferred from the computational materials alone.

### AI-assisted writing disclosure

This internal draft was prepared with AI-assisted language and code-support tools. The authors must independently verify all scientific statements, numerical values, citations, authorship assignments, and venue-specific disclosure requirements before submission.

## References

Hong, J.; Pancheshnyi, S.; Tam, E.; Lowke, J. J.; Prawer, S.; Murphy, A. B. *Kinetic modelling of NH3 production in N2–H2 non-equilibrium atmospheric-pressure plasma catalysis.* **Journal of Physics D: Applied Physics** 2017, *50*, 154005. https://doi.org/10.1088/1361-6463/aa6229.

Hong, J.; Pancheshnyi, S.; Tam, E.; Lowke, J. J.; Prawer, S.; Murphy, A. B. *Corrigendum: Kinetic modelling of NH3 production in N2–H2 non-equilibrium atmospheric-pressure plasma catalysis.* **Journal of Physics D: Applied Physics** 2018, *51*, 109501. https://doi.org/10.1088/1361-6463/aaa988.

Rouwenhorst, K. H. R.; Engelmann, Y.; van ’t Veer, K.; Postma, R. S.; Bogaerts, A.; Lefferts, L. *Plasma-driven catalysis: green ammonia synthesis with intermittent electricity.* **Green Chemistry** 2020, *22*, 6258–6287. https://doi.org/10.1039/D0GC02058C.

Rouwenhorst, K. H. R.; Kim, H.-H.; Lefferts, L. *Vibrationally Excited Activation of N2 in Plasma-Enhanced Catalytic Ammonia Synthesis: A Kinetic Analysis.* **ACS Sustainable Chemistry & Engineering** 2019, *7*, 17515–17522. https://doi.org/10.1021/acssuschemeng.9b04997.

van ’t Veer, K.; Reniers, F.; Bogaerts, A. *Zero-dimensional modeling of unpacked and packed bed dielectric barrier discharges: the role of vibrational kinetics in ammonia synthesis.* **Plasma Sources Science and Technology** 2020, *29*, 045020. https://doi.org/10.1088/1361-6595/ab7a8a.

## Internal self-review checklist (remove before submission)

| Review dimension | Evidence or response in this draft | Status and required action |
|---|---|---|
| Contribution | The paper contributes a P0-gated CSTR workflow, a 108-point N2–E/N direct-NH-source map, and a falsifiable no-crossover result under explicit conditions. It does **not** claim that the result is universal or experimentally demonstrated. | **Needs literature strengthening.** Expand the related-work corpus before submission and compare the study's scope with the closest state-resolved kinetic analyses. |
| Writing clarity | The governing CSTR equation, excluded species, P0 criteria, scan values, ROP window, reaction families, retry treatment, and diagnostic branches are specified. Terminology distinguishes *direct NH-source share*, *final NH3 concentration*, and *fixed-electron-density boundary*. | **Pass for an internal draft.** A domain expert should still check mechanism nomenclature against the final kinetic input and journal notation. |
| Experimental strength | The manuscript reports deterministic model evidence over 108 accepted cases, causal reaction disabling, Rydberg-loss sensitivity, and two robustness endpoints. No experimental validation, spatially resolved calculation, or charge-self-consistent discharge calculation is presented. | **Needs new evidence for a high-impact claim.** Retain the paper as a conditional computational mechanism study unless independent diagnostics or experimentally anchored validation are added. |
| Evaluation completeness | The manuscript documents CSTR versus closed-batch behaviour, dense N2–E/N coverage, pathway disabling, a lifetime sensitivity, and temperature/electron-density endpoints. It does not contain dense temperature, electron-density, frequency, or duty-cycle scans; rejected pulse results are excluded. | **Needs new experiments for a pulse or broad-robustness paper.** Do not add a pulse-pathway conclusion without a validated afterglow electron boundary. |
| Method design soundness | The CSTR boundary removes closed-batch accumulation from steady-state pathway comparisons, but it remains a 0D, average-residence-time approximation with imposed electron density. These limits are stated in Methods and Discussion. | **Conditionally pass.** The title, abstract, conclusions, and figure captions must retain the CSTR/CW/fixed-ne qualification during all future revisions. |
