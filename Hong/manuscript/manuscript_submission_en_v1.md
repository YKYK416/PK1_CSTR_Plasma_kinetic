# Electronically Excited H2 Sustains Direct NH Entry while Downstream Turnover Controls NH3 in a Plasma CSTR

**Authors:** TBD  
**Affiliations:** TBD  
**Corresponding author:** TBD  
**Manuscript status:** English submission draft, version 1.7 (five-module Methods, literature-grounded introduction and publication-figure refresh), 6 September 2026

![Graphical abstract: bounded gas-phase CSTR result, Fe(110) structural alternatives, and the field-transfer comparison gate](../analysis/publication_figures_20260903_v2/Graphical_abstract_scope.png)

## Graphical abstract

The graphical abstract summarizes both the quantitative boundary and the forward validation logic. A continuous-wave N2/H2 plasma continuous stirred-tank reactor (CSTR) is scanned over feed composition and reduced electric field. Electronically excited H2* supplies the leading direct NH-entry family throughout the numerically accepted map in a literature-traceable reconstruction of the corrected Hong gas-phase mechanism, whereas the final NH3 response separates into a productive intermediate-field window and a high-field turnover-loss regime. A source-validated Fe(110) branch expands the structural surface alternatives, but a gas–surface flux comparison is intentionally withheld until local-field transfer and a shared reactor boundary are independently established.

## Abstract

Low-temperature plasma ammonia synthesis couples electron-impact, excited-state, radical, ionic, and transport chemistry, so the pathway that initiates N–H bond formation need not determine final NH3 production. Here, a literature-traceable reconstruction of the corrected Hong N2/H2 gas-phase mechanism was integrated in a continuous-wave, continuously fed and exhausted reactor over 108 numerically accepted states spanning N2 mole fractions of 0.1–0.9 and reduced electric fields of 20–240 Td. Terminal-window reaction-of-progress integration, direct-NH source-family classification, diagnostic family disabling, symmetric local kinetic perturbations, a targeted residence-time screen, a deterministic reaction-scenario envelope, and flux-weighted reaction-hypergraph analysis were used to separate NH entry from NH3 turnover. Within this explicitly bounded CSTR model, the electronically excited H2* family was the largest integrated direct source of NH at every accepted map point, contributing 84.892–99.664% of direct NH production. At 120 Td and x_N2 = 0.3333, disabling this family decreased NH3 by 83.31%, whereas the tested Rydberg-labelled and excited-N alternatives produced sub-percent changes. In 20 additional P0-accepted cases spanning 0.5–30 ms at four contrasting states, the named-H2* direct-NH share remained 76.832–99.596% and its dominance margin remained above one. In 14 additional P0-accepted 1/3–3× rate-group scenarios at two reference states, the minimum dominance margin was 10.724. No direct-NH crossover was therefore observed in the composition–field map or in these targeted transport and reaction-scenario extensions. NH3 performance nevertheless followed a distinct downstream control layer: inventory and outlet productivity maximized at 140 Td and x_N2 = 0.1, whereas at 240 Td and x_N2 = 0.9 gross formation and destruction nearly cancelled (D/P = 0.999436) through proton-driven ammonia loss. A source-validated Fe(110) DFT–microkinetic branch is then used only to define explicit Eley–Rideal and Langmuir–Hinshelwood surface alternatives and a falsifiable boundary-harmonization protocol; it does not supply a cross-boundary flux ranking. The resulting 26 figures support a model-specific two-layer mechanism and a disciplined route toward gas–surface discrimination. The result is a reproducible baseline for falsification by self-consistent electron kinetics, pulsed operation, surfaces, and experiment; it is not a universal pathway claim for all plasma-ammonia reactors.

**Keywords:** plasma ammonia synthesis; N2/H2 plasma; chemical kinetics; reaction-of-progress analysis; excited hydrogen; reaction networks; CSTR; kinetic control

## 1. Introduction

Ammonia is one of the chemical foundations of modern society. Its primary use remains the manufacture of nitrogen fertilisers, yet its high hydrogen content, carbon-free molecular composition, established storage infrastructure, and ease of liquefaction also make it a prospective energy carrier and a chemical store for renewable electricity. These latter roles have made the route by which nitrogen is fixed increasingly important. The conventional Haber–Bosch process is exceptionally mature, productive, and difficult to displace at large scale. Nevertheless, its customary coupling to centralised hydrogen production, high-pressure synthesis loops, and continuously operated thermal equipment poses a mismatch with geographically distributed renewable electricity whose availability changes on much shorter time scales. This mismatch does not make Haber–Bosch obsolete; it motivates complementary synthesis routes that may operate at a smaller scale, at milder bulk conditions, and with rapid electrical controllability [1,2]. Plasma-assisted nitrogen fixation is therefore best viewed not as a simple replacement claim, but as an effort to identify circumstances in which electrical activation can create a useful chemical degree of freedom that conventional thermal operation does not possess.

The central chemical obstacle is the extraordinary stability of molecular nitrogen. In a thermal catalytic synthesis loop, the gas temperature and catalyst surface jointly determine the population of molecules able to surmount the N2 activation barrier, the coverage of nitrogen-containing intermediates, the rates of successive hydrogenation, and the degree to which the reverse reaction erodes the product yield. The resulting requirement for elevated temperature and pressure is a system-level compromise among kinetics, equilibrium, separation, heat integration, and catalyst durability. Non-thermal plasma changes the nature of that compromise. A relatively small electron population can acquire an energy distribution that is far from the translational temperature of the neutral gas. Electron-impact collisions can then populate electronic, vibrational, and dissociative channels of N2 and H2 while most heavy species remain comparatively cool. In principle, this permits a pathway to begin with activated reactants rather than with a uniformly heated gas [1,3,17]. The appeal is real, but it should not be confused with a guarantee of high energy efficiency, high single-pass conversion, or a universal mechanistic advantage.

The phrase “plasma activation” conceals a hierarchy of coupled processes. The electron energy distribution function (EEDF) determines the relative probabilities of excitation, ionisation, dissociation, attachment, and superelastic collisions. The products of those processes include vibrationally excited molecules, electronically excited molecules, radicals, ions, electrons, and photons. These species react with one another during a discharge, persist or decay in the afterglow, diffuse to walls, and, in packed reactors, encounter surfaces whose temperature, charging, adsorption state, roughness, pore geometry, dielectric character, and local electric field may all differ from those of the bulk gas. The same nominal input power can consequently represent very different microscopic states in a filamentary dielectric-barrier discharge (DBD), a diffuse discharge, a radio-frequency plasma, a glow discharge, or a pulsed plasma. A reactor-level NH3 concentration is an indispensable performance measure, but it is the integrated outcome of this entire hierarchy rather than a direct measurement of the first elementary N–H bond-forming event [3,17–19].

This distinction is particularly consequential for ammonia. A molecule detected at the outlet may have been formed through gas-phase radical chemistry, through a plasma-modified surface cycle, through both in sequence, or through a route whose decisive intermediate was produced in a short microdischarge and converted only in the afterglow. Conversely, a large instantaneous rate of an NH-forming elementary reaction does not necessarily imply a large contribution to outlet NH3. The NH intermediate can be hydrogenated, dissociated, ionised, recycled to molecular nitrogen, converted to other nitrogen–hydrogen species, or transported out of the active volume before it reaches ammonia. Product-level observations alone therefore cannot establish a dominant elementary pathway without a reaction network, a stated reactor boundary, and a time-resolved or flux-resolved interpretation. This is one reason why apparently different mechanistic accounts in the plasma-ammonia literature can each be plausible within their own physical and modelling boundary.

Early and continuing experimental work has established that ammonia can be produced in plasma-only reactors and in plasma–catalyst systems under bulk conditions far milder than those conventionally used in thermal synthesis. DBD reactors dominate much of this literature because they can be operated at atmospheric pressure, tolerate frequent power modulation, and are readily combined with packed beds. Reactor studies have examined the effects of feed composition, flow rate, pressure, applied voltage, power, residence time, packing material, and catalyst identity. They have also demonstrated that an added catalyst is not a passive replacement for an empty discharge volume: it can change the discharge mode, local electric field, filament propagation, gas heating, surface area, adsorption capacity, and the spectrum of reactive species available near the solid. The resulting comparison between “plasma only” and “plasma catalysis” is therefore mechanistically informative only when those concurrent changes are considered [1,11,17,18].

This experimental diversity has generated several recurring observations, but not a single universal route. In some DBD studies, ammonia productivity responds strongly to the catalyst and support, suggesting that the treatment of activated nitrogen and hydrogen at the surface is central. In others, appreciable ammonia is observed in the absence of a catalyst, establishing that gas-phase chemistry cannot be assumed negligible. Packed-bed systems often exhibit changes in selectivity or energy yield that cannot be explained by thermal catalysis alone, yet those changes may arise from combinations of adsorption, plasma-induced desorption, dielectric enhancement, microdischarge redistribution, and surface reactions. Radio-frequency and other plasma sources introduce another set of electron densities, gas temperatures, characteristic residence times, and power-coupling modes. Elevated-pressure experiment–model work has further indicated that nitrogen electronic excitation and vibrationally excited hydrogen can influence different portions of the network under the same broad N2/H2 feed chemistry [8,12,13]. The appropriate conclusion is not that one observation invalidates another; it is that pathway labels must always be accompanied by the discharge, surface, pressure, and state-space definitions from which they were inferred.

A useful way to organise the proposed mechanisms is by the physical location and electronic state of the first activation event. The most direct gas-phase family begins with electron-impact dissociation of N2 and H2 to N and H atoms, followed by nitrogen–hydrogen association reactions such as N + H2 -> NH + H or N + H + M -> NH + M, depending on the mechanism. Excited N2 can modify this family by altering the availability of N atoms or by participating in state-specific reactions. Vibrationally excited H2 can, in turn, affect the balance between H-atom formation, hydrogen abstraction, and reactions of NHx intermediates. Other models include electronically excited hydrogen states as explicit intermediates, allowing a reactive H2* pool to alter the direct formation of NH. In a mechanism that resolves these states, a label such as “gas-phase route” is still insufficient: the identity of the state, its source term, its quenching channels, and the cross section or rate coefficient used to represent it can determine whether that route is important.

The corresponding surface families are normally grouped as Langmuir–Hinshelwood (LH), Eley–Rideal (ER), and related radical-assisted sequences. In LH chemistry, reactants adsorb before reacting, and a surface nitrogen species may undergo successive hydrogenation to NH*, NH2*, and NH3*. In ER chemistry, a gas-phase radical or excited reactant collides directly with an adsorbed partner. Plasma-created N atoms, H atoms, vibrationally excited N2, and excited hydrogen species can therefore change the surface reaction landscape even when the mean gas temperature is low. The distinction matters because the same catalyst can be optimised differently for a thermal LH cycle, a vibrationally assisted N2 dissociation process, and a radical-rich ER process. Microkinetic and density-functional-theory studies have consequently predicted that catalyst ranking, the apparent relevance of N2 vibrational excitation, and the relative importance of ER and LH contributions can shift as the assumed plasma species densities and vibrational distributions are changed [4,5,10,16,20,21].

The mechanistic importance of plasma-generated radicals has also been framed in terms of how plasma activation modifies an equilibrium-limited catalytic reaction. A minimal plasma-enabled microkinetic model, supported by packed-bed observations, showed that plasma-activated nitrogen can enable ammonia yields above the corresponding bulk thermal-equilibrium limit when the thermal route is kinetically restricted by N2 dissociation [20]. This result is conceptually important: a non-equilibrium plasma can supply activated reactants to a catalyst rather than merely acting as a heat source. At the same time, the result does not identify the detailed mechanism of every DBD reactor. It rests on a specified simplified representation of plasma activation and a thermal microkinetic framework. It therefore motivates, rather than removes, the need to determine which activated species are actually present and which elementary steps receive their flux in a particular discharge.

The time structure of atmospheric-pressure plasma operation adds another layer of difficulty. In filamentary DBDs, power is not deposited uniformly in a steady homogeneous volume. Short-lived microdischarges can produce electron-impact excitation and dissociation in spatially concentrated channels, after which reactive neutrals and surface intermediates evolve in a longer afterglow. A detailed model of a plasma-catalytic DBD found that ammonia could be decomposed during microdischarges while net formation occurred in the afterglow, illustrating why a temporally averaged input power or a time-averaged rate alone can obscure the chemical sequence [22]. The implication extends to pulsed operation. Frequency, duty cycle, pulse width, rise time, repetition rate, and the inter-pulse relaxation period can each change the balance among EEDF-controlled production, radical recombination, vibrational relaxation, wall loss, adsorption, and desorption. A continuous-wave calculation can provide a useful limiting case, but it is not automatically a surrogate for a pulsed discharge without an explicit temporal model.

Accordingly, the reported dependence of ammonia synthesis on N2 fraction, reduced electric field E/N, power, electron density, temperature, and residence time should not be read as a collection of independently interchangeable knobs. The reduced field is closely tied to electron heating and hence to electron-impact rate coefficients, but its effect depends on the gas composition, pressure, and EEDF model. Electron density controls the magnitude of electron-driven source terms only together with the corresponding state-dependent rate coefficients. Gas temperature affects density, heavy-particle kinetics, equilibrium tendencies, quenching, diffusion, and surface processes. Feed composition changes dilution, electron-energy loss channels, vibrational populations, radical pools, and the availability of both N- and H-containing collision partners. Residence time determines whether fast plasma production, slower neutral chemistry, and product removal are observed within the same reactor window. In a packed bed, the surface temperature and local field can also be decoupled from their nominal bulk values. It follows that a parametric map is most informative when it varies selected variables under an explicit set of held-fixed assumptions, rather than presenting a response as though it were an intrinsic property of nitrogen alone.

These issues have encouraged modelling across several levels of description. Global or zero-dimensional plasma-chemical models are valuable because they can include many electron-impact and heavy-particle reactions, scan broad parameter spaces, and connect reaction-of-progress analysis to outlet species concentrations. They require, however, a closure for the EEDF, a representation of power deposition or electron density, a reactor model, and rate coefficients for all included state-specific processes. Spatially resolved fluid, kinetic, and particle models can represent microdischarge structure and local fields more directly, but are computationally more demanding and still require chemical data and boundary conditions. Surface microkinetic models add adsorption, coverage, and catalyst-specific elementary steps; coupled plasma–surface models then require an interface between gas-phase fluxes, surface kinetics, and the discharge physics. Recent reviews of plasma-ammonia modelling identify precisely this multiscale closure problem: the literature contains increasingly capable models, but their predictions remain sensitive to reactor geometry, plasma source, catalyst representation, and the chemical state space selected by the modeller [14,18,19].

The choice of state space is especially important when the scientific question concerns a pathway crossover. A conventional species-only network may represent N2, H2, N, H, NH, NH2, and NH3, whereas a state-resolved network distinguishes vibrational manifolds, metastables, electronically excited molecules, ions, and electronically excited hydrogen states. Adding such states may reveal physically plausible routes, but it also introduces uncertainties in excitation, quenching, pooling, wall-loss, and reaction data. Conversely, omitting a state may force its chemical influence into an effective rate coefficient or remove a route altogether. This is not a reason to avoid mechanistic modelling. It is a reason to interpret “dominant pathway” as a conditional statement: dominant within a declared reaction set, kinetic data set, transport model, and operating domain. A rigorous model study should expose those conditions rather than relying on the apparent precision of a colourful pathway diagram.

Several recurrent limitations follow from this observation. First, final NH3 concentration or energy yield cannot by itself distinguish an upstream NH-entry route from a downstream conversion or loss bottleneck. Two operating points may form NH through the same reaction family but produce different NH3 amounts because their NH hydrogenation, NH3 destruction, or residence-time terms differ. Second, ranking the largest instantaneous rate at one time and one state can be misleading in a discharge with transients or in a continuous-flow reactor that has not reached an acceptably stationary terminal window. Integrated positive fluxes, source–sink budgets, and numerical-convergence tests are needed before an apparent pathway ordering is treated as a physical result. Third, a pathway drawn from a gas-phase model and a pathway inferred from a surface-sensitive experiment are not independent confirmations unless the shared assumptions and missing processes are made explicit. Fourth, changes in E/N, electron density, power, catalyst, and temperature frequently alter more than one physical mechanism at once, limiting the interpretation of one-factor trends [14,17–19].

A further limitation is that mechanistic analysis has often been presented as a list of selected reaction rates rather than as a reproducible network problem. Large plasma mechanisms contain parallel reactions, state-specific branches, loops, and competing sinks. A reaction-of-progress table remains essential because it preserves the signed kinetic flux, but it is difficult to use such a table alone to distinguish a genuinely persistent route from a visually prominent yet low-flux branch. Reaction hypergraphs and flux-weighted graphs offer a complementary representation: reactions can be retained as first-class objects, species can be connected by the reactions that consume and create them, and metrics such as weighted degree, betweenness, path coverage, and community structure can be evaluated alongside integrated fluxes. Graph theory does not create causal evidence by itself; its value is to make the topology of a stated chemical model auditable and to identify where targeted kinetic-control tests should be placed [15].

The need for explicit mechanistic boundaries is sharpened by the availability of influential, literature-traceable reaction mechanisms. Hong et al. assembled a widely used atmospheric-pressure N2/H2 plasma-catalysis mechanism that couples electron-driven activation with a detailed heavy-particle chemistry [6]. The subsequent corrigendum made material changes to the implementation, including the rate expression for N + H2(v) -> NH + H, and excluded selected ionic reactions because of numerical instability and minor importance under the originally reported conditions [7]. The corrected mechanism is therefore a valuable and reproducible basis for a controlled gas-phase analysis, but it is not a universal description of every plasma-ammonia reactor. It does not, without a dedicated extension, establish the behaviour of a fully state-resolved EEDF, a catalyst-specific surface network, a spatially resolved microdischarge, or an arbitrary pulsed waveform. Treating it as a defined boundary is both more conservative and more scientifically useful than treating it as an all-purpose mechanistic truth.

Against this background, the present study asks a deliberately narrow but consequential question: within a literature-traceable reconstruction of the corrected Hong gas-phase N2/H2 mechanism and an explicitly defined plasma continuous stirred-tank reactor (CSTR), does the leading direct route to NH switch as nitrogen fraction and E/N are varied over a dense operating map? Direct NH formation is an appropriate focal point because it is the first nitrogen–hydrogen bond-forming entry into the NHx manifold. Its identity constrains the species and excitation channels that must be supplied upstream, while its subsequent fate determines whether that entry becomes outlet NH3. The question is thus framed as two coupled but non-identical layers: an NH-entry layer, defined by direct NH source families, and an NH3-turnover layer, defined by formation, recycling, loss, and transport competition downstream of NH.

The analysis is designed to reduce several common ambiguities in pathway attribution. First, every accepted continuous-wave operating point is subjected to the same terminal-window numerical-convergence gate before reaction pathways are compared. This prevents an apparent crossover from being generated by unequal transient relaxation across the map. Second, direct NH sources are classified into reaction families and integrated over the terminal window across a dense x_N2–E/N grid rather than being inferred from a single instantaneous rate. Third, the apparent contribution of electronically excited H2 is tested by disabling the relevant reaction family at a reference state, thereby separating a flux correlation from a controlled mechanistic intervention within the model. Fourth, NH-entry topology is kept distinct from NH3-level source–sink accounting, local kinetic controls, and product turnover. Finally, reaction-resolved hypergraphs and graph metrics are used as structured visual evidence alongside, not in place of, reaction-of-progress ranking and integrated flux analysis.

There is also an important evidentiary distinction between quantifying ammonia and assigning a route to ammonia. Accurate product analysis requires attention to calibration, adsorption and desorption, material blanks, nitrogen source controls, moisture, background contamination, and the time at which the product is sampled. These experimental safeguards establish whether a reported NH3 signal is genuine and quantitatively reliable; they do not on their own reveal which state-resolved elementary reaction first produced NH. Conversely, a model can resolve every reaction it contains without proving that its assumed excited-state density, surface coverage, wall-loss probability, or EEDF closure is realised experimentally. The most convincing mechanistic programme therefore combines disciplined product quantification, electrical and optical diagnostics, state-sensitive measurements where feasible, and models whose assumptions are accessible to falsification. Reviews centred on diagnostic tools and multiscale modelling have emphasised that correlating discharge properties, reactive species, and surface processes remains a central unresolved need in DBD plasma catalysis [17–19].

For this reason, the pathway question cannot be reduced to a binary choice between “gas phase” and “surface.” A radical generated in the plasma can adsorb and enter an ER reaction; an adsorbed intermediate can desorb and be processed in the gas; a packing can reshape the discharge while simultaneously providing a reaction surface. The experimentally relevant object is often a coupled gas–surface reaction network. Yet retaining that insight should not encourage over-interpretation of a model that contains only one side of the interface. A gas-phase CSTR model can identify which gas reactions are internally responsible for its NH and NH3 fluxes. A surface microkinetic model can compare elementary surface cycles under its supplied gas-phase boundary. Neither calculation, without a shared flux boundary and a consistent representation of local fields and transport, supplies a quantitative ranking between gas and surface branches. The present work therefore uses the gas-phase result to specify what a future cross-boundary comparison must reproduce, rather than using a qualitative surface alternative as retrospective proof of a gas-phase pathway.

Dense parametric mapping is valuable in this setting for a different reason than simple process optimisation. A small number of chosen operating points may identify a promising productivity window while missing a pathway boundary, an extinction region, a numerical transient, or a change in the dominant sink. Sampling N2 fraction coarsely can obscure the trade-off between the supply of nitrogen-containing precursors and hydrogen-rich collision partners. Sampling E/N coarsely can hide shifts in the EEDF-controlled excitation hierarchy. The same caution applies to frequency, duty cycle, temperature, electron density, and residence time: a response surface is meaningful only after the physical meaning of each scan variable and the quantities held fixed are stated. The present composition–field map is thus not presented as an optimised operating envelope for a particular device. It is a deliberately controlled test of whether the leading direct NH source changes within the corrected reaction mechanism when its two most immediate gas-phase control variables are scanned systematically.

Finally, the use of a two-layer mechanism has implications for how negative results are reported. A finding that no direct-NH crossover occurs over an accepted grid is not an absence of chemistry and is not a statement that all operating variables are unimportant. It specifies that, for the tested reaction system, the observed changes in NH3 cannot be attributed to a switch among the classified direct NH-entry families. The explanatory effort must then move downstream, toward the balance of NH consumption, NHx conversion, ammonia formation, ammonia destruction, and flow removal. This separation guards against a common circular interpretation in which the reaction feeding a named intermediate is automatically assigned responsibility for the final product trend. It also provides a compact falsification strategy: if experiments, self-consistent EEDF calculations, pulsed simulations, or surface-coupled models show a different NH-entry ordering, the discrepancy can be traced to a defined change in the model boundary rather than to an undefined disagreement about “plasma effects.”

The scope is intentionally limited. The present calculations are not an experimental validation of a universal plasma-ammonia pathway, do not prescribe a catalyst design, and do not claim that electronically excited H2 is the dominant NH precursor in every plasma source or material environment. They also do not equate a continuous-wave CSTR with a filamentary or pulsed DBD. Instead, they provide a reproducible mechanistic baseline: under a corrected and explicitly bounded gas-phase mechanism, which NH-entry family persists across a dense composition–field map, which pathways do not exhibit a resolved crossover, and which downstream processes control whether upstream NH formation becomes NH3. Such a baseline is valuable because it defines the specific chemical, kinetic, temporal, or surface assumptions that must be changed before an alternative regime can be credibly claimed.

## 2. Methods

### 2.1. Kinetic model and scope

The calculations used a literature-traceable reconstruction of the corrected Hong \(\mathrm{N_2/H_2}\) non-equilibrium atmospheric-pressure plasma mechanism and its subsequent corrigendum [6,7]. The mechanism couples electron-impact excitation and dissociation to a state-resolved heavy-particle network. It contains ground-state \(\mathrm{N_2}\) and \(\mathrm{H_2}\), atoms, \(\mathrm{NH_x}\) intermediates, ions, electronically excited species, selected vibrationally excited states, electron-impact reactions, binary reactions, three-body reactions, and relaxation or loss processes. Electron-impact rate coefficients were generated using the ZDPlasKin–BOLSIG+ workflow for each prescribed reduced electric field \(E/N\) and inlet composition.

This framework resolves the internal kinetic consequences of prescribed \(E/N\) and electron-density boundaries. It does not solve a self-consistent electron energy distribution function, electrical circuit, power deposition, sheath structure, or spatial discharge dynamics. The electron density was imposed as an external parameter rather than predicted from charge balance or transport.

The present implementation is deliberately restricted to the gas-phase kinetic network. It contains no catalyst-specific adsorption, desorption, surface diffusion, surface-site balance, Langmuir–Hinshelwood pathway, Eley–Rideal pathway, wall-material-specific reaction probability, axial gradient, microdischarge channel, or local surface-field solution. All pathway rankings therefore refer only to the reconstructed gas-phase mechanism under the CSTR boundary defined below; they cannot establish the relative contribution of gas and surface pathways in a packed-bed reactor.

The term “named \(\mathrm{H_2^*}\)” is an input-mechanism label for the four electronically labelled precursors \(\mathrm{H_2(B3SIG)}\), \(\mathrm{H_2(B1SIG)}\), \(\mathrm{H_2(C3PI)}\), and \(\mathrm{H_2(A3SIG)}\). These species participate in the direct \(\mathrm{NH}\)-forming reactions

\[
\mathrm{N + H_2(state) \rightarrow NH + H}.
\]

The label identifies a reproducible set of elementary reaction records; it is not a new spectroscopic assignment, a measured excited-state density, or a universal claim regarding plasma ammonia reactors.

### 2.2. CSTR model, imposed plasma parameters, and operating map

The gas-phase mechanism was integrated in a zero-dimensional, perfectly mixed continuous stirred-tank reactor (CSTR) with continuous feed and exhaust. For each heavy-particle species \(i\), the number-density balance was

\[
\frac{dn_i}{dt}
=
\sum_r\nu_{ir}R_r
+
\frac{n_{i,\mathrm{in}}-n_i}{\tau},
\]

where \(n_i\) is number density, \(\nu_{ir}\) is the signed stoichiometric coefficient of reaction \(r\), \(R_r\) is its directional reaction-of-progress rate, \(n_{i,\mathrm{in}}\) is the inlet number density, and \(\tau\) is residence time. The first term represents chemical production and consumption, while the second represents inlet replenishment and outlet removal. Feed species relax toward their imposed inlet densities, whereas products, radicals, excited species, and ions absent from the feed are removed through the outlet after formation. Electron density was prescribed independently and was not solved as a feed or outlet species balance.

The primary boundary was

\[
T_{\mathrm{g}}=300~\mathrm{K}, \qquad
\tau=10~\mathrm{ms}, \qquad
n_{\mathrm{e}}=1.17\times10^8~\mathrm{cm^{-3}}.
\]

Calculations used continuous-wave forcing: within each case, the electron-density boundary and electron-impact rate-coefficient set were held fixed. Thus, continuous wave does not represent a resolved voltage waveform, pulse frequency, duty cycle, afterglow period, or microdischarge sequence. The imposed \(E/N\) changes electron-impact rate coefficients, but it is not an applied voltage, absorbed power, or local surface field. Likewise, the imposed \(n_{\mathrm{e}}\) scales electron-driven source terms but is not a predicted discharge property. The model therefore supports controlled pathway comparisons but does not close a power balance; reported \(\mathrm{NH_3}\) concentration, outlet productivity, conversion, and nitrogen utilisation are material-balance descriptors rather than energy-performance metrics.

**Table 1. Primary computational boundary and parameter scan.**

| Quantity | Value or treatment |
|---|---|
| Reactor model | Zero-dimensional, perfectly mixed CSTR with continuous feed and exhaust |
| Gas temperature | \(T_{\mathrm{g}}=300~\mathrm{K}\) |
| Residence time | \(\tau=10~\mathrm{ms}\) |
| Electron density | \(n_{\mathrm{e}}=1.17\times10^8~\mathrm{cm^{-3}}\), imposed |
| Forcing | Continuous wave; no physical pulse waveform in the primary map |
| \(\mathrm{N_2}\) inlet mole fraction | \(x_{\mathrm{N_2}}=0.1\)–\(0.9\) in increments of \(0.1\) |
| Reduced electric field | \(E/N=20\)–\(240~\mathrm{Td}\) in increments of \(20~\mathrm{Td}\) |
| Primary grid | \(108\) CSTR states |
| Primary exclusions | Self-consistent power/EEDF evolution, pulse frequency or duty-cycle effects, and catalyst surface chemistry |

The map combines nine nitrogen fractions with twelve \(E/N\) values, producing

\[
9\times12=108
\]

primary CSTR states. Its purpose is to test whether a direct-\(\mathrm{NH}\) source crossover occurs within the fixed gas-phase model, not to fit an experimental response surface or optimise a practical reactor. Temperature and electron-density endpoints were evaluated only as targeted robustness checks; they do not constitute a complete uncertainty quantification or replace self-consistent electron kinetics.

### 2.3. Numerical integration and terminal-window acceptance

The stiff kinetic system was integrated with the generated ZDPlasKin module and its associated implicit solver. Each calculation retained terminal species densities, reaction-resolved ROP records, time increments, numerical-window diagnostics, input parameters, and run-status metadata. The pulse-capable driver was operated in continuous-wave mode. Thus, the term “window” denotes a repeated numerical analysis interval, not a physical voltage pulse, microdischarge, or afterglow.

All pathway comparisons were gated by the terminal-window acceptance criterion \(P0\). For consecutive windows \(k-1\) and \(k\), the maximum relative state residual was calculated as

\[
\varepsilon_{\mathrm{state},k}=\max_i\left[\frac{|n_{i,k}-n_{i,k-1}|}{\tfrac12|n_{i,k}+n_{i,k-1}|+n_{\mathrm{floor}}}\right],
\]

where \(n_{\mathrm{floor}}\) is a small numerical floor that prevents near-zero species from producing undefined ratios. An analogous residual was evaluated for \(\mathrm{NH_3}\),

\[
\varepsilon_{\mathrm{NH_3},k}
=
\frac{|n_{\mathrm{NH_3},k}-n_{\mathrm{NH_3},k-1}|}
{\tfrac12|n_{\mathrm{NH_3},k}+n_{\mathrm{NH_3},k-1}|+n_{\mathrm{floor}}}.
\]

The \(P0\) gate required the state residual, \(\mathrm{NH_3}\) residual, and stored window-energy diagnostic to satisfy

\[
\varepsilon\leq10^{-3}
\]

for three consecutive terminal windows. Only the final three accepted windows were used for ROP integration and source-sink analysis.

The stored window-energy quantity is retained as an execution diagnostic. Under fixed continuous-wave forcing it is identically zero within output precision and is not an energy-conservation residual, an energy balance, or an efficiency calculation. \(P0\) establishes numerical comparability of terminal windows; it does not validate the chemical database, the prescribed electron-density boundary, or the physical completeness of the zero-dimensional model. All \(108\) primary states passed this gate, including two documented numerical retries. Cases that did not pass \(P0\) were excluded rather than replaced by neighbouring conditions or transient ROP values.

### 2.4. ROP integration, direct-NH classification, NH3 turnover, and robustness tests

Each ROP record contains a reaction identifier, rate \(R_r\), time increment \(\Delta t\), and numerical-window index. Integrated flux for family \(f\) was evaluated over the final three \(P0\)-accepted windows \(\mathcal W\) by

\[
\Phi_f=\sum_{j\in\mathcal W}\sum_{r\in f}R_r(t_j)\Delta t_j.
\]

Records with missing or non-finite time, rate, or time increment were rejected, as were records with non-positive \(\Delta t\). Family assignment was based on pre-specified, normalised reaction selectors rather than visual inspection of pathway diagrams.

Six direct-\(\mathrm{NH}\) families were tracked: vibrational-\(\mathrm{H_2}\) reactions \(\mathrm{N+H_2(V1,V2,V3)\rightarrow H+NH}\); the four named-\(\mathrm{H_2^*}\) reactions; \(\mathrm{N+H_2(RYDBERG\_SUM)\rightarrow H+NH}\); \(\mathrm{N(^2D)+H_2\rightarrow H+NH}\); \(\mathrm{N(^2P)+H_2\rightarrow H+NH}\); and termolecular association \(\mathrm{H+N+M\rightarrow NH+M}\), with \(M=\mathrm{N_2}\) or \(\mathrm{H_2}\). The direct-\(\mathrm{NH}\) share was

\[
S_f=\frac{\Phi_f}{\sum_g\Phi_g}.
\]

The named-\(\mathrm{H_2^*}\) dominance margin was defined as

\[
M_{\mathrm{H2^*}}=\frac{\Phi_{\mathrm{named\ H2^*}}}{\max_{g\ne\mathrm{named\ H2^*}}\Phi_g}.
\]

A direct-\(\mathrm{NH}\) crossover required a secondary family to exceed the named-\(\mathrm{H_2^*}\) family in the same \(P0\)-accepted terminal interval, i.e.

\[
M_{\mathrm{H_2^*}}<1.
\]

Shannon entropy and the corresponding effective number of families were computed only as descriptive summaries of the same six integrated source shares; they were not treated as independent kinetic evidence.

Direct \(\mathrm{NH}\) entry was distinguished from downstream \(\mathrm{NH_x}\) turnover. For \(s\in\{\mathrm{NH},\mathrm{NH_2},\mathrm{NH_3}\}\), the signed reaction contribution was integrated using the corresponding stoichiometric coefficient. Positive and negative contributions were accumulated separately as gross reaction formation \(P_s\) and destruction \(D_s\):

\[
P_s=\sum_{j\in\mathcal W}\sum_r\max\left(\nu_{sr}R_r(t_j),0\right)\Delta t_j,
\]

\[
D_s=-\sum_{j\in\mathcal W}\sum_r\min\left(\nu_{sr}R_r(t_j),0\right)\Delta t_j.
\]

The ammonia destruction-to-formation ratio was

\[
\frac{D}{P}=\frac{D_{\mathrm{NH_3}}}{P_{\mathrm{NH_3}}},
\]

while CSTR outlet removal was represented separately by \(n_{\mathrm{NH_3}}/\tau\). A near-unity \(D/P\) therefore indicates reaction-network cancellation, not zero outlet flux or thermodynamic equilibrium. This two-layer accounting prevents a direct-\(\mathrm{NH}\) source from being assigned automatically as the cause of a final \(\mathrm{NH_3}\) trend.

Three independent reaction-family disabling calculations were performed at \(E/N=120~\mathrm{Td}\) and \(x_{\mathrm{N_2}}=0.3333\): removal of the four named-\(\mathrm{H_2^*}\)-to-\(\mathrm{NH}\) reactions, removal of the Rydberg-\(\mathrm{H_2^*}\)-to-\(\mathrm{NH}\) reaction, and removal of the two excited-\(\mathrm{N}\)-to-\(\mathrm{NH}\) reactions. Each case was generated from a fresh copy of the validated input, rebuilt independently, and compared with the unmodified continuous-wave CSTR baseline. The purpose was to test whether the ROP-identified family materially sustains \(\mathrm{NH_3}\) within this mechanism. Such disabling is a defined numerical intervention, not proof of an invariant rate-limiting step or an experimental pathway attribution.

Local kinetic controls were determined by symmetric finite perturbation. Selected groups included direct-NH families, two NH3-forming association groups, and two ion-mediated NH3-loss groups. For response y, the central logarithmic coefficient was

\[
C_y=\frac{\ln[y(k=1.10)]-\ln[y(k=1/1.10)]}{\ln(1.10)-\ln(1/1.10)}.
\]

Positive values indicate that increasing the selected group locally increases the response, whereas negative values indicate the opposite. Perturbed calculations were independently required to meet \(P0\). Detailed controls were evaluated at the model productivity maximum, \(E/N=140~\mathrm{Td}\) and \(x_{\mathrm{N_2}}=0.1\), and at the loss-dominated state, \(E/N=240~\mathrm{Td}\) and \(x_{\mathrm{N_2}}=0.9\). A \(30\)-point atlas using five nitrogen fractions and six fields mapped selected controls across the domain. These coefficients are local sensitivity diagnostics, not global uncertainty indices.

Residence-time robustness was tested at four contrasting states: \(E/N=60~\mathrm{Td}\) and \(x_{\mathrm{N_2}}=0.1\); \(E/N=140~\mathrm{Td}\) and \(x_{\mathrm{N_2}}=0.1\); \(E/N=140~\mathrm{Td}\) and \(x_{\mathrm{N_2}}=0.5\); and \(E/N=240~\mathrm{Td}\) and \(x_{\mathrm{N_2}}=0.9\). Residence time was varied over

\[
0.5~\mathrm{ms}\leq\tau\leq30~\mathrm{ms}.
\]

To avoid terminating long-residence-time calculations before sufficient flushing, the numerical horizon was at least \(12\tau\), while retaining the \(P0\) requirement. A complementary deterministic scenario envelope independently multiplied the named-\(\mathrm{H_2^*}\)-to-\(\mathrm{NH}\) group, named-\(\mathrm{H_2^*}\) wall-relaxation group, or Rydberg-\(\mathrm{H_2^*}\)-to-\(\mathrm{NH}\) group by \(1/3\) and \(3\) at the productive and loss-dominated states. These are transparent scenario tests; no probability distributions, correlations, or confidence intervals were assigned to rate coefficients.

### 2.5. Network analysis and Fe(110) boundary audit

Reaction records were represented as flux-weighted species-reaction bipartite hypergraphs. Reactions were retained as nodes, preserving their full sets of co-reactants and products instead of reducing each reaction to pairwise arrows. A nitrogen-lineage projection was used for weighted-betweenness and pathway-hierarchy summaries. Displayed records were checked for \(\mathrm{N}\) and \(\mathrm{H}\) atom balance. Graph metrics were used to organise the active flux network and identify structural alternatives; integrated ROP flux remained the quantitative basis for direct-pathway ranking [15].

A public Fe(110) DFT-microkinetic implementation from Shao and Mesbah [16] was audited as a distinct surface branch and was not merged with the CSTR. Its source contains \(468\) active records, including \(46\) explicit surface reactions. Reproduction of four native field-temperature cases gave a maximum absolute relative difference of \(0.06345\%\) for the \(4000~\mathrm{s}\) \(\mathrm{NH_3}\) number density. This verifies portability at the source boundary, but it does not independently validate the gas-phase direct-\(\mathrm{NH}\) ranking because the two inputs share neutral direct-\(\mathrm{NH}\) reactions.

The CSTR uses E/N, whereas the Fe(110) branch uses a local surface field. Their dimensional relation is

\[
E_{\mathrm{bulk}}=\frac{(E/N)p}{k_{\mathrm{B}}T},\qquad F_s=\beta E_{\mathrm{bulk}},
\]

where \(\beta\) represents geometry, sheath, dielectric, packing, material, and local-field enhancement. It was neither set to unity nor fitted to \(\mathrm{NH_3}\). No gas-surface flux ranking is therefore reported without a common reactor boundary, independently calibrated local field, compatible surface-area and site-density normalisation, and simultaneous gas- and surface-resolved ROP calculation.

For every accepted result, the project archive retains input provenance, run parameters, \(P0\) status, terminal species series, ROP records, and figure-source CSV files. The primary claims are limited to the corrected gas-phase continuous-wave CSTR mechanism. They do not establish pulse-frequency or duty-cycle effects, device energy efficiency, self-consistent discharge behaviour, or universal gas-versus-surface pathway dominance.


## 3. Results

### 3.1. A common P0 gate establishes a comparable continuous-wave CSTR map

Figure 1 establishes the numerical basis of the paper. At the reference condition, the continuous-wave CSTR calculation satisfied the terminal-window P0 criterion while the matched closed-system regression continued to accumulate NH3 over the same finite interval. The terminal CSTR state, NH3 residual, and energy residual all crossed the \(10^{-3}\) P0 criterion. The full grid contained 108 accepted points, including two visible retry cases. Thus, the pathway maps below compare equal, numerically converged terminal windows rather than arbitrary transients with different degrees of accumulation.

The maximum relative-state residual among accepted grid points was \(8.40\times10^{-4}\), and the maximum NH3 residual was \(8.33\times10^{-4}\). The stored window-energy diagnostic is zero within the output precision because the continuous-wave forcing is fixed; it is not evidence of energy closure. Together, these diagnostics give the numerical meaning of “accepted” in every result that follows. They establish a common analysis basis for the continuous-wave CSTR only; they are not evidence of a pulsed-discharge periodic state.

![Figure 1](../analysis/publication_figures_20260903_v2/Figure_1_CSTR_P0_validation.png)

### 3.2. Persistent direct NH entry emerges as a map-level result

Figure 2 defines the mechanism boundary, the CSTR mass-balance boundary, and the direct-NH source-family ledger. Figure 3 then reports the central mapping result. The named H2* family was the largest integrated direct source of NH at every accepted x_N2–E/N condition. Its share ranged from 84.892% to 99.664%. The lowest share occurred at 240 Td and \(x_{\mathrm{N_2}}=0.9\), yet the named-H2* contribution remained larger than each competing direct-NH family. There was consequently no direct-NH pathway crossover within the accepted map.

This conclusion is deliberately narrower than a statement about total NH3 chemistry. It says that the first direct entry into NH remains electronically H2*-dominated in this corrected mechanism. It does not say that the named family always controls NH3 inventory, that it is the only active branch, or that it would remain dominant after introducing surfaces, spatial gradients, pulsing, or different excited-state chemistry. The apparent absence of crossover is nevertheless robust to the composition and E/N variation represented here.

![Figure 2](../analysis/publication_figures_20260903_v2/Figure_2_CSTR_mechanism_map.png)

![Figure 3](../analysis/publication_figures_20260903_v2/Figure_3_dense_EN_N2_CSTR_maps.png)

### 3.3. Family disabling establishes a strong causal role for the named H2* route

At 120 Td and \(x_{\mathrm{N_2}}=0.3333\), the causal diagnostic shown in Figure 2C distinguishes a visible ROP branch from a productive causal gate. Disabling the named electronic-H2* family decreased NH3 by 83.31% relative to the baseline. In contrast, disabling the tested Rydberg alternative changed NH3 by −0.60%, and disabling the tested excited-N alternative changed NH3 by −0.32%. These small signed responses do not demonstrate that those branches are chemically irrelevant; they show that, at this state and under this implementation, they do not sustain the bulk of the NH3 inventory. The result is consistent with the terminal ROP ranking and supplies a causal test that a pathway share alone cannot provide.

**Table 2. Isolated diagnostic perturbations at 120 Td and \(x_{\mathrm{N_2}}=0.3333\).** A disabled branch is a diagnostic counterfactual, not a replacement mechanism.

| Diagnostic branch | Final NH3 (cm\(^{-3}\)) | Change relative to baseline | Interpretation |
|---|---:|---:|---|
| Baseline | \(1.534\times10^{14}\) | — | Unmodified reference mechanism |
| Disable Rydberg-H2* to NH | \(1.524\times10^{14}\) | −0.60% | Secondary direct source at this reference state |
| Disable named H2* to NH | \(2.559\times10^{13}\) | −83.31% | Strong causal contribution to the NH3 network |
| Disable N(2D)/N(2P) to NH | \(1.529\times10^{14}\) | −0.32% | Minor direct contribution at this reference state |

![Figure 4](../analysis/publication_figures_20260903_v2/Figure_4_productivity_window.png)

### 3.4. NH3 exhibits an intermediate-field productive window despite persistent NH-entry dominance

The product response differs sharply from the direct-NH pathway map. The final NH3 inventory reaches its maximum of \(2.557\times10^{14}\) cm\(^{-3}\) at 140 Td and \(x_{\mathrm{N_2}}=0.1\), corresponding to an outlet productivity of \(2.557\times10^{16}\) cm\(^{-3}\) s\(^{-1}\). Both lower-field nitrogen-rich conditions and the strongest-field conditions yield less product, despite the named-H2* family retaining a large direct-NH share. Figure 4 shows the resulting intermediate-field productive window.

The comparison reveals a key separation of mechanism layers. A high named-H2* fraction identifies how NH first enters the NHx network. It does not establish whether NH is efficiently hydrogenated to NH2 and NH3, whether product loss remains modest, or whether a CSTR residence time permits net accumulation. Product-level interpretation must therefore retain the downstream turnover network.

### 3.5. Temperature and electron-density endpoints do not expose a low-vibrational-H2 crossover

Figure 5 reports the endpoint robustness calculations. They do not identify a condition at which a low-vibrational-H2 route overtakes the named electronic-H2* direct-NH family. This does not establish global robustness to all electron densities or temperatures; the endpoint exercise samples only a limited subset of possible model boundaries. Its value is more modest and more precise: the observed lack of direct-NH crossover is not an artefact of the particular nominal endpoint selected for the primary CSTR map.

![Figure 5](../analysis/publication_figures_20260903_v2/Figure_5_temperature_and_electron_robustness.png)

### 3.6. Flux-weighted hypergraphs preserve the chemical arity of the active network

Figure 6 visualizes three representative terminal states: a low-field point (20 Td, \(x_{\mathrm{N_2}}=0.1\)), the productivity maximum (140 Td, \(x_{\mathrm{N_2}}=0.1\)), and a high-field nitrogen-rich point (240 Td, \(x_{\mathrm{N_2}}=0.9\)). The flux-weighted bipartite hypergraphs retain reaction nodes, reactant co-requirements, and product branching. This representation avoids the misleading implication that a three-body association or ion–molecule process can be reduced to an ordinary pairwise arrow.

The named electronic-H2* entry branch is prominent in every representative state. The high-field graph develops stronger association and secondary-excitation structure upstream and stronger ion-mediated activity downstream, but these changes do not overturn the integrated direct-NH hierarchy. The graph therefore supports a mechanistic interpretation in which network texture changes with field even though the leading direct NH-entry family persists.

![Figure 6](../analysis/publication_figures_20260903_v2/Figure_6_flux_weighted_reaction_hypergraphs.png)

### 3.7. Topological diagnostics support, but do not replace, the ROP ranking

Figure 7 combines the dense map of named-H2* flux divided by the largest secondary flux with reaction- and family-resolved rankings at the representative conditions. The dominance margin remains above unity throughout the accepted map, quantitatively excluding a direct-NH crossover by the chosen criterion. The parsed hypergraph records number 434, 468, and 468 at the low-field, productive-window, and high-field representative states, respectively. All three record sets have zero N- and H-atom imbalance.

The nitrogen-lineage weighted-betweenness projection identifies intermediates that structurally connect the active nitrogen network. This is useful for deciding where to inspect flux and sensitivity data, but it does not provide a state-dependent rate by itself. The agreement between topology and ROP is therefore treated as convergence of two forms of evidence rather than as independent proof of kinetics.

![Figure 7](../analysis/publication_figures_20260903_v2/Figure_7_graph_theoretical_path_hierarchy.png)

### 3.8. NH3 accumulation is controlled by formation–destruction competition

Figure 8 moves downstream from NH entry to NH3 turnover. At the productivity maximum, the terminal NH3 formation and destruction fluxes are \(3.562\times10^{16}\) and \(1.024\times10^{16}\) cm\(^{-3}\) s\(^{-1}\), respectively, giving \(D/P=0.287\). This imbalance permits a substantial positive outlet productivity. At 240 Td and \(x_{\mathrm{N_2}}=0.9\), formation and destruction are both large but nearly cancel: \(1.459\times10^{15}\) and \(1.458\times10^{15}\) cm\(^{-3}\) s\(^{-1}\), respectively, corresponding to \(D/P=0.999436\) and an outlet productivity of only \(8.308\times10^{11}\) cm\(^{-3}\) s\(^{-1}\).

The leading high-field NH3 loss is \(\mathrm{H^+ + NH_3 \rightarrow NH_3^+ + H}\), followed by \(\mathrm{NH_3^+ + NH_3 \rightarrow NH_4^+ + NH_2}\). The high field therefore does not simply “destroy ammonia”; it shifts the system into a high-throughput, nearly cancelling production–loss state. This distinction is essential when interpreting a low outlet inventory.

![Figure 8](../analysis/publication_figures_20260903_v2/Figure_8_NHx_source_sink_turnover.png)

### 3.9. Material-balance descriptors locate the productive region without invoking energy efficiency

Figure 9 maps the CSTR outlet productivity, N2 conversion, nitrogen-atom utilization to NH3, and the productivity–conversion frontier. The productivity maximum coincides with the 140 Td, \(x_{\mathrm{N_2}}=0.1\) inventory maximum. In contrast, maximum N2 conversion is 0.2527% at 240 Td and \(x_{\mathrm{N_2}}=0.1\), whereas nitrogen utilization peaks at 7.567% at 120 Td and \(x_{\mathrm{N_2}}=0.1\). These distinct extrema demonstrate that conversion, utilization, and product outlet cannot be treated as interchangeable objectives.

The figure intentionally reports only material-balance descriptors. The fixed-electron-density model does not establish absorbed power, electrical waveform, or device-level energy deposition, so assigning an energy-efficiency ranking would exceed the evidence.

![Figure 9](../analysis/publication_figures_20260903_v2/Figure_9_CSTR_productivity_and_conversion.png)

### 3.10. Product-level local control shifts downstream of the direct NH gate

At the productivity optimum, local perturbation analysis identifies the association \(\mathrm{NH+H_2+M\rightarrow NH_3+M}\) as a strong positive control, with a normalized coefficient of 1.035. At high field, the control structure shifts. The three-body reaction \(\mathrm{H+NH_2+M\rightarrow NH_3+M}\) has a positive coefficient of 0.974, while proton-driven NH3 ionization contributes a strong negative coefficient of −0.856. Thus, the rate levers for product inventory can lie downstream of the persistent direct-NH source gate.

![Figure 10](../analysis/publication_figures_20260903_v2/Figure_10_local_kinetic_control.png)

### 3.11. The kinetic-control atlas separates formation-controlled and loss-controlled regions

Figure 11 extends the local-control panel to 30 targeted operating states and 90 sensitivity runs. Of the 30 leading positive controls, 22 are NH-association controls and seven are downstream hydrogenation controls; a proton-related contribution appears only in a weak low-field case. In contrast, proton-driven ionization is the leading negative control at every 240 Td targeted condition. The field therefore introduces a coherent control transition: the principal positive levers remain in NHx formation, but the principal negative lever becomes ammonia ionization and consequent ionic turnover.

![Figure 11](../analysis/publication_figures_20260903_v2/Figure_11_kinetic_control_switching_atlas.png)

### 3.12. A two-layer regime map reconciles persistent NH entry with variable product outcome

Figure 12 brings the two mechanism layers together. The first layer classifies direct NH entry by the dominance of the named H2* family. The second layer classifies NH3 turnover by the balance between formation and destruction and by the leading local controls. The combination explains why the model can retain an electronically H2*-dominated NH entry over the full map while displaying an intermediate-field productive region and a high-field loss-dominated region.

The two-layer map is also a guard against over-interpretation. It prevents a direct pathway label from being used as a surrogate for overall ammonia performance, and it identifies which physical additions—such as self-consistent electron kinetics, pulse-resolved afterglow chemistry, or surfaces—would have to modify either the entry layer or the turnover layer to alter the conclusion.

![Figure 12](../analysis/publication_figures_20260903_v2/Figure_12_two_layer_mechanistic_regime_map.png)

### 3.13. The direct-NH dominance margin quantifies the absence of crossover

Figure 13 strengthens the all-grid source-share result by reporting the direct dominance margin rather than a percentage alone. The numerator is the integrated named-H2* direct-NH source and the denominator is the largest integrated secondary direct-NH source at each accepted state. This comparison is stricter than checking whether the named family exceeds the sum of selected minor sources, because it tests the only competitor capable of defining a direct-pathway crossover.

The margin remains greater than one everywhere in the accepted map. The closest approach occurs at the high-field, nitrogen-rich boundary, but it does not reach parity. The conclusion is therefore not merely that the named family usually dominates; it dominates with a nonzero margin over the strongest direct competitor throughout the scanned CSTR domain.

![Figure 13](../analysis/publication_figures_20260903_v2/Figure_13_direct_NH_dominance_margin.png)

### 3.14. The NH-entry-to-NH3 cascade exposes the location of product decoupling

Figure 14 traces the terminal flux cascade from direct NH entry through NH2 and NH3 formation to NH3 destruction. At the productive state, the NH-entry flux is effectively transmitted through the hydrogenation sequence and the final NH3 destruction branch remains subdominant. At high field, direct NH entry remains substantial, but the transmitted product signal is diminished by downstream competition and by the ion-mediated destruction loop.

This cascade makes the central decoupling visible. The electronically excited H2* branch sets the main upstream entrance to NH. The NH3 inventory is decided later, at the point where formation, recycling, and ionic loss divide the downstream flux. The resulting mechanistic interpretation is stronger than a single reaction-path diagram because every branch is weighted by the terminal ROP accumulated at its own operating state.

![Figure 14](../analysis/publication_figures_20260903_v2/Figure_14_NH_entry_to_NH3_turnover_cascade.png)

### 3.15. Product-turnover control phase space identifies formation and loss regimes

Figure 15 overlays product-turnover descriptors with the leading local controls. It separates a formation-favoured region around the intermediate-field productive window from high-field states in which the destruction fraction approaches one and proton-driven loss becomes the salient negative lever. The graphic is particularly useful because it does not conflate a high gross formation flux with a favourable outlet condition. A state can have vigorous ammonia formation and still be a poor producer if nearly all formed NH3 is turned over before exit.

The phase-space view also provides a testable design hypothesis. Efforts to improve the intermediate-field region should prioritize the downstream association and hydrogenation reactions carrying positive control. Efforts to rescue the high-field region should instead suppress the ion-mediated NH3-loss loop or alter the plasma state that supplies it. Whether those interventions are physically realizable requires a self-consistent discharge and, for catalytic designs, a surface-resolved model.

![Figure 15](../analysis/publication_figures_20260903_v2/Figure_15_product_turnover_control_phase_space.png)

### 3.16. The central conclusion is jointly numerical, causal, and topological

Figure 16 consolidates the evidence hierarchy. Numerical reliability is supplied by the P0 acceptance gate and the explicit retry record. Mechanistic ranking is supplied by the integrated direct-NH map and the dominance margin. Causal relevance is supplied by the family-disabling test. Product-level interpretation is supplied by the formation–destruction turnover and the local-control atlas. The hypergraph and graph-theoretical analyses give a transparent network representation and atom-balance audit. Figure 16 is an evidence-synthesis figure, not an additional independent observation.

No individual panel would justify a broad pathway claim. Together, the panels support a more disciplined statement: within the corrected Hong mechanism and the stated CSTR conditions, the named electronic-H2* family is a persistent dominant direct NH-entry route, while NH3 performance is controlled by a separate downstream turnover layer that changes across the operating map.

![Figure 16](../analysis/publication_figures_20260903_v2/Figure_16_numerical_mechanistic_robustness_overview.png)

### 3.17. A targeted residence-time screen preserves direct-NH dominance while changing turnover

Figure 17 tests whether the two-layer interpretation is an artefact of the nominal 10 ms CSTR residence time. Four contrasting continuous-wave states were evaluated from 0.5 to 30 ms: a low-field state (60 Td, \(x_{\mathrm{N_2}}=0.1\)), the productivity maximum (140 Td, \(x_{\mathrm{N_2}}=0.1\)), an intermediate-composition state (140 Td, \(x_{\mathrm{N_2}}=0.5\)), and the high-field loss-dominated state (240 Td, \(x_{\mathrm{N_2}}=0.9\)). All 20 calculations passed P0 after a run horizon of at least 12 residence times.

The direct-NH hierarchy persists throughout this targeted transport screen. The named-H2* source share ranges from 76.832% to 99.596%, and its dominance margin over the largest secondary direct source ranges from 4.885 to 383.46. The smallest margin occurs at 240 Td, \(x_{\mathrm{N_2}}=0.9\), and 30 ms, where the high-field loss loop is most fully expressed; it nevertheless remains above the crossover criterion of one. Residence time therefore changes how much product can accumulate and how closely formation and destruction cancel, but does not overturn the leading direct NH-entry family in this panel.

The result reinforces the two-layer interpretation rather than replacing the full composition–field map. At the high-field state, productivity falls from \(2.127\times10^{15}\) to \(9.256\times10^{10}\) cm\(^{-3}\) s\(^{-1}\) as residence time increases from 0.5 to 30 ms, while \(D/P\) rises from 0.428 to 0.999953. The transport extension is deliberately targeted: it supports robustness to residence time within the stated mechanism, but it is not a full reactor response surface or a kinetic uncertainty quantification.

![Figure 17](../analysis/publication_figures_20260903_v2/Figure_17_residence_time_robustness.png)

### 3.18. A deterministic H2* reaction-scenario envelope does not induce crossover

Figure 18 asks a different robustness question from Figure 17. Instead of changing transport, it perturbs individually the reaction groups most directly implicated in the named-H2* interpretation. At the productivity maximum and the high-field loss-dominated state, the four named H2* to NH channels, the four named-H2* wall-relaxation channels, and the Rydberg-H2* to NH channel were each independently scaled by \(1/3\) and 3. All 14 scenario calculations passed P0.

The direct-NH hierarchy remained intact in every scenario. The most adverse case is the high-field state with the named-H2* to NH group reduced to one third: its named-H2* share is 84.877% and its dominance margin is 10.724. The productive state remains farther from crossover, with margins of 60.259–65.044 across the named-H2* to NH scenarios. The Rydberg perturbations produce only small changes in the source share and margin over this range.

This result strengthens, but does not universalize, the mechanistic claim. It demonstrates that the two selected reference states do not cross over under these explicit independent rate-group scenarios. It does not provide a probability that crossover is absent, quantify uncertainty in the BOLSIG+ cross sections, or cover coupled changes in several rate coefficients. These distinctions are retained in the caption and discussion so Figure 18 is not misread as formal statistical uncertainty analysis.

![Figure 18](../analysis/publication_figures_20260903_v2/Figure_18_H2star_scenario_envelope.png)

**Table 3. Claim–evidence–boundary map for the central two-layer mechanism.**

| Claim | Primary evidence | Counterfactual or cross-check | Boundary that remains |
|---|---|---|---|
| The accepted continuous-wave map is numerically comparable | Figure 1; P0 residuals for all 108 cases | Matched closed-system regression fails the same gate | Numerical convergence is not experimental validation |
| Named H2* is the largest direct NH-entry family | Figures 3 and 13; terminal-window integrated ROP | Direct dominance margin is positive at every accepted state | Fixed reconstructed mechanism, 300 K, 10 ms, imposed electron density |
| The named H2* family is causally relevant at the reference state | Figure 2C and Table 2 | Family disabling decreases NH3 by 83.31%; tested alternatives are sub-percent | Local counterfactual at 120 Td and x_N2 = 0.3333, not a global rate-limiting-step proof |
| NH3 outcome is governed downstream of direct NH entry | Figures 4, 8–12, 14, and 15 | Source–sink balance and local controls reveal formation- and loss-dominated regimes | No self-consistent power, pulse waveform, surface chemistry, or experimental calibration |
| Targeted transport variation does not induce a direct-NH crossover | Figure 17; 20 P0-accepted cases from 0.5–30 ms | Minimum named-H2* margin is 4.885 at the high-field, long-residence-time endpoint | Four representative states only; not a full transport response surface |
| Selected H2* reaction scenarios do not induce a direct-NH crossover | Figure 18; 14 P0-accepted 1/3–3× scenarios | Minimum margin is 10.724 after reducing named-H2* to NH at the high-field reference state | Independent deterministic scenarios only; not statistical UQ or a cross-section uncertainty envelope |
| Direct-NH family distribution remains concentrated | Figure 26; all 108 accepted states | \(H/\ln(6)=0.013812\)–0.328297 and \(N_{\mathrm{eff}}=1.025057\)–1.800795; minimum named-H2*/runner-up ratio = 10.054666 | Descriptive transform of the same integrated ROP shares used in Figures 3 and 13; not independent kinetic evidence or UQ |
| Network topology is chemically consistent with the flux result | Figures 6, 7, and 16 | Atom-balance audit and agreement with ROP ranking | Graph metrics organize evidence; they do not independently establish kinetic causality |

### 3.19. Contrasting gas-CSTR states define a prospective discrimination atlas

Figure 19 converts the accepted gas-CSTR evidence into a falsifiable state-selection ledger. It retains a low-field entry-limited state (20 Td, \(x_{\mathrm{N_2}}=0.1\)), the productive reference state (140 Td, \(x_{\mathrm{N_2}}=0.1\)), and a high-field nitrogen-rich loss-stressed state (240 Td, \(x_{\mathrm{N_2}}=0.9\)). The selection is not based on NH3 alone. It jointly preserves the direct-NH hierarchy, the dominance margin, and the contrast in downstream turnover control. Thus, the three states provide intentionally different tests of an upstream entry response and a downstream formation–loss response.

The atlas does not report independent validation. Rather, it states before an experiment or cross-model fit which observables would be needed to challenge the gas-only interpretation: a calibrated NH3 outlet measurement, time-resolved electrical observables, state-sensitive NH/NHx diagnostics, feed-nitrogen closure, and a technically validated constraint on the excited-state or electron boundary. A disagreement in a single NH3 trace would not diagnose the route, whereas a reproducible multi-observable departure from the predicted ordering would falsify the corresponding bounded model explanation.

![Figure 19](../analysis/publication_figures_20260903_v2/Figure_19_graph_informed_model_discrimination_atlas.png)

### 3.20. A source-validated Fe(110) branch broadens the structural alternatives but not the gas-phase evidence

Figure 20 distinguishes mechanism-family correspondence from numerical comparability. The source-validated Fe(110) branch retains a Hong-derived gas backbone while adding 46 explicit surface records. Its ten neutral \(N/N(2D)/N(2P)+H_2^*\rightarrow H+NH\) entries exactly overlap the corresponding local gas-phase reactions. The Fe(110) model is therefore not an independent validation of the named-H2* conclusion. Its scientific value here is different: it supplies a documented set of adsorption, dissociation, Eley–Rideal and surface-hydrogenation alternatives that a gas-only CSTR cannot represent.

The author-boundary reproduction and porting audit in Figure 21 confirm that the external source can be reproduced without modifying its reaction network, but they do not authorize cross-boundary flux comparison. The gas model imposes \(E/N\) and a CSTR inflow/outflow operator, whereas the Fe(110) model uses a local field and source-native surface geometry and site normalization. A reaction label in common is consequently insufficient to establish a common kinetic or reactor state.

![Figure 20](../analysis/publication_figures_20260903_v2/fe110_and_design/figure20/Figure_20_mechanism_ensemble_readiness.png)

![Figure 21](../analysis/publication_figures_20260903_v2/fe110_and_design/figure21/Figure_21_external_model_scope_and_porting_gate.png)

### 3.21. Dimensional field auditing prevents an invalid gas–surface pathway comparison

Figure 22 makes the necessary field bridge explicit. At 1 atm and 300 K, the 20, 140 and 240 Td gas states correspond to bulk fields of 0.000049, 0.000342 and 0.000587 V Å\(^{-1}\), respectively. Reaching the 0.06–0.11 V Å\(^{-1}\) reference local fields of the Fe(110) model would require \(\beta\) values spanning 102–2248 across these anchors. This is a dimensional requirement, not an inferred enhancement factor for the present device.

The consequence is substantive rather than semantic. A surface rate-of-progress is physically interpretable against the gas-CSTR result only after an independent electrostatic calculation or calibrated field diagnostic constrains \(\beta(E/N,\mathrm{geometry},\mathrm{material})\), and after pressure, temperature, mixture, electron boundary, transport, area-to-volume ratio, roughness, site density, inlet and residence time are shared. Until then, the allowed result is a comparison protocol and a pre-registered set of gas/surface first-N–H fractions and independent observables, not a numerical claim that one branch dominates the other.

![Figure 22](../analysis/publication_figures_20260903_v2/fe110_and_design/figure22/Figure_22_boundary_harmonization_protocol.png)

### 3.22. Surface first-N–H alternatives and graph symmetries are structural constraints, not kinetic rankings

The frozen Fe(110) input encodes five reactions that create `NHSurf` by joining N-bearing and H-bearing precursors (Figure 23). Four are Eley–Rideal alternatives: \(N+H\mathrm{Surf}\), \(N(2D)+H\mathrm{Surf}\), \(N(2P)+H\mathrm{Surf}\), and \(H+N\mathrm{Surf}\). The fifth is the Langmuir–Hinshelwood step \(N\mathrm{Surf}+H\mathrm{Surf}\rightarrow NH\mathrm{Surf}+\mathrm{Surf}\). The same source contains 13 \(N_2\)-state dissociative-adsorption records, 8 \(H_2\)-state dissociative-adsorption records and multiple downstream `NHSurf`–`NH3` release routes. These record counts establish that the surface first-N–H decision is structurally multibranched; they do not assign a probability or flux to any branch.

Figure 24 formalizes the same distinction through automorphism orbits of the unweighted directed species–reaction graph. Twelve excited or vibrational \(N_2\) labels, four named electronic \(H_2\) labels, three low-vibrational \(H_2\) labels, and the three active-N labels \(N\), \(N(2D)\), and \(N(2P)\) each admit an explicitly verified topology-preserving cyclic permutation, including their incident reaction nodes. The resulting graph symmetry is broken before any kinetic ranking by state-specific DFT energetics, cross sections, populations, wall loss and local field. No state labels were therefore merged and no degeneracy multiplier was introduced. The group-theoretic result identifies the conditions that a valid reduced mechanism would have to demonstrate, rather than supplying a shortcut to reaction-rate equivalence.

![Figure 23](../analysis/publication_figures_20260903_v2/fe110_and_design/figure23/Figure_23_Fe110_surface_NH_entry_graph.png)

![Figure 24](../analysis/publication_figures_20260903_v2/fe110_and_design/figure24/Figure_24_Fe110_surface_symmetry_orbits.png)

### 3.23. A prospective factorial design makes the gas–surface hypothesis experimentally falsifiable

Figure 25 translates the structural audit into a pre-registered validation design. A blocked \(2^3\) factorial crosses a Fe(110)-representative catalyst with a matched inert surface, independently calibrated low and high local fields, and \(^{14}\mathrm{N_2}/\mathrm{H_2}\) versus \(^{15}\mathrm{N_2}/\mathrm{H_2}\) feeds. The displayed schedule has three independent reactor/catalyst-reconditioning blocks and eight randomized arms per block. It is a reproducible pilot and variance-estimation layout, not a formal sample-size justification.

The independent unit is a complete block/run, not an individual point in a time trace. The proposed measurement panel consequently combines waveform-resolved V–I–Q, NH/NHx and active-N diagnostics, isotope-resolved product with blank and memory controls, and surface characterization. A Fe-versus-inert contrast can test a material-conditioned response after electrical matching; an isotope contrast can test feed-nitrogen provenance; and a calibrated local-field contrast can test field response. None alone identifies an Eley–Rideal or Langmuir–Hinshelwood flux. The design reserves pathway ranking for a subsequent hybrid model that satisfies the Figure 22 comparison gates.

![Figure 25](../analysis/publication_figures_20260903_v2/fe110_and_design/figure25/Figure_25_fe110_pathway_discrimination_design.png)

### 3.24. Direct-NH pathway concentration remains low even at the most diverse accepted state

Figure 26 converts the same six direct-NH family shares used for the primary source map into pathway entropy and an effective family count. The normalized entropy is only 0.013812–0.328297 across all 108 accepted points, corresponding to \(N_{\mathrm{eff}}=1.025057\)–1.800795 rather than a broadly distributed six-family mixture. The most diverse point is 240 Td and \(x_{\mathrm{N_2}}=0.9\), where the named-H2* share is still 84.892324% and association is the runner-up. The smallest named-H2*/runner-up ratio anywhere in the map is 10.054666.

The identity of the runner-up changes systematically: N(2P) + H2 is second at 59 points, Rydberg-H2* at 46 points, and association at three high-field nitrogen-rich points. This supports a useful refinement of the no-crossover statement: the secondary competitor changes with operating condition, but the direct-NH distribution never approaches an equal-share mixture in the accepted map. Because Figure 26 is calculated from the same terminal integrated ROP family totals as Figures 3 and 13, it is a compact quantitative description of the existing evidence, not an independent validation layer or an uncertainty analysis.

![Figure 26](../analysis/publication_figures_20260903_v2/Figure_26_direct_nh_pathway_concentration.png)

## 4. Discussion

The first major result is a conditional absence-of-crossover finding. In a field where individual studies can emphasize electronically excited species, vibrationally excited species, atomic radicals, or surface routes, it is scientifically tempting to frame the literature as contradictory. The present result illustrates a more useful interpretation: different pathway reports often refer to different reaction-network boundaries. The reconstructed Hong mechanism contains a particular representation of electronically excited H2 and a corrected representation of H2(v)-mediated NH formation. With this representation, imposed electron density, 300 K CSTR boundary, and 10 ms residence time, the electronic-H2* family remains the largest direct source of NH from 20 to 240 Td and from 0.1 to 0.9 N2 fraction. “No crossover” therefore means that no competitor exceeds the named-H2* direct flux in the accepted terminal windows of this map; it does not mean that crossover is impossible in another mechanism or reactor.

The claim does not conflict with studies in which vibrational N2 assists surface activation [4,5] or in which H2(v) is important for downstream NH consumption [13]. Those studies ask different questions and include different physical ingredients. Catalyst-containing systems can introduce adsorption, site competition, Eley–Rideal reactions, and altered electric fields [8–11]. Higher-pressure systems can change collision frequencies and the relative survival of excited states [13]. A comparison without preserving these boundaries would convert a useful mechanism calculation into an unjustified universal assertion.

The second major result is the separation of NH entry from NH3 productivity. The named H2* route remains prominent at both the productive intermediate-field state and the high-field nitrogen-rich state. Yet the latter state has a very small outlet productivity because its gross NH3 formation is nearly cancelled by destruction. This is an important practical message for kinetic interpretation: a reaction family that creates the first N–H bond is not necessarily the rate lever that maximizes final ammonia. The local-control atlas captures this distinction. At favourable conditions, positive controls lie in association and downstream hydrogenation. At high field, proton-driven NH3 ionization becomes the strongest negative control.

The third result concerns evidence architecture. Terminal-window ROP integration is preferable to an instantaneous-rate snapshot when intermediate pools can lag the active chemistry. The P0 gate prevents terminal comparisons from being contaminated by unequal transient accumulation. Family disabling tests add a causal stress test. The pathway-concentration transform shows that the no-crossover result is not merely a ranking artifact created by nearly equal source families: the same ROP totals correspond to an effective direct-NH family count below 1.81 everywhere. Graph and hypergraph analysis add visual accountability: they show whether a purported route is embedded in a chemically plausible, atom-balanced reaction context and whether the same topology persists across representative states. This combination follows the principle that graph theory can expose network connectivity and candidate intervention points, whereas kinetic fluxes retain responsibility for quantitative pathway attribution [15]. Formal group-theoretical equivalence is not imposed here because the electronic H2 labels in the mechanism have different energetics, cross sections, and loss channels; treating them as symmetry-equivalent would obscure physically meaningful distinctions.

The cross-mechanism audit extends this evidence architecture without blurring its layers. The reproduced Fe(110) DFT–microkinetic branch demonstrates that an explicit catalytic surface contains qualitatively different first-N–H alternatives, including both Eley–Rideal and Langmuir–Hinshelwood steps [16]. It does not provide an experimental or independent gas-phase validation, because part of its gas backbone overlaps the Hong direct-NH records. More importantly, the field variables are not interchangeable: a bulk reduced electric field cannot be relabelled as a catalyst-local field. The required \(F_s=\beta E_{\mathrm{bulk}}\) bridge and the recovered source-native area-to-volume and site parameters convert a vague limitation into a testable model-interface problem. This is a stronger statement than simply noting that surfaces are absent: it specifies the physical information that must be supplied before a surface pathway can either preserve or overturn the gas-phase ranking.

Several limitations define the next research steps. First, electron density was fixed rather than calculated with the electric field, chemistry, and waveform. This makes the map a controlled chemical-kinetic study, not a predictive plasma-device simulation. Second, the forcing is continuous-wave; pulse frequency, duty cycle, on-time E/N, and afterglow evolution were not independently varied. A preliminary pulse branch did not pass the same P0 acceptance gate and is therefore intentionally excluded from the present results. Third, the excited-state aggregation, imported BOLSIG+ cross-section mapping, and electronic-state loss treatment are mechanism choices. The present endpoint checks do not constitute a kinetic uncertainty quantification of those choices. Fourth, although a Fe(110) surface branch is now source-validated and structurally audited, it has not been boundary-harmonized with the gas CSTR; the paper therefore cannot decide whether a catalyst would preserve, divert, or bypass the gas-phase electronic-H2* NH-entry gate. Fifth, the study is not experimentally calibrated. Recent experiment–model work demonstrates the value of using measured plasma conditions and product diagnostics to constrain zero-dimensional mechanisms [12,13].

These limitations also define a concrete falsification program. A targeted residence-time screen now shows that the direct-NH margin remains positive from 0.5 to 30 ms at four contrasting states, and independent 1/3–3× scenarios for three selected H2*-relevant reaction groups likewise preserve the margin at two reference states. A full transport response surface and a formal uncertainty analysis for electron-impact excitation, N + H2* to NH conversion, excited-state loss, and correlated cross-section choices are still needed; the deterministic scenarios in Figure 18 are not probability-weighted uncertainty bounds. A pulse-resolved model with self-consistent electron kinetics should be admitted only after it passes the same numerical acceptance gate; it can then test whether the dominance margin collapses during the afterglow or changes with duty cycle. Surface-resolved simulations should test whether adsorption and Eley–Rideal hydrogenation introduce a competing N–H entry branch. Finally, time-resolved diagnostics of NHx-related species and product concentration would determine whether the predicted formation–loss transition can be observed experimentally. The present paper supplies the baseline against which each extension can be compared rather than claiming to have completed those extensions.

## 5. Conclusions

Using a literature-traceable reconstruction of the corrected Hong N2/H2 gas-phase mechanism in a continuous-wave plasma CSTR, this work maps direct NH formation and NH3 turnover over 108 numerically accepted x_N2–E/N conditions. The named electronically excited H2* family is the largest integrated direct NH source at every accepted point, contributing 84.892–99.664% of direct NH formation. A direct-NH crossover was not observed, and the flux margin relative to the largest secondary source remains positive across the entire accepted map. In a separate 20-case, 0.5–30 ms targeted residence-time screen, the named-H2* share remains 76.832–99.596% and the minimum direct dominance margin is 4.885. Across 14 accepted independent 1/3–3× H2*-relevant rate-group scenarios at two reference states, the minimum margin is 10.724. At a representative 120 Td state, disabling the named H2* family reduces NH3 by 83.31%, establishing that its high ROP share is causally important within this local model counterfactual.

The final NH3 response is governed by a distinct downstream layer. The product inventory and outlet productivity maximize at 140 Td and \(x_{\mathrm{N_2}}=0.1\), whereas high-field nitrogen-rich conditions exhibit nearly complete formation–destruction cancellation. Proton-driven ionization is the principal high-field negative control. Figures 1–18 jointly establish numerical convergence, pathway hierarchy, causal relevance, reaction-network context, turnover competition, and kinetic-control switching. Figure 26 quantitatively describes the concentration of that same direct-NH family distribution, while Figures 19–25 provide the source-validated cross-mechanism audit, dimensional comparison gate, surface-topology and symmetry analyses, and a prospective falsification design. These later figures expand the explanatory and validation scope without being recast as independent gas–surface performance results.

The central conclusion is intentionally bounded: electronically excited H2* is the persistent dominant *direct NH-entry* family within the reconstructed corrected Hong mechanism and the stated CSTR boundary, not necessarily in all plasma or plasma-catalytic ammonia reactors. The 0.5–30 ms targeted transport screen and the selected deterministic rate-group scenarios strengthen this statement without turning it into a statistical uncertainty claim. The Fe(110) audit shows that a catalyst can encode distinct surface first-N–H alternatives, but it also establishes why a gas–surface ranking cannot be made from unmatched \(E/N\) and \(F_s\) models. The two-layer separation of NH entry from NH3 turnover provides a specific, falsifiable framework for formal kinetic uncertainty, pulse-resolved, surface-resolved, and experimentally constrained tests.

## Figure captions

**Figure 1. Terminal-window numerical-convergence validation for the continuous-wave plasma CSTR.** The P0 gate compares normalized state, NH3, and energy residuals for the CSTR/CW calculation with a matched closed-system regression. All 108 accepted grid points satisfy the terminal three-window \(10^{-3}\) residual requirement. In continuous-wave operation, the repeated window is a numerical analysis unit, not a physical pulse period; retry cases are retained explicitly.

**Figure 2. Model provenance, boundary, and selected causal diagnostics.** A literature-traceable reconstruction of the corrected Hong gas-phase mechanism is embedded in a continuously fed and exhausted zero-dimensional reactor. Panels A and B summarize the sparse CSTR evidence, panel C gives the reaction-family disabling calculation, and panel D gives the Rydberg-loss sensitivity. The diagram distinguishes imposed boundary variables from outputs and identifies the named electronic-H2* route investigated in this work; the dense Figure 3 map supplies the global N2–E/N conclusion.

**Figure 3. Dense direct-NH pathway map over composition and reduced electric field.** Integrated terminal-window direct-NH source shares are reported for \(x_{\mathrm{N_2}}=0.1\)–0.9 and E/N = 20–240 Td. The named electronic-H2* family is largest at every accepted point.

**Figure 4. Intermediate-field window for NH3 inventory and outlet productivity.** Product-level response is non-monotonic with E/N even though the named electronic-H2* direct-NH share remains high. The inventory and productivity maximum occurs at 140 Td and \(x_{\mathrm{N_2}}=0.1\).

**Figure 5. Endpoint robustness of direct-NH pathway ranking.** Temperature and imposed-electron-density endpoint tests do not reveal a low-vibrational-H2 direct-NH family that exceeds the named electronic-H2* family. These calculations are robustness checks, not a full uncertainty analysis.

**Figure 6. Flux-weighted reaction hypergraphs at three representative CSTR states.** Species and reaction nodes preserve the arity of multireactant elementary processes. Edge widths are proportional to terminal three-window integrated ROP. Low-field, productivity-maximizing, and high-field nitrogen-rich states show different network texture while retaining the leading named-H2* NH-entry branch.

**Figure 7. ROP hierarchy and graph-theoretical pathway context.** The dominance ratio, family- and reaction-resolved fluxes, and nitrogen-lineage weighted-betweenness diagnostics are presented for representative states. All parsed reaction-record sets are N- and H-atom balanced.

**Figure 8. NH3 source–sink turnover distinguishes productive and cancelling regimes.** Gross NH3 formation and destruction are compared at the productivity maximum and at the high-field nitrogen-rich state. High-field proton-driven ionization yields almost complete formation–loss cancellation.

**Figure 9. CSTR material-balance performance descriptors.** Outlet productivity, N2 conversion, N utilization to NH3, and the productivity–conversion frontier identify distinct extrema. No energy-efficiency claim is made because electron density is imposed and absorbed power is not self-consistently modeled.

**Figure 10. Local kinetic control of NH3.** Positive formation controls and negative loss controls are evaluated at representative states. The leading controls move downstream of the persistent direct-NH entry gate.

**Figure 11. Kinetic-control switching atlas.** The targeted 30-state sensitivity map resolves the transition from association and hydrogenation controls to strong high-field proton-driven NH3-loss control.

**Figure 12. Two-layer mechanistic regime map.** Direct NH entry and NH3 turnover are classified separately and then combined, reconciling persistent named-H2* NH entry with strongly variable product performance.

**Figure 13. Direct-NH dominance margin.** The ratio of integrated named-H2* direct-NH flux to the largest secondary direct-NH source remains above unity throughout the accepted map, supplying a quantitative no-crossover criterion.

**Figure 14. Flux-resolved cascade from NH entry to NH3 turnover.** The diagram follows terminal flux through NH, NH2, and NH3 formation and loss, showing the downstream location at which direct-NH dominance becomes decoupled from product accumulation.

**Figure 15. Product-turnover control phase space.** Formation-favoured and loss-dominated regions are distinguished using gross turnover, destruction fraction, and local kinetic controls. The high-field regime is characterized by strong ion-mediated NH3 loss.

**Figure 16. Numerical and mechanistic evidence hierarchy.** P0 convergence, direct-NH hierarchy, causal disabling, endpoint checks, atom-balance validation, turnover analysis, and kinetic-control results are assembled into an evidence hierarchy for the central conditional conclusion. This synthesis panel organizes the supporting evidence rather than adding an independent data set.

**Figure 17. Targeted residence-time robustness of the two-layer mechanism.** Four contrasting continuous-wave CSTR states are evaluated from 0.5 to 30 ms, and all 20 cases satisfy P0. Residence time changes NH3 productivity and the formation–destruction balance, particularly in the high-field state, while the named electronic-H2* direct-NH family remains above the largest secondary source throughout the panel. This is a targeted transport screen, not a full reactor response surface.

**Figure 18. Deterministic H2*-reaction scenario envelope.** The named-H2* to NH channels, named-H2* wall relaxation, and Rydberg-H2* to NH are independently scaled by \(1/3\) and 3 at the productivity maximum and the high-field loss-dominated state. All 14 calculations satisfy P0, and the named electronic-H2* family remains above the largest secondary direct-NH source. The dashed line marks the crossover criterion. This is an explicit deterministic scenario envelope, not statistical uncertainty quantification.

**Figure 19. Graph-informed model-discrimination atlas.** Three contrasting accepted continuous-wave CSTR states are selected using direct-NH hierarchy, dominance margin, and downstream turnover contrast. The figure pre-specifies observables and a falsification logic; it is a verification design derived from accepted model evidence, not an experimental result.

**Figure 20. Mechanism-family and ensemble-readiness audit.** The corrected Hong and Fe(110) DFT–microkinetic inputs are compared at reaction-family level. Structural correspondence is distinguished from kinetic and reactor-boundary equivalence. The Fe(110) branch is source-validated but is not a cross-boundary flux comparator.

**Figure 21. Fe(110) source scope, author-boundary reproduction, and porting gate.** The public Shao–Mesbah input is reproduced at its author boundary and its gas/backbone versus explicit-surface roles are separated. This audit establishes source traceability, not an independent validation of the Hong gas-phase NH-entry rate or a surface-flux result.

**Figure 22. Boundary-harmonization and pathway-identifiability protocol for gas-phase Hong and Fe(110) DFT–microkinetic models.** Gas \(E/N\) and local surface \(F_s\) occupy different physical scales and must be linked by \(F_s=\beta E_{\mathrm{bulk}}\). The five gates, shared-boundary ledger, and pre-registered flux/observable metrics prevent an invalid gas–surface comparison. This is a comparison guardrail, not a surface-flux result.

**Figure 23. Graph-theoretic audit of the explicit Fe(110) surface branch.** The source contains 46 explicit surface records and five structural first-surface-N–H entries: four Eley–Rideal and one Langmuir–Hinshelwood route. Record counts and connectivity describe encoded alternatives only; no reaction rate or route fraction is implied.

**Figure 24. Automorphism-orbit analysis of the source-locked Fe(110) surface reaction topology.** The unweighted directed bipartite graph contains verified \(N_2\), electronic-\(H_2\), vibrational-\(H_2\), and active-N topology orbits. State-specific energetics, populations, local fields and transport break these symmetries before flux ranking, so no kinetic degeneracy factor is introduced.

**Figure 25. Prospective blocked factorial validation design for gas–surface pathway discrimination.** A \(2^3\) design crosses Fe(110)-representative versus inert surface, independently calibrated low/high local field, and \(^{14}\mathrm{N_2}\)/\(^{15}\mathrm{N_2}\) feed within three randomized independent blocks. The design defines independent units and measurement gates; it is not an experimental dataset or a claimed surface-flux result.

**Figure 26. Direct-NH pathway concentration and runner-up identity across the accepted continuous-wave CSTR map.** Shannon entropy and its effective-family transform are calculated from the same six terminal integrated direct-NH family shares used in Figures 3 and 13. The runner-up changes among N(2P) + H2, Rydberg-H2*, and association, but the named electronic-H2* family remains at least 10.054666 times the runner-up. This is a descriptive re-expression of existing ROP evidence, not an independent validation, kinetic uncertainty analysis, or a surface-pathway result.

## Data and code availability

The source data, analysis scripts, figure-generation workflows, frozen external-input provenance ledgers, reconstructed mechanism inputs, and local reproducibility manifests underlying Figures 1–26 and the graphical abstract are retained with this project. A citable public archive, versioned release, and permanent identifier must be created before external submission; the final manuscript will replace this statement with the archive DOI and release version. The results reported here exclude any calculation that did not satisfy the declared P0 acceptance criterion.

## Declarations

### Ethics statement

This computational study did not involve human participants, animals, clinical data, or identifiable personal information; ethics approval was not required.

### Author contributions

Author order and CRediT roles are to be completed after confirmation by the authors. The final submission must identify contributions to conceptualization, methodology, software, validation, formal analysis, visualization, writing, supervision, and funding acquisition where applicable.

### Funding

TBD. The final submission must state either the applicable funding agencies and grant numbers or that the work received no external funding.

### Conflict of interest

TBD. The authors must declare all financial and non-financial competing interests before submission.

### AI-assisted writing disclosure

Generative AI was used as a drafting and language-editing aid under author supervision. The authors remain responsible for the scientific content, numerical results, source verification, and final wording. This disclosure should be adapted to the selected journal’s current policy before submission.

## References

[1] Rouwenhorst K H R, Engelmann Y, van ’t Veer K, Postma R S, Bogaerts A and Lefferts L 2020 Plasma-driven catalysis: green ammonia synthesis with intermittent electricity *Green Chemistry* **22** 6258–87 https://doi.org/10.1039/D0GC02058C

[2] Winter L R and Chen J G 2021 N2 fixation by plasma-activated processes *Joule* **5** 300–15 https://doi.org/10.1016/j.joule.2020.11.009

[3] Qu Z, Zhou R, Sun J, Gao Y, Li Z, Zhang T, Zhou R, Liu D, Tu X, Cullen P and Ostrikov K 2023 Plasma-assisted sustainable nitrogen-to-ammonia fixation: mixed-phase, synergistic processes and mechanisms *ChemSusChem* **17** e202300783 https://doi.org/10.1002/cssc.202300783

[4] Mehta P, Barboun P, Herrera F A, Kim J, Rumbach P, Go D B, Hicks J C and Schneider W F 2018 Overcoming ammonia synthesis scaling relations with plasma-enabled catalysis *Nature Catalysis* **1** 269–75 https://doi.org/10.1038/s41929-018-0045-1

[5] Rouwenhorst K H R, Kim H-H and Lefferts L 2019 Vibrationally excited activation of N2 in plasma-enhanced catalytic ammonia synthesis: a kinetic analysis *ACS Sustainable Chemistry & Engineering* **7** 17515–22 https://doi.org/10.1021/acssuschemeng.9b04997

[6] Hong J, Pancheshnyi S, Tam E, Lowke J J, Prawer S and Murphy A B 2017 Kinetic modelling of NH3 production in N2–H2 non-equilibrium atmospheric-pressure plasma catalysis *Journal of Physics D: Applied Physics* **50** 154005 https://doi.org/10.1088/1361-6463/aa6229

[7] Hong J, Pancheshnyi S, Tam E, Lowke J J, Prawer S and Murphy A B 2018 Corrigendum: Kinetic modelling of NH3 production in N2–H2 non-equilibrium atmospheric-pressure plasma catalysis (2017 *J. Phys. D: Appl. Phys.* **50** 154005) *Journal of Physics D: Applied Physics* **51** 109501 https://doi.org/10.1088/1361-6463/aaa988

[8] van ’t Veer K, Reniers F and Bogaerts A 2020 Zero-dimensional modeling of unpacked and packed bed dielectric barrier discharges: the role of vibrational kinetics in ammonia synthesis *Plasma Sources Science and Technology* **29** 045020 https://doi.org/10.1088/1361-6595/ab7a8a

[9] Sun J, Chen Q, Zhao X, Lin H and Qin W 2022 Kinetic investigation of plasma catalytic synthesis of ammonia: insights into the role of excited states and plasma-enhanced surface chemistry *Plasma Sources Science and Technology* **31** 094009 https://doi.org/10.1088/1361-6595/ac8e2c

[10] Chen Z, Koel B E and Sundaresan S 2022 Plasma-assisted catalysis for ammonia synthesis in a dielectric barrier discharge reactor: key surface reaction steps and potential causes of low energy yield *Journal of Physics D: Applied Physics* **55** 055202 https://doi.org/10.1088/1361-6463/ac2f12

[11] Wang Y, Craven M, Yu X, Ding J, Bryant P, Huang J and Tu X 2019 Plasma-enhanced catalytic synthesis of ammonia over a Ni/Al2O3 catalyst at near-room temperature: insights into the importance of the catalyst surface on the reaction mechanism *ACS Catalysis* **9** 10780–93 https://doi.org/10.1021/acscatal.9b02538

[12] Vervloedt S C L and von Keudell A 2024 Ammonia synthesis by plasma catalysis in an atmospheric RF helium plasma *Plasma Sources Science and Technology* **33** 045005 https://doi.org/10.1088/1361-6595/ad38d6

[13] Sun J, Wang W, Lu C and Tu X 2025 Elevated pressure effects on plasma-driven ammonia synthesis: insights from experiments and kinetic modeling *ACS Sustainable Chemistry & Engineering* **13** 15576–87 https://doi.org/10.1021/acssuschemeng.5c06251

[14] Vodlan K, Likozar B and Huš M 2025 Modeling of plasma-activated ammonia synthesis *Chemical Engineering Journal* **509** 161459 https://doi.org/10.1016/j.cej.2025.161459

[15] Holmes T D, Rothman R H and Zimmerman W B 2021 Graph theory applied to plasma chemical reaction engineering *Plasma Chemistry and Plasma Processing* **41** 531–57 https://doi.org/10.1007/s11090-021-10152-z

[16] Shao K and Mesbah A 2024 A study on the role of electric field in low-temperature plasma catalytic ammonia synthesis via integrated density functional theory and microkinetic modeling *JACS Au* **4** 525–44 https://doi.org/10.1021/jacsau.3c00654

[17] Gong F, Jing Y and Xiao R 2024 Plasma-assisted ammonia synthesis under mild conditions for hydrogen and electricity storage: mechanisms, pathways, and application prospects *Frontiers in Energy* **18** 418–35 https://doi.org/10.1007/s11708-024-0949-1

[18] Zheng Y, Meng S, Hao Y, Yi Y, Yang D and Cui Z 2025 Dielectric barrier discharge plasma-catalytic ammonia synthesis: from tools to mechanism *Plasma Processes and Polymers* **22** e70075 https://doi.org/10.1002/ppap.70075

[19] Gorbanev Y, Fedirchyk I and Bogaerts A 2024 Plasma catalysis in ammonia production and decomposition: use it, or lose it? *Current Opinion in Green and Sustainable Chemistry* **47** 100916 https://doi.org/10.1016/j.cogsc.2024.100916

[20] Mehta P, Barboun P M, Engelmann Y, Go D B, Bogaerts A, Schneider W F and Hicks J C 2020 Plasma-catalytic ammonia synthesis beyond the equilibrium limit *ACS Catalysis* **10** 6726–34 https://doi.org/10.1021/acscatal.0c00684

[21] Engelmann Y, van ’t Veer K, Gorbanev Y, Neyts E C, Schneider W F and Bogaerts A 2021 Plasma catalysis for ammonia synthesis: a microkinetic modeling study on the contributions of Eley–Rideal reactions *ACS Sustainable Chemistry & Engineering* **9** 13151–63 https://doi.org/10.1021/acssuschemeng.1c02713

[22] van ’t Veer K, Engelmann Y, Reniers F and Bogaerts A 2020 Plasma-catalytic ammonia synthesis in a DBD plasma: role of microdischarges and their afterglows *The Journal of Physical Chemistry C* **124** 22871–83 https://doi.org/10.1021/acs.jpcc.0c05110
