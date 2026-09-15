# Hong 2017/2018 corrected kinetic-input reconstruction

This directory is a fresh, literature-traceable reconstruction of the kinetic model in:

- J. Hong *et al.*, *J. Phys. D: Appl. Phys.* **50** (2017) 154005, DOI 10.1088/1361-6463/aa6229.
- J. Hong *et al.*, corrigendum, *J. Phys. D: Appl. Phys.* **51** (2018) 109501, DOI 10.1088/1361-6463/aaa988.

`kinet_hong_2017_2018_corrected.inp` is the main ZDPlasKin input. It contains all species listed in corrected Table 1; electron-impact reactions through BOLSIG+; the W11-W16 vibration model; W17-W140 gas-phase/wall reactions; and S1-S18 surface reactions for the **metal** column of Table 5.

`bolsigdb.dat` is included because the paper delegates W1-W3, W9-W10, W22 and W130 to BOLSIG+ rather than printing numerical rate coefficients. It is a runtime dependency, not a claim that every cross section in it is uniquely specified by the two papers.

The result is a reconstruction, not an author-supplied archive: neither paper publishes its exact `kinet.inp`, BOLSIG+ database revision, initial surface coverage, or the numerical convention that turns the incident normal energy `Ez` in equation (38) into a gas-condition value. See `sources_and_assumptions.md` before interpreting numerical results.
