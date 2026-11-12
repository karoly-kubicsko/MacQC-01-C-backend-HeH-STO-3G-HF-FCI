# Basis sets

The Gaussian basis sets of ΨLab, in Gaussian94 format (`.gbs`), exactly as written by the
[Basis Set Exchange](https://www.basissetexchange.org) (BSE, data version 0.12), downloaded on
2026-11-12 with

```
https://www.basissetexchange.org/api/basis/<name>/format/gaussian94/?elements=<Z,...>
```

ΨLab reads these files at run time (`backend/src/basis.cpp`); they are not edited by hand. The
functions are used as spherical harmonics (5d, 7f, 9g, 11h), as the correlation-consistent sets are
designed. Contracted functions are normalised; the contraction coefficients refer to normalised
primitives (Gaussian94 convention).

| File | Basis set | H | Li | O | highest l |
|---|---|---|---|---|---|
| `aug-cc-pvdz.gbs`  | aug-cc-pVDZ  | [3s2p] | [4s3p2d] | [4s3p2d] | d |
| `aug-cc-pvtz.gbs`  | aug-cc-pVTZ  | [4s3p2d] | [5s4p3d2f] | [5s4p3d2f] | f |
| `aug-cc-pvqz.gbs`  | aug-cc-pVQZ  | [5s4p3d2f] | [6s5p4d3f2g] | [6s5p4d3f2g] | g |
| `aug-cc-pv5z.gbs`  | aug-cc-pV5Z  | [6s5p4d3f2g] | [7s6p5d4f3g2h] | [7s6p5d4f3g2h] | h |
| `aug-cc-pcvdz.gbs` | aug-cc-pCVDZ | – | [5s4p2d] | [5s4p2d] | d |
| `aug-cc-pcvtz.gbs` | aug-cc-pCVTZ | – | [7s6p4d2f] | [7s6p4d2f] | f |
| `aug-cc-pcvqz.gbs` | aug-cc-pCVQZ | – | [9s8p6d4f2g] | [9s8p6d4f2g] | g |
| `aug-cc-pcv5z.gbs` | aug-cc-pCV5Z | – | – | [11s10p8d6f4g2h] | h |

**Element coverage.**

- The core-valence sets aug-cc-pCVXZ are not defined for hydrogen, which has no core electrons.
  As is customary, ΨLab uses the matching aug-cc-pVXZ functions for H (e.g. aug-cc-pVDZ for H in
  an aug-cc-pCVDZ calculation) and says so in the output.
- aug-cc-pCV5Z is defined in the BSE for B–Ar only, not for Li. ΨLab therefore refuses LiH in
  aug-cc-pCV5Z with an explanation (use aug-cc-pCVQZ or aug-cc-pV5Z).

Number of basis functions per molecule (spherical):

| | aug-cc-pVDZ | TZ | QZ | 5Z | aug-cc-pCVDZ | CVTZ | CVQZ | CV5Z |
|---|---|---|---|---|---|---|---|---|
| H2  | 18 | 46 | 92 | 160 | 18 | 46 | 92 | 160 |
| LiH | 32 | 69 | 126 | 207 | 36 | 82 | 155 | – |
| H2O | 41 | 92 | 172 | 287 | 45 | 105 | 201 | 341 |

## References

Please cite the BSE and the original publications when you use these basis sets:

- B. P. Pritchard, D. Altarawy, B. Didier, T. D. Gibson, T. L. Windus, *A New Basis Set
  Exchange: An Open, Up-to-date Resource for the Molecular Sciences Community*,
  J. Chem. Inf. Model. **59**, 4814 (2019), doi:10.1021/acs.jcim.9b00725.
- D. Feller, J. Comput. Chem. **17**, 1571 (1996), doi:10.1002/(SICI)1096-987X(199610)17:13<1571::AID-JCC9>3.0.CO;2-P;
  K. L. Schuchardt et al., J. Chem. Inf. Model. **47**, 1045 (2007), doi:10.1021/ci600510j.
- H, O (cc-pVXZ): T. H. Dunning, Jr., J. Chem. Phys. **90**, 1007 (1989), doi:10.1063/1.456153.
- H, O (diffuse "aug" functions): R. A. Kendall, T. H. Dunning, Jr., R. J. Harrison,
  J. Chem. Phys. **96**, 6796 (1992), doi:10.1063/1.462569.
- Li (aug-cc-pVXZ): B. P. Prascher, D. E. Woon, K. A. Peterson, T. H. Dunning, Jr.,
  A. K. Wilson, Theor. Chem. Acc. **128**, 69 (2011), doi:10.1007/s00214-010-0764-0.
- O (core-valence functions, cc-pCVXZ): D. E. Woon, T. H. Dunning, Jr., J. Chem. Phys. **103**,
  4572 (1995), doi:10.1063/1.470645.
- Li (aug-cc-pCVXZ, X = D, T, Q): from the original BSE (Feller 1996; Schuchardt et al. 2007),
  diffuse functions from Prascher et al. (2011).
