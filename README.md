# multiplex-pcr-thermo

Thermodynamic calculations for multiplex PCR primer design.

Supporting script for

> Kalendar R., Kairov U. **A Practical Guide to Thermodynamically Grounded Primer Design for Multiplex and Ultra-Multiplex PCR.** *ACS Measurement Science Au* (2026).

`SI_thermo_calculations.py` reproduces every calculated value in Supporting Information Tables S1–S6 of the paper and the worked examples in its main text. It needs only the Python standard library.

## Quick start

```bash
git clone https://github.com/rkalendar/multiplex-pcr-thermo.git
cd multiplex-pcr-thermo
python3 SI_thermo_calculations.py > output.txt
```

Python 3.8 or later. No packages to install. The run takes a few seconds and prints about 130 lines; the expected result is in [`SI_thermo_calculations_output.txt`](SI_thermo_calculations_output.txt). Random sampling uses fixed seeds, so the output is identical on every run.

## What it calculates

| Output section | Content | In the paper |
|---|---|---|
| Table S2 | Tm of the 28 study primers at 5–1000 nM, with ΔH°, ΔS° and ΔG°37 | Primer concentration |
| Table S3 | Shift in Tm with free Mg²⁺ from 0.5 to 2.0 mM; working Tm against design Tm | Multiplex PCR methods |
| Table S4 | Linguistic complexity (LC) and purine–pyrimidine complexity (YR), broken down by word length | Linguistic complexity |
| Table S5 | Average against worst-case Tm of degenerate primers (806R, 515F-Y) | Dimer formation |
| Table S6 | ΔG°37 against Tm for short duplexes at three concentrations | Dimer thresholds |
| Annealing above Tm | Template occupancy, mismatch discrimination, association half-times | Annealing temperature |

The first section also checks the script against the published Table S1: Tm agrees to within ±0.05 °C, ΔG°37 to within ±0.06 kcal/mol, and LC and YR exactly.

## Model

**Melting temperature.** Unified nearest-neighbor parameters for DNA duplexes (SantaLucia 1998), with initiation terms and the symmetry correction:

```
Tm(1 M NaCl) = ΔH° / (ΔS° + R·ln(C/x))
```

`C` is the concentration of each strand. `x = 2` for two different strands at equal concentration (the same as C_T/4 with C_T = 2C) and `x = 1` for a self-complementary duplex. The result is then corrected for salt with the Mg²⁺ equation of Owczarzy et al. (2008).

**Design conditions** (the defaults of FastPCR and PCRpanel, and of every function here):

| Quantity | Value |
|---|---|
| Monovalent cations | 55 mM (50 mM K⁺ + about 5 mM Tris⁺) |
| Free Mg²⁺ | 1.0 mM |
| Strand concentration | 0.2 µM each |

The Mg²⁺ value is the **free** concentration, not the total. dNTPs bind Mg²⁺ about 1:1, so a reaction with 1.5 mM total Mg²⁺ and 0.8 mM total dNTPs has about 0.7 mM free Mg²⁺.

**Free energy.** ΔG°37 is the standard free energy at 37 °C and 1 M NaCl. It has no concentration term.

**Complexity.** For a sequence of length *N* over an alphabet of *a* letters, with *U_k* distinct words of length *k*:

```
LC(%) = 100 × Σ U_k / Σ min(a^k, N − k + 1),   k = 1 … K
```

| Score | Alphabet | K (≤ 60 nt) | K (> 60 nt) |
|---|---|---|---|
| LC | A, C, G, T (a = 4) | 3 | 4 |
| YR | purine/pyrimidine (a = 2) | 6 | 8 |

A homopolymer gives the minimum, which is not zero: 8.3% for 18 nt.

**Degenerate bases.** IUPAC codes are handled by averaging ΔH° and ΔS° over the allowed bases. This gives the average stability of the mixture, not that of its most stable member (Table S5).

## Use as a module

```python
import SI_thermo_calculations as S

p = "ccgaagtcgcctcccacttt"            # primer 5800

S.tm(p)                               # 65.7  °C at design conditions
S.tm(p, C=1e-8)                       # 61.4  °C at 10 nM
S.tm(p, C=1e-8, mg=0.0007)            # 60.6  °C at 10 nM, 0.7 mM free Mg2+
S.dG37(p)                             # -27.0 kcal/mol
S.dH_dS(p)                            # (-157.6, -420.8)  kcal/mol, cal/mol/K

S.complexity(p)[:3]                   # (89.5, 34, 38)  LC %, observed, maximum
S.complexity(p, binary=True)[:3]      # (80.3, 49, 61)  YR %

S.tm("GGACTACNVGGGTWTCTAAT")          # 56.6  °C, averaged over 24 variants
S.occupancy(p, T_c=62, C=1e-8)        # 0.56  fraction of matched sites bound
```

Concentrations are in mol/L (`C` per strand, `mon` monovalent, `mg` free Mg²⁺). Sequences may be upper or lower case.

## Scope

The script is a transparent reference for the calculations in the paper, not a primer design tool.

- Only fully matched Watson–Crick duplexes are scored. Mismatches, dangling ends, hairpins and G-quadruplexes are not modelled.
- Duplex formation is treated as two-state, with ΔH° and ΔS° independent of temperature.
- The salt correction covers monovalent cations and Mg²⁺ only; dNTP binding, DMSO and other additives are left to the user through the free Mg²⁺ value.

For primer and panel design use [FastPCR](https://primerdigital.com/tools/fastpcr.html) or [PCRpanel](https://github.com/rkalendar/PCRpanel).

## References

1. SantaLucia J. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proc. Natl. Acad. Sci. U.S.A.* 1998, 95, 1460–1465.
2. Allawi H. T., SantaLucia J. Thermodynamics and NMR of internal G·T mismatches in DNA. *Biochemistry* 1997, 36, 10581–10594.
3. Owczarzy R., Moreira B. G., You Y., Behlke M. A., Walder J. A. Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations. *Biochemistry* 2008, 47, 5336–5353.
4. von Ahsen N., Wittwer C. T., Schütz E. Oligonucleotide melting temperatures under PCR conditions. *Clin. Chem.* 2001, 47, 1956–1961.

## Contact

Ruslan Kalendar — ruslan.kalendar@helsinki.fi

## License

GNU General Public License v3.0, see [LICENSE](LICENSE).
