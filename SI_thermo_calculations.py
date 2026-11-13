#!/usr/bin/env python3
"""
Supporting Information script for
"A Practical Guide to Thermodynamically Grounded Primer Design for Multiplex
and Ultra-Multiplex PCR" (Kalendar & Kairov).

Reproduces every worked calculation in the manuscript and in Tables S2-S6:
  * linguistic complexity (LC%) and purine-pyrimidine complexity (YR%)
  * nearest-neighbor Tm and dG37 under the design conditions used by FastPCR /
    PCRpanel, and their dependence on primer and free Mg2+ concentration
  * average versus worst-case stability of a degenerate primer
  * the relationship between dG37 and Tm thresholds for primer dimers

Python 3.8+, standard library only.  Run:  python3 SI_thermo_calculations.py

Model (identical to the software defaults):
  - Unified DNA/DNA nearest-neighbor parameters (SantaLucia 1998; Allawi &
    SantaLucia 1997), including initiation terms and the symmetry correction.
  - dH and dS refer to 1 M NaCl.  Tm(1 M) = dH / (dS + R ln(C/x)), where C is
    the concentration of each strand and x = 2 for two non-self-complementary
    strands at equal concentration (equivalently C_T/4 with C_T = 2C),
    x = 1 for a self-complementary duplex.
  - Salt: Owczarzy et al. 2008 Mg2+ correction (their eq. 16, with the
    monovalent-dependent coefficients a, d, g for 0.22 <= R < 6), applied for
    55 mM monovalent cation (50 mM K+ + ~5 mM Tris+) and 1.0 mM FREE Mg2+.
  - Design concentration C = 0.2 uM.
  - dG37 reported is the standard free energy at 37 C and 1 M NaCl, from the
    unified dG37 table (it contains no concentration term).
With these settings the Tm values of all 28 primers in Table S1 are
reproduced to within +/-0.05 C and the dG37 values to within +/-0.06 kcal/mol.
"""
import itertools
import math
import random
import statistics as st

R = 1.987  # cal mol-1 K-1

# ---------------------------------------------------------------- parameters
NN = {  # (dH kcal/mol, dS cal/mol/K), 5'->3' top-strand dinucleotide
    'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2), 'AT': (-7.2, -20.4),
    'TA': (-7.2, -21.3), 'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4), 'CT': (-7.8, -21.0),
    'AG': (-7.8, -21.0), 'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4), 'GG': (-8.0, -19.9),
    'CC': (-8.0, -19.9)}
G37 = {'AA': -1.00, 'TT': -1.00, 'AT': -0.88, 'TA': -0.58, 'CA': -1.45,
       'TG': -1.45, 'GT': -1.44, 'AC': -1.44, 'CT': -1.28, 'AG': -1.28,
       'GA': -1.30, 'TC': -1.30, 'CG': -2.17, 'GC': -2.24, 'GG': -1.84,
       'CC': -1.84}
INIT = {'GC': (0.1, -2.8, 0.98), 'AT': (2.3, 4.1, 1.03)}  # dH, dS, dG37
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT',
         'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
         'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}

DESIGN = dict(C=2e-7, mon=0.055, mg=0.001)  # FastPCR / PCRpanel defaults


def rc(s):
    return s.upper()[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def self_compl(s):
    return s.upper() == rc(s)


# ------------------------------------------------------- nearest-neighbor
def dH_dS(seq):
    """dH (kcal/mol) and dS (cal/mol/K) at 1 M NaCl; IUPAC codes averaged."""
    s = seq.upper()
    H = S = 0.0
    for i in range(len(s) - 1):
        combos = [(a, b) for a in IUPAC[s[i]] for b in IUPAC[s[i + 1]]]
        H += sum(NN[a + b][0] for a, b in combos) / len(combos)
        S += sum(NN[a + b][1] for a, b in combos) / len(combos)
    for end in (s[0], s[-1]):
        bases = IUPAC[end]
        for b in bases:
            h, e, _ = INIT['GC' if b in 'GC' else 'AT']
            H += h / len(bases)
            S += e / len(bases)
    if set(s) <= set('ACGT') and self_compl(s):
        S += -1.4
    return H, S


def dG37(seq):
    """Standard free energy at 37 C, 1 M NaCl (kcal/mol), unified dG37 table."""
    s = seq.upper()
    g = sum(G37[s[i:i + 2]] for i in range(len(s) - 1))
    for end in (s[0], s[-1]):
        g += INIT['GC' if end in 'GC' else 'AT'][2]
    if self_compl(s):
        g += 0.43
    return g


def f_gc(seq):
    s = seq.upper()
    return sum(sum(b in 'GC' for b in IUPAC[c]) / len(IUPAC[c]) for c in s) / len(s)


def owczarzy2008(T1, n_bp, fgc, mon, mg):
    """Owczarzy et al. 2008 Mg2+ correction of a 1 M NaCl Tm (K -> K)."""
    a, b, c, d = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5
    e, f, g = -4.82e-4, 5.25e-4, 8.31e-5
    ratio = math.sqrt(mg) / mon
    if ratio < 0.22:  # monovalent ions dominate: Owczarzy 2004
        ln = math.log(mon)
        return 1 / (1 / T1 + (4.29 * fgc - 3.95) * 1e-5 * ln + 9.40e-6 * ln * ln)
    if ratio < 6.0:
        lm = math.log(mon)
        a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mon) * lm)
        d = 1.42e-5 * (1.279 - 4.03e-3 * lm - 8.03e-3 * lm * lm)
        g = 8.31e-5 * (0.486 - 0.258 * lm + 5.25e-3 * lm ** 3)
    lmg = math.log(mg)
    inv = (1 / T1 + a + b * lmg + fgc * (c + d * lmg)
           + (e + f * lmg + g * lmg * lmg) / (2 * (n_bp - 1)))
    return 1 / inv


def tm(seq, C=DESIGN['C'], mon=DESIGN['mon'], mg=DESIGN['mg']):
    """Tm (C) of a duplex whose strands are each at concentration C (M)."""
    H, S = dH_dS(seq)
    x = 1 if (set(seq.upper()) <= set('ACGT') and self_compl(seq)) else 2
    T1 = 1000 * H / (S + R * math.log(C / x))
    return owczarzy2008(T1, len(seq), f_gc(seq), mon, mg) - 273.15


def occupancy(seq, T_c, C, ddG37=0.0, mon=DESIGN['mon'], mg=DESIGN['mg']):
    """Fraction of template sites bound by primer (primer in excess, conc C)
    at T_c (C).  Salt-corrected dS is the value that reproduces tm(seq, C).
    ddG37 > 0 destabilises the duplex (e.g. a mismatch), taken as entropic."""
    H, _ = dH_dS(seq)
    Tm_K = tm(seq, C, mon, mg) + 273.15
    S_eff = 1000 * H / Tm_K - R * math.log(C / 2)      # cal/mol/K
    S_eff -= 1000 * ddG37 / 310.15
    T = T_c + 273.15
    K = math.exp(-(1000 * H - T * S_eff) / (R * T))
    return K * C / (1 + K * C)


# ------------------------------------------------------- complexity
def complexity(seq, binary=False, K=None):
    """LC% (4-letter, K=3) or YR% (binary R/Y, K=6) as implemented.
    LC = 100 * sum_{k=1..K} U_k / sum_{k=1..K} min(a^k, N-k+1),
    U_k = number of distinct words of length k observed, N = length,
    a = 4 (LC) or 2 (YR).  K = 3 (LC) and 6 (YR) for oligonucleotides up to
    60 nt; 4 and 8 for longer sequences."""
    s = seq.upper()
    a = 4
    if binary:
        s = ''.join('R' if c in 'AG' else 'Y' for c in s)
        a = 2
    N = len(s)
    if K is None:
        K = (3 if N <= 60 else 4) * (2 if binary else 1)
    rows, obs, exp = [], 0, 0
    for k in range(1, K + 1):
        if N - k + 1 <= 0:
            break
        u = len({s[i:i + k] for i in range(N - k + 1)})
        m = min(a ** k, N - k + 1)
        rows.append((k, u, m))
        obs += u
        exp += m
    return 100 * obs / exp, obs, exp, rows


# ------------------------------------------------------- data
TABLE_S1 = [  # id, sequence, Tm and dG37 as listed in Table S1
    ('5820', 'cgaagtcgcctcccacttttggt', 68.0, -30.9, 93, 78),
    ('5853', 'AGACCGGCTGTCCAcgaagtcgcctcccacttttggt', 80.1, -53.1, 87, 76),
    ('5854', 'tccccgaagtcgcctcccacttttggt', 74.2, -37.7, 87, 73),
    ('5800', 'ccgaagtcgcctcccacttt', 65.7, -27.0, 89, 80),
    ('5832', 'cctaagtccccgaagtcgcctcccacttt', 73.6, -39.4, 79, 73),
    ('5828', 'AGACCGGCTGTCCAccgaagtcgcctcccacttt', 79.7, -49.2, 85, 72),
    ('5802', 'cgcagcttctctgtggagacc', 66.1, -28.2, 90, 86),
    ('5834', 'gtgctcttcgcagcttctctgtggagacc', 72.7, -39.5, 79, 76),
    ('5830', 'AGACCGGCTGTCCAcgcagcttctctgtggagacc', 79.4, -50.4, 81, 88),
    ('5848', 'cccgaagtcgcctcccactt', 67.9, -27.8, 87, 84),
    ('5849', 'cgccagtataggccgaccct', 67.2, -27.5, 87, 90),
    ('5850', 'cttctctgtggagaccggtcagc', 67.6, -30.5, 90, 82),
    ('5857', 'gcagcttctctgtggagaccggtcagc', 73.3, -37.7, 87, 84),
    ('5851', 'ctggtggcatccgcacagc', 67.8, -27.3, 84, 85),
    ('5852', 'ccgggtcagaattggcggct', 68.8, -28.4, 89, 75),
    ('5801', 'actgtaggtcacaggctctacgt', 65.4, -29.2, 85, 75),
    ('5833', 'ACTTTTTCCTTTTCactgtaggtcacaggctctacgt', 71.3, -46.1, 76, 64),
    ('5829', 'AGACCGGCTGTCCAactgtaggtcacaggctctacgt', 77.7, -51.0, 85, 73),
    ('5803', 'cggttagctggtggattgctaca', 65.4, -29.5, 85, 78),
    ('5835', 'ctatttcaggcggttagctggtggattgctaca', 71.3, -42.3, 86, 84),
    ('5831', 'AGACCGGCTGTCCAcggttagctggtggattgctaca', 78.4, -51.7, 87, 71),
    ('5821', 'tcagagccacatgtatcttccggttt', 66.4, -32.4, 98, 89),
    ('5822', 'gcccacaatggttaagcagaatgcaa', 66.8, -33.0, 84, 92),
    ('5823', 'gggagaagacacaggtggcaataaca', 67.2, -32.8, 82, 58),
    ('5824', 'cgttgcccatgctggtatcaaacc', 67.2, -31.5, 93, 87),
    ('5825', 'tggcaaagaagtatccacactcgtcc', 67.1, -33.0, 91, 79),
    ('5826', 'tgcttctgtggtgaaggagctttcat', 67.0, -32.7, 82, 82),
    ('5827', 'cctggaaaggccactttgtaagctca', 67.6, -33.1, 93, 90),
]
CONCS_nM = [1000, 500, 200, 100, 50, 30, 15, 10, 5]


def section(title):
    print('\n' + '=' * 78 + '\n' + title + '\n' + '=' * 78)


def main():
    # ---------------------------------------------------- S2: Tm vs [primer]
    section('Table S2. Tm (C) versus primer concentration; 55 mM monovalent, '
            '1.0 mM free Mg2+')
    dev_tm, dev_dg, dev_lc, dev_yr = [], [], [], []
    print('ID\tnt\tdH\tdS\tdG37\t' + '\t'.join(f'{c} nM' for c in CONCS_nM)
          + '\td(500->50)\td(500->10)\td(200->10)')
    for pid, seq, tm_s1, dg_s1, lc_s1, yr_s1 in TABLE_S1:
        H, S = dH_dS(seq)
        t = {c: tm(seq, c * 1e-9) for c in CONCS_nM}
        dev_tm.append(t[200] - tm_s1)
        dev_dg.append(dG37(seq) - dg_s1)
        dev_lc.append(round(complexity(seq)[0]) - lc_s1)
        dev_yr.append(round(complexity(seq, True)[0]) - yr_s1)
        print(f'{pid}\t{len(seq)}\t{H:.1f}\t{S:.1f}\t{dG37(seq):.1f}\t'
              + '\t'.join(f'{t[c]:.1f}' for c in CONCS_nM)
              + f'\t{t[500]-t[50]:.2f}\t{t[500]-t[10]:.2f}\t{t[200]-t[10]:.2f}')
    print(f'\nReproduction of Table S1: Tm max |dev| = {max(map(abs, dev_tm)):.3f} C;'
          f' dG37 max |dev| = {max(map(abs, dev_dg)):.3f} kcal/mol;'
          f' LC max |dev| = {max(map(abs, dev_lc))}; YR max |dev| = {max(map(abs, dev_yr))}')
    for a, b in [(500, 50), (500, 10), (200, 10), (200, 15), (50, 10), (100, 30)]:
        d = [tm(s, a * 1e-9) - tm(s, b * 1e-9) for _, s, *_ in TABLE_S1]
        print(f'Tm shift {a} -> {b} nM: mean {st.mean(d):.2f} C, '
              f'range {min(d):.2f}-{max(d):.2f} C (n = {len(d)})')
    H, _ = dH_dS('ccgaagtcgcctcccacttt')
    T = tm('ccgaagtcgcctcccacttt', 5e-7) + 273.15
    print(f'Analytical check, 5800: dTm/dlog10C ~ ln(10) R Tm^2/|dH| = '
          f'{math.log(10) * R * T * T / (-1000 * H):.2f} C per decade')

    # ---------------------------------------------------- S3: free Mg2+
    section('Table S3. Effect of the free Mg2+ concentration (C = 200 nM, 55 mM '
            'monovalent); dTm relative to 1.0 mM free Mg2+')
    total_mg, dntp = 1.5, 0.8
    print(f'Total Mg2+ {total_mg} mM, total dNTP {dntp} mM (0.2 mM each), '
          f'1:1 binding -> free Mg2+ ~ {total_mg - dntp:.1f} mM')
    for mg in [0.5, 0.7, 1.0, 1.5, 2.0]:
        d = [tm(s, 2e-7, 0.055, mg * 1e-3) - tm(s) for _, s, *_ in TABLE_S1]
        print(f'free Mg2+ {mg:.1f} mM: mean {st.mean(d):+.2f} C, '
              f'range {min(d):+.2f} to {max(d):+.2f} C')
    d = [tm(s, 2e-7, 0.050, 1e-3) - tm(s) for _, s, *_ in TABLE_S1]
    print(f'monovalent 50 vs 55 mM: mean {st.mean(d):+.2f} C')
    # design Tm 60-62 C -> working conditions of the multiplex experiments
    random.seed(7)
    shifts = []
    while len(shifts) < 400:
        L = random.randint(18, 26)
        s = ''.join(random.choice('ACGT') for _ in range(L))
        if 60 <= tm(s) <= 62:
            for C in (10e-9, 15e-9):
                shifts.append(tm(s) - tm(s, C, 0.055, 0.0007))
    print(f'Primers with design Tm 60-62 C (18-26 nt, n = 200): working Tm at '
          f'10-15 nM and 0.7 mM free Mg2+ is lower by {st.mean(shifts):.1f} C '
          f'(range {min(shifts):.1f}-{max(shifts):.1f} C)')

    # ---------------------------------------------------- S4: LC / YR
    section('Table S4. Worked examples of LC% (K = 3) and YR% (K = 6)')
    examples = [
        ('poly(A), 18 nt', 'A' * 18), ('poly(A), 20 nt', 'A' * 20),
        ('(CA)10, alternating Y/R', 'CA' * 10),
        ('(AG)8C, ISSR-type', 'AG' * 8 + 'C'),
        ('(ACGT)5', 'ACGT' * 5),
        ('polypurine tract', 'GAGGAAGAAGAGAAGGAGGA'),
        ('primer 5821', 'tcagagccacatgtatcttccggttt'),
        ('primer 5822', 'gcccacaatggttaagcagaatgcaa'),
        ('primer 5823', 'gggagaagacacaggtggcaataaca'),
        ('primer 5833 (tailed)', 'ACTTTTTCCTTTTCactgtaggtcacaggctctacgt'),
    ]
    for name, s in examples:
        lc, o, e, rows = complexity(s)
        yr, oy, ey, rowsy = complexity(s, True)
        print(f'{name:26s} N={len(s):2d}  LC {lc:5.1f}% = {o}/{e} '
              f'[{"; ".join(f"k{k}:{u}/{m}" for k, u, m in rows)}]   '
              f'YR {yr:5.1f}% = {oy}/{ey} '
              f'[{"; ".join(f"k{k}:{u}/{m}" for k, u, m in rowsy)}]')
    for N in (18, 20, 25, 30):
        print(f'Homopolymer lower bound, N = {N}: LC = {complexity("A" * N)[0]:.1f}%,'
              f' YR = {complexity("A" * N, True)[0]:.1f}%')
    random.seed(1)
    for N in (20, 25):
        v = [complexity(''.join(random.choice('ACGT') for _ in range(N)))[0]
             for _ in range(100000)]
        print(f'Random {N}-mers (n = 100000): mean LC {st.mean(v):.1f}%; '
              f'<70%: {100 * sum(x < 70 for x in v) / len(v):.1f}%; '
              f'<75%: {100 * sum(x < 75 for x in v) / len(v):.1f}%; '
              f'>=90%: {100 * sum(x >= 90 for x in v) / len(v):.1f}%; '
              f'=100%: {100 * sum(x > 99.99 for x in v) / len(v):.2f}%')

    # ---------------------------------------------------- S5: degenerate
    section('Table S5. Average versus worst-case stability of a degenerate primer')
    for name, s in [('806R', 'GGACTACNVGGGTWTCTAAT'),
                    ('515F-Y', 'GTGYCAGCMGCCGCGGTAA')]:
        iso = [''.join(p) for p in itertools.product(*[IUPAC[c] for c in s])]
        n = len(iso)
        t_avg = tm(s)
        t200 = [tm(x) for x in iso]
        t_own = [tm(x, 2e-7 / n) for x in iso]
        print(f'{name} {s}: {n} variants')
        print(f'  averaged parameters, 200 nM total: Tm {t_avg:.1f} C')
        print(f'  variants at 200 nM: min {min(t200):.1f} ({iso[t200.index(min(t200))]}), '
              f'mean {st.mean(t200):.1f}, max {max(t200):.1f} ({iso[t200.index(max(t200))]})')
        print(f'  variants at own concentration ({200 / n:.1f} nM): '
              f'{min(t_own):.1f}-{max(t_own):.1f} C')

    # ---------------------------------------------------- S6: dimers
    section('Table S6. dG37 (1 M NaCl) and Tm of short duplexes (55 mM '
            'monovalent, 1.0 mM free Mg2+)')
    print(f'Identity: at 1 M NaCl, C = 200 nM (C/2 = 0.1 uM), a duplex with '
          f'dG37 = {R * 310.15 * math.log(1e-7) / 1000:.2f} kcal/mol has Tm = 37 C '
          f'for any sequence.')
    duplexes = ['GCCGC', 'GCCGCG', 'GCCGCGG', 'GCGCGC', 'ACGTCA', 'GACGTCA',
                'GACGTCAG', 'TGACGTCAG', 'ATGACGTCAG', 'AAATTTAA',
                'AAAAAAAAAAAAA', 'TTAGCTAAT']
    print('duplex\tbp\tGC%\tdG37\tTm(1M,200nM)\tTm 200nM\tTm 50nM\tTm 10nM')
    for s in duplexes:
        H, S = dH_dS(s)
        x = 1 if self_compl(s) else 2
        t1 = 1000 * H / (S + R * math.log(2e-7 / x)) - 273.15
        print(f'{s}{" (self)" if self_compl(s) else ""}\t{len(s)}\t'
              f'{100 * f_gc(s):.0f}\t{dG37(s):.2f}\t{t1:.1f}\t{tm(s, 2e-7):.1f}\t'
              f'{tm(s, 5e-8):.1f}\t{tm(s, 1e-8):.1f}')
    random.seed(3)
    for C in (2e-7, 5e-8, 1e-8):
        for target in (37, 20):
            vals = []
            for L in range(6, 15):
                for _ in range(400):
                    s = ''.join(random.choice('ACGT') for _ in range(L))
                    if not self_compl(s) and abs(tm(s, C) - target) < 0.5:
                        vals.append(dG37(s))
            print(f'C = {C * 1e9:g} nM: duplexes (6-14 bp) with Tm = {target} +/- 0.5 C '
                  f'have dG37 from {min(vals):.1f} to {max(vals):.1f} kcal/mol (n = {len(vals)})')
    # A13 example quoted in the response to reviewers
    print(f'A13 at 100 nM: dG37 {dG37("A" * 13):.2f} kcal/mol, Tm {tm("A" * 13, 1e-7):.1f} C')

    # ---------------------------------------------------- annealing above Tm
    section('Annealing above Tm: template occupancy, mismatch discrimination, '
            'association time')
    s = 'ccgaagtcgcctcccacttt'  # 5800
    for C in (2e-7, 1e-8):
        print(f'5800 at {C * 1e9:g} nM: working Tm {tm(s, C):.1f} C')
        for T in (55, 60, 62, 64, 66, 68, 72):
            th = occupancy(s, T, C)
            thm = occupancy(s, T, C, ddG37=3.0)
            print(f'  {T} C: matched site occupancy {th:.3f}; site with a '
                  f'3 kcal/mol mismatch {thm:.4f}; discrimination {th / thm:.1f}x')
    for kon in (1e5, 1e6):
        print(f'k_on = {kon:.0e} M-1 s-1: half-time of association '
              + ', '.join(f'{c} nM {math.log(2) / (kon * c * 1e-9):.0f} s'
                          for c in (200, 50, 15, 10, 5)))
    # Rychlik et al. 1990: Ta,opt = 0.3 Tm(primer) + 0.7 Tm(product) - 14.9
    for L in (400, 1500, 3000):
        tp = 81.5 + 16.6 * math.log10(0.05) + 0.41 * 50 - 675 / L
        print(f'Rychlik Ta,opt, 50% GC product of {L} bp, primer Tm 67 C: '
              f'product Tm {tp:.1f} C, Ta,opt {0.3 * 67 + 0.7 * tp - 14.9:.1f} C')


if __name__ == '__main__':
    main()
