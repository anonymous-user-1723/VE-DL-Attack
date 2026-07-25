# Verification Experiments for "Key-Recovery-Friendly Differential-Linear Cryptanalysis for ARX Ciphers"

This repository contains the code of the two verification experiments reported in our rebuttal response (Reviewer B, Q3):

1. **Splitting multiple modular subtractions** (`splitting-modular-subtraction/`): the proportion of valid samples when recovering a borrow bit by simultaneously splitting two (resp. three) modular subtractions, matching Lemma 4 (resp. Lemma 5) of the paper.
2. **Two-round partial decryption on LEA** (`lea-partial-decryption/`): a reduced 7-round key-recovery attack (5-round differential-linear distinguisher + the identical 2-round partial decryption of Figure 6), verifying the partial decryption procedure of Section 6.1 in an attack setting.

Requirements: a C++17 compiler (`g++` or `clang++`); `pthread` for the second experiment.

---

## Experiment 1: Splitting multiple modular subtractions

Directory: `splitting-modular-subtraction/`

For a chain of two (resp. three) consecutive modular subtractions as in Figure 4 of the paper, we measure the proportion of valid samples (those for which the target borrow bit can be determined) when the subtractions are split simultaneously, and compare it with the theoretical value of Lemma 4 (resp. Lemma 5). A sample is counted as valid only when the borrow bit is actually recovered by the attack-style (constructive) CDB procedure; each constructively valid sample is additionally sanity checked by key-perturbation invariance.

### Build and run

```bash
cd splitting-modular-subtraction
g++ -O2 -std=c++17 -o verify_two_modsub verify_two_modsub.cpp
g++ -O2 -std=c++17 -o verify_three_modsub verify_three_modsub.cpp

./verify_two_modsub --T 10000 --seed 123
./verify_three_modsub --T 10000 --seed 123
```

Both programs run a small grid of splitting parameters with `N = 100000` random samples per cell.

### Expected results (the rows reported in the rebuttal)

- two subtractions, `t=1, m=4`: `p_hat = 0.8111` vs. theory `(1-(1+2^-1*4)*2^-4) = 0.8125`;
- three subtractions, `s=1, t=1, m=4`: `p_hat = 0.6553` vs. theory `1-(1+4*2^-1*(1+(4+2*1-1)*2^-2))*2^-4 = 0.65625`.

All grid cells pass with `abs_err` well below the threshold. Lemmas 2, 4 and 5 actually state lower bounds on the expected number of valid samples: extra valid samples are also found in the experiments, making our estimates conservative (this does not affect our attacks; the `p_obo` column gives the one-by-one splitting baseline for comparison).

---

## Experiment 2: Two-round partial decryption on LEA

Directory: `lea-partial-decryption/`

We mount a reduced 7-round attack on LEA:

- **5-round differential-linear distinguisher**: input difference `Delta_in = (0, [30], 0, 0)`, output mask `Gamma_out = ([5,6,9,18], [3,4,24,27], [6,29], [0,29])` (identical to the paper), correlation `cor = 2^-2.72`;
- **2-round partial decryption** `x^7 -> x^5`: identical to Figure 6 of the paper (the same 85-bit key-guessing space, the same 7 borrow bits to recover, and the same valid-sample rule), implemented by known/unknown three-value propagation in `tristate.cpp` / `partial_dec.cpp`.

Key guessing: only **6 of the 85 key bits** are guessed (64 candidates); the other 79 bits are fixed to their correct values. The 6 guessed bit positions (canonical enumeration of the 85 bits) are `18, 57, 8, 64, 47, 23`.

Attack parameters from the paper's formula (1) with filtering strength `a = 30` and success rate `P_S = 0.9`:

| quantity | value |
|---|---|
| `N_v = ((Phi^-1(0.9)+Phi^-1(1-2^-30))/cor)^2` | 2308 |
| `N = N_v * 2^14.19` (pairs per attack) | 43130000 |
| statistic | `Q(k) = T(k)/N_v(k)` |
| threshold `t = Phi^-1(1-2^-30)/sqrt(N_v)` | 0.1251 |

For every valid pair of the correct candidate, the statistic is additionally compared against a full two-round decryption under the real round keys (embedded sanity check).

### Build and run (100 attacks)

```bash
cd lea-partial-decryption
g++ -O2 -std=c++17 -o attack attack.cpp partial_dec.cpp tristate.cpp lea.cpp -pthread
./attack --bits 18,57,8,64,47,23 --trials 100 --seed 100 --threads 8
```

Each attack uses a fresh random master key (the pair streams are regenerated deterministically from the seed, so no ciphertexts are stored; memory usage is below 20 MB). Running time: about 15 minutes on 8 threads; `--threads 1` also works but is proportionally slower.

### Expected results (as reported in the rebuttal)

- **90 of 100 attacks succeeded** (`summary: 90/100 trials succeeded`), i.e., the correct key guess passes the threshold `Q(k*) > t`;
- **no wrong key guess survived** (`wrong_pass = 0` in all 100 attacks, 6300 wrong candidates in total);
- the correct candidate ranks first in all 100 attacks;
- the embedded sanity check reports `sanity_mismatch = 0` everywhere, i.e., the statistic of every valid pair computed by the partial decryption exactly matches the full two-round decryption.
