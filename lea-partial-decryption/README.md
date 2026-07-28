# Two-round partial decryption on LEA

We mount a reduced 7-round attack on LEA:

- **5-round DL distinguisher**: input difference `Delta_in = (0, [30], 0, 0)`, output mask `Gamma_out = ([5,6,9,18], [3,4,24,27], [6,29], [0,29])` (identical to the paper), correlation `cor = 2^-2.72`;
- **2-round partial decryption** `x^7 -> x^5`: identical to Figure 6 of the paper (the same 85-bit key-guessing space, the same 7 borrow bits to recover, and the same valid-sample rule), implemented by known/unknown three-value propagation in `tristate.cpp` / `partial_dec.cpp`.

Key guessing: only **6 of the 85 key bits** are guessed (64 candidates); the other 79 bits are fixed to their correct values. The default guessed bit positions (canonical enumeration of the 85 bits) are `18, 57, 8, 64, 47, 23` (override with `--bits`).

Attack parameters from the paper's formula (1) with filtering strength `a = 30` and success rate `P_S = 0.9`:

| quantity | value |
|---|---|
| `N_v = ((Phi^-1(0.9)+Phi^-1(1-2^-30))/cor)^2` | 2308 |
| `N = N_v * 2^14.19` (pairs per attack) | 43130000 |
| key recovery statistic | observed correlation of the DL distinguisher |
| threshold `t = Phi^-1(1-2^-30)/sqrt(N_v)` | 0.1251 |

For every valid pair of the correct key guess, the statistic is additionally compared against a full two-round decryption under the real round keys (embedded sanity check).

## Build and run (100 attacks)

```bash
g++ -O2 -std=c++17 -o attack attack.cpp partial_dec.cpp tristate.cpp lea.cpp -pthread
./attack --trials 100 --seed 100 --threads 8
```

Each attack uses a fresh random master key. Running time: about 15 minutes on 8 threads; `--threads 1` also works but is proportionally slower.

## Expected results (as reported in the rebuttal)

- **90 of 100 attacks succeeded** (`summary: 90/100 trials succeeded`), i.e., the correct key guess passes the threshold `t`;
- **no wrong key guess survived** (`wrong_pass = 0` in all 100 attacks, 6300 wrong candidates in total);
- the correct candidate ranks first in all 100 attacks;
- the embedded sanity check reports `sanity_mismatch = 0` everywhere, i.e., the statistic of every valid pair computed by the partial decryption exactly matches the full two-round decryption.
