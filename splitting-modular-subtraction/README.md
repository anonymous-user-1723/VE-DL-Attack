# Splitting multiple modular subtractions

For a chain of two (resp. three) consecutive modular subtractions as in Figure 4 of the paper, we measure the proportion of valid samples (those for which the target borrow bit can be determined) when the subtractions are split simultaneously, and compare it with the theoretical value of Lemma 4 (resp. Lemma 5). A sample is counted as valid only when the borrow bit is actually recovered by the attack-style (constructive) CDB procedure; each constructively valid sample is additionally sanity checked by key-perturbation invariance.

## Build and run

```bash
g++ -O2 -std=c++17 -o verify_two_modsub verify_two_modsub.cpp
g++ -O2 -std=c++17 -o verify_three_modsub verify_three_modsub.cpp

./verify_two_modsub --T 10000 --seed 123
./verify_three_modsub --T 10000 --seed 123
```

Both programs run a small grid of splitting parameters with `N = 100000` random samples per cell.

## Expected results (the rows reported in the rebuttal)

- two subtractions, `t=1, m=4`: `p_hat = 0.8111` vs. theory `(1-(1+2^-1*4)*2^-4) = 0.8125`;
- three subtractions, `s=1, t=1, m=4`: `p_hat = 0.6553` vs. theory `1-(1+4*2^-1*(1+(4+2*1-1)*2^-2))*2^-4 = 0.65625`.

All grid cells pass with `abs_err` well below the threshold. Lemmas 2, 4 and 5 actually state lower bounds on the expected number of valid samples: extra valid samples are also found in the experiments, making our estimates conservative (this does not affect our attacks; the `p_obo` column gives the one-by-one splitting baseline for comparison).
