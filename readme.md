# Verification Experiments for "Key-Recovery-Friendly Differential-Linear Cryptanalysis for ARX Ciphers"

This repository contains the code of three verification experiments:

1. **13-round Speck32/64 key-recovery attack** (`13r-speck32-attack/`): a full single-key differential-linear attack recovering the 64-bit master key of 13-round Speck32/64, with 500 recorded attack trials and a success-rate parsing script.
2. **Splitting multiple modular subtractions** (`splitting-modular-subtraction/`): measures the proportion of valid samples when recovering a borrow bit by simultaneously splitting two (resp. three) modular subtractions, matching Lemma 4 (resp. Lemma 5) of the paper.
3. **Two-round partial decryption on LEA** (`lea-partial-decryption/`): a reduced 7-round key-recovery attack (5-round differential-linear distinguisher + the identical 2-round partial decryption of Figure 6), verifying the partial decryption procedure of Section 6.1 in an attack setting.
