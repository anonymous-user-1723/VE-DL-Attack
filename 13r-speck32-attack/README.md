# 13-round Speck32/64 key-recovery attack

This directory implements the differential-linear attack on 13-round Speck32/64 proposed in the paper. The objective is to recover the complete 64-bit master key of Speck32/64 in the single-key setting.

The code is written in C++ and is compiled via `g++`. Run `make` to produce the executable `key_recovery_attack`.

## Running the key recovery attack

The attack is implemented in `key_recovery_attack.cpp` as a multi-thread program using C++ `<thread>`.

The macro `ATTACK_THREAD_NUM` is the number of threads used to conduct attacks.

To run attacks: `./key_recovery_attack [n]`, where `n` is the total number of attack trials. The program assigns these `n` attacks evenly to `ATTACK_THREAD_NUM` threads.

The attack record of each thread is saved in `./attack_records/`, which must be created before running the attack program.

## Showing existing attack results

We have run the 13-round attack 500 times using 50 threads. The attack records are provided in `./attack_records_500times/`.

Run `./parse_results.py` to show the success rate of the 500-trial experiment.
