# Improved Differential-Linear Attack on 13-round Speck32/64

This code repository implements the differential-linear attack on 13-round Speck32/64 proposed in the paper "Towards Key-Recovery-Friendly Differential-Linear Cryptanalysis for ARX Ciphers". The objective of this attack code is to recover the complete 64-bit master key of Speck32/64 in the single-key setting.

The code is written in C++ and is compiled via g++. Run `make` to produce the executable file `key_recovery_attack`.

## Running the key recovery attack

The attack code is implemented in `key_recovery_attack.cpp`, which is a multi-thread program using the C++ lib \<thread\> .

The macro definition `ATTACK_THREAD_NUM` is the number of threads used to conduct attacks.

To conduct attacks, run the following command: `./key_recovery_attack [n]`, where the parameter `n` is the total number of attack trials. The program will assign these `n` attacks evenly to `ATTACK_THREAD_NUM` threads.

The attack record of each thread will be saved in the folder `./attack_records/`, which one has to create before running the attack program.

## Showing existing attack results

We have run the 13-round attack 500 times using 50 thread. The attack records are provided in the folder `./attack_records_500times/`.

One can run the script `./parse_results.py` to show the success rate of our 500-time attack experiment.
