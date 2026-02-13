#pragma once
#include <stdint.h>
#include <utility>
#include <stdlib.h>
#include <vector>
#include "speck32.h"
#include "rand_gen.h"
using namespace std;

void cal_mean_std(const double data[], const uint32_t& length, double& mean, double& sigma);
word extract_linear_mask_ID1(const block& c0, const block& c1);
void generate_plaintext_structure(const block& diff, const vector<uint32_t>& neutral_bits, block p0[], block p1[], RandomGenerator* rng);
void generate_plaintext_structure(const block& diff, const vector<uint32_t>& neutral_bits, block p0[], block p1[]);
word extract_target_linear_value(const word& x, const vector<linear_constraint>& xor_info);
void generate_plaintext_pairs_for_13r_attack(block p0[], block p1[], const uint32_t& N);