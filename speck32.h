#pragma once
#include <stdint.h>
#include <utility>
#include <stdlib.h>
#include <vector>
#include "rand_gen.h"
using namespace std;

typedef uint16_t word;
typedef pair<word, word> block;

#define WORD_SIZE 16
#define BLOCK_SIZE 32
#define ALPHA 7
#define BETA 2
#define MASK_VAL 0xffff
#define M 4
#define MAX_NR 50
#define RAND_WORD (global_random_generator.random_16bits())
#define RAND_BLOCK {RAND_WORD, RAND_WORD}
#define RAND_WORD_X(x) (x->random_16bits())
#define RAND_BLOCK_X(x) {RAND_WORD_X(x), RAND_WORD_X(x)}

struct linear_constraint {
    uint32_t xor_bit_pos[2];
    uint32_t num_bits;
    uint32_t xor_value;
};

struct data_subset {
    vector<block> ly0;
    vector<block> ly1;
};

struct simplified_data_subset_stage2 {
    block ly0[256];
    block ly1[256];
    uint32_t length;
};

word rol(const word& a, const uint32_t& b);
word ror(const word& a, const uint32_t& b);
void enc_one_round(const block& p, const word& k, block& c);
void dec_one_round(const block& c, const word& k, block& p);
void dec_one_round(const word& cl, const word& cr, const word& k, word& pl, word& pr);
void expand_key(const word mk[], word keys[], const uint32_t& nr);
void encrypt(const block& p, const word keys[], const uint32_t& nr, block& c);
void decrypt(const block& c, const word keys[], const uint32_t& nr, block& p);
bool check_testvector();

void generate_one_user_key(word user_key[M]);
void generate_one_user_key(word user_key[M], RandomGenerator *rand_engine);
block generate_one_plaintext();
block generate_one_plaintext(RandomGenerator *rand_engine);
void make_test_set(const uint32_t& n, const block& diff, const uint32_t& num_rounds, block c0[], block c1[], bool Y[]);
void make_target_diff_samples(const uint32_t& n, const block& diff, const uint32_t& num_rounds, block c0[], block c1[], bool whether_positive);
void collect_ciphertext_structure(const uint32_t& n, const uint32_t& attack_nr, const block p0[], const block p1[], block c0[], block c1[], const word user_round_keys[]);
block gen_rand_block_with_linear_constraint(const vector<linear_constraint>& constraints, RandomGenerator* rand_engine);
word gen_rand_word_with_linear_constraint(const vector<linear_constraint>& constraints, RandomGenerator* rand_engine);