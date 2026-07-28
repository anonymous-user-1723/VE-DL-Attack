#include <stdio.h>
#include <time.h>
#include "speck32.h"

word rol(const word& a, const uint32_t& b) {
    return (a << b) | (a >> (WORD_SIZE - b));
}

word ror(const word& a, const uint32_t& b) {
    return (a >> b) | (a << (WORD_SIZE - b));
}

void enc_one_round(const block& p, const word& k, block& c) {
    word& cl = c.first = p.first;
    word& cr = c.second = p.second;
    cl = ror(cl, ALPHA);
    cl += cr;
    cl ^= k;
    cr = rol(cr, BETA);
    cr ^= cl;
}

void dec_one_round(const block& c, const word& k, block& p) {
    word& pl = p.first = c.first;
    word& pr = p.second = c.second;
    pr ^= pl;
    pr = ror(pr, BETA);
    pl ^= k;
    pl -= pr;
    pl = rol(pl, ALPHA);
}

void dec_one_round(const word& cl, const word& cr, const word& k, word& pl, word& pr) {
    pr = ror(cr ^ cl, BETA);
    pl = rol((cl ^ k) - pr, ALPHA);
}

void expand_key(const word mk[], word keys[], const uint32_t& nr) {
    keys[0] = mk[M - 1];
    word l[MAX_NR];
    block inner_state;
    for (int i = 0, j = M - 2; i < M - 1; i++, j--) {
        l[i] = mk[j];
    }
    for (int i = 0, j = M - 1, k = 1; i < nr - 1; i++, j++, k++) {
        enc_one_round(block(l[i], keys[i]), i, inner_state);
        l[j] = inner_state.first;
        keys[k] = inner_state.second;
    }
}

void encrypt(const block& p, const word keys[], const uint32_t& nr, block& c) {
    c = p;
    for (int i = 0; i < nr; i++) {
        enc_one_round(c, keys[i], c);
    }
}

void decrypt(const block& c, const word keys[], const uint32_t& nr, block& p) {
    p = c;
    for (int i = nr - 1; i >= 0; i--) {
        dec_one_round(p, keys[i], p);
    }
}

bool check_testvector() {
    word mk[4] = {0x1918, 0x1110, 0x0908, 0x0100};
    block p = {0x6574, 0x694c};
    word keys[22];
    expand_key(mk, keys, 22);
    block c;
    encrypt(p, keys, 22, c);
    if (c == block(0xa868, 0x42f2)) {
        printf("Testvector verified.\n");
        return true;
    } else {
        printf("Testvector not verified.\n");
        return false;
    }
}

void generate_one_user_key(word user_key[M]) {
    for (int i = 0; i < M; i++) {
        user_key[i] = RAND_WORD;
    }
}

void generate_one_user_key(word user_key[M], RandomGenerator *rand_engine) {
    for (int i = 0; i < M; i++) {
        user_key[i] = RAND_WORD_X(rand_engine);
    }
}

block generate_one_plaintext() {
    return RAND_BLOCK;
}

block generate_one_plaintext(RandomGenerator *rand_engine) {
    return RAND_BLOCK_X(rand_engine);
}

void make_test_set(const uint32_t& n, const block& diff, const uint32_t& num_rounds, block c0[], block c1[], bool Y[]) {
    word ks[50];
    word mk[M];
    block p0, p1;
    for (uint32_t i = 0; i < n; i++) {
        generate_one_user_key(mk);
        expand_key(mk, ks, num_rounds);
        p0 = generate_one_plaintext();
        Y[i] = RAND_WORD & 1;
        if (Y[i]) {
            p1.first = p0.first ^ diff.first;
            p1.second = p0.second ^ diff.second;
        } else {
            p1 = generate_one_plaintext();
        }
        encrypt(p0, ks, num_rounds, c0[i]);
        encrypt(p1, ks, num_rounds, c1[i]);
    }
}

void make_target_diff_samples(const uint32_t& n, const block& diff, const uint32_t& num_rounds, block c0[], block c1[], bool whether_positive) {
    word ks[50];
    word mk[M];
    block p0, p1;
    for (uint32_t i = 0; i < n; i++) {
        generate_one_user_key(mk);
        expand_key(mk, ks, num_rounds);
        p0 = generate_one_plaintext();
        if (whether_positive) {
            p1.first = p0.first ^ diff.first;
            p1.second = p0.second ^ diff.second;
        } else {
            p1 = generate_one_plaintext();
        }
        encrypt(p0, ks, num_rounds, c0[i]);
        encrypt(p1, ks, num_rounds, c1[i]);
    }
}

void collect_ciphertext_structure(const uint32_t& n, const uint32_t& attack_nr, const block p0[], const block p1[], block c0[], block c1[], const word user_round_keys[]) {
    for (uint32_t i = 0; i < n; i++) {
        encrypt(p0[i], user_round_keys, attack_nr, c0[i]);
        encrypt(p1[i], user_round_keys, attack_nr, c1[i]);
    }
}

block gen_rand_block_with_linear_constraint(const vector<linear_constraint>& constraints, RandomGenerator* rand_engine) {
    uint32_t rand_int = RAND_WORD_X(rand_engine);
    rand_int = (rand_int << 16) | RAND_WORD_X(rand_engine);
    uint32_t tmp;
    for (auto x : constraints) {
        tmp = x.xor_value;
        for (int i = 0; i < x.num_bits; i++) tmp ^= (rand_int >> x.xor_bit_pos[i]);
        rand_int ^= (tmp & 0x1) << x.xor_bit_pos[0];
    }
    return {uint16_t(rand_int >> 16), uint16_t(rand_int & 0xFFFFu)};
}

word gen_rand_word_with_linear_constraint(const vector<linear_constraint>& constraints, RandomGenerator* rand_engine) {
    word rand_word = RAND_WORD_X(rand_engine);
    word tmp;
    for (auto x : constraints) {
        tmp = x.xor_value;
        for (int i = 0; i < x.num_bits; i++) tmp ^= (rand_word >> x.xor_bit_pos[i]);
        rand_word ^= (tmp & 0x1) << x.xor_bit_pos[0];
    }
    return rand_word;
}