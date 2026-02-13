#include <stdint.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "rand_gen.h"
#include "speck32.h"
using namespace std;

void cal_mean_std(const double data[], const uint32_t& length, double& mean, double& sigma) {
    double sum = 0;
    for (uint32_t i = 0; i < length; i++) sum += data[i];
    mean = sum / length;
    sum = 0;
    for (uint32_t i = 0; i < length; i++) {
        sum += pow(data[i] - mean, 2);
    }
    sigma = sqrt(sum / (length - 1));
}

// dl[2, 8], dr[0, 2, 7, 8]
word extract_linear_mask_ID1(const block& c0, const block& c1) {
    word l_mask = 0x104u, r_mask = 0x185u;
    word dl = l_mask & (c0.first ^ c1.first), dr = r_mask & (c0.second ^ c1.second);
    word res = (dl >> 2) ^ (dl >> 8) ^ (dr) ^ (dr >> 2) ^ (dr >> 7) ^ (dr >> 8);
    return res & 0x1u;
}

void generate_plaintext_structure(const block& diff, const vector<uint32_t>& neutral_bits, block p0[], block p1[]) {
    uint32_t structure_size = 1;
    p0[0] = RAND_BLOCK;
    for (auto x : neutral_bits) {
        word dl = 0, dr = 0;
        if (x < WORD_SIZE) {
            dr |= 1 << x;
        } else {
            dl |= 1 << (x - WORD_SIZE);
        }
        for (int i = 0, j = structure_size; i < structure_size; i++, j++) {
            p0[j].first = p0[i].first ^ dl;
            p0[j].second = p0[i].second ^ dr;
        }
        structure_size <<= 1;
    }
    for (int i = 0; i < structure_size; i++) {
        p1[i].first = p0[i].first ^ diff.first;
        p1[i].second = p0[i].second ^ diff.second;
    }
}

void generate_plaintext_structure(const block& diff, const vector<uint32_t>& neutral_bits, block p0[], block p1[], RandomGenerator* rng) {
    uint32_t structure_size = 1;
    p0[0] = RAND_BLOCK_X(rng);
    for (auto x : neutral_bits) {
        word dl = 0, dr = 0;
        if (x < WORD_SIZE) {
            dr |= 1 << x;
        } else {
            dl |= 1 << (x - WORD_SIZE);
        }
        for (int i = 0, j = structure_size; i < structure_size; i++, j++) {
            p0[j].first = p0[i].first ^ dl;
            p0[j].second = p0[i].second ^ dr;
        }
        structure_size <<= 1;
    }
    for (int i = 0; i < structure_size; i++) {
        p1[i].first = p0[i].first ^ diff.first;
        p1[i].second = p0[i].second ^ diff.second;
    }
}

word extract_target_linear_value(const word& x, const vector<linear_constraint>& xor_info) {
    word res = 0;
    for (auto& lc : xor_info) {
        word xor_val = 0;
        for (uint32_t i = 0; i < lc.num_bits; i++) {
            xor_val ^= (x >> lc.xor_bit_pos[i]);
        }
        res = (res << 1) | (xor_val & 0x1u);
    }
    return res;
}

// Generate plaintext pairs ((l,r), (l',r')), such that:
// l ^ l' = 0x3800, r ^ r' = 0x10
// l[11] ^ r[4] = 0, l[12] ^ r[4] = 1, r[5] ^ r[6] = 0
void generate_plaintext_pairs_for_13r_attack(block p0[], block p1[], const uint32_t& N) {
    assert(N <= (1u << 28));
    uint32_t *p0_uint32 = new uint32_t[1u << 28];
    uint32_t index = 0;
    for (uint32_t p0_0_3 = 0; p0_0_3 < (1u << 4); p0_0_3++) {
        for (uint32_t p0_7_26 = 0; p0_7_26 < (1u << 20); p0_7_26++) {
            for (uint32_t p0_29_31 = 0; p0_29_31 < (1u << 3); p0_29_31++) {
                for (uint32_t p0_5 = 0; p0_5 < 2; p0_5++) {
                    uint32_t p0_val = p0_0_3 | (p0_7_26 << 7) | (p0_29_31 << 29) | (p0_5 << 5) | (p0_5 << 6) | (1u << 28);
                    p0_uint32[index++] = p0_val;
                }
            }
        }
    }
    // Shuffle the generated plaintexts
    srand((unsigned)time(NULL));
    random_shuffle(p0_uint32, p0_uint32 + index);
    // Select the first N plaintexts
    for (uint32_t i = 0; i < N; i++) {
        p0[i].first = (word)((p0_uint32[i] >> 16) & 0xffffu);
        p0[i].second = (word)(p0_uint32[i] & 0xffffu);
        p1[i].first = p0[i].first ^ 0x3800u;
        p1[i].second = p0[i].second ^ 0x10u;
    }
    delete[] p0_uint32;
}