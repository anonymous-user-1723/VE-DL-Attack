#include <vector>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <chrono>
#include <algorithm>
#include <thread>
#include "speck32.h"
#include "util.h"
using namespace std;

#define ATTACK_THREAD_NUM 50

uint32_t calculate_counter_index_stage1(const block& pt0, const block& pt1) {
    // Counter index: pt0.first[10~9,0] || pt0.second[10~9,0] || pt1.first[10~9,0] || pt1.second[10~9,0] || xor((pt0.first^pt1.first)[1,11],(pt0.second^pt1.second)[0,1,2,7,8,11])
    //                12~10             || 9~7                || 6~4               || 3~1                || 0
    uint32_t index = pt0.first;
    index = ((index & 0x600) << 2) | ((index & 0x1) << 10);
    index |= ((pt0.second & 0x600) >> 1) | ((pt0.second & 0x1) << 7);
    index |= ((pt1.first & 0x600) >> 4) | ((pt1.first & 0x1) << 4);
    index |= ((pt1.second & 0x600) >> 7) | ((pt1.second & 0x1) << 1);   
    word xor_l = pt0.first ^ pt1.first, xor_y = pt0.second ^ pt1.second;
    word xor_val = (xor_l >> 1) ^ (xor_l >> 11) ^ (xor_y) ^ (xor_y >> 1) ^ (xor_y >> 2) ^ (xor_y >> 7) ^ (xor_y >> 8) ^ (xor_y >> 11);
    index |= (xor_val & 0x1);
    return index;
}

void decrypt_one_round(const block& ct, const word& kg, block& pt) {
    word &pl = pt.first = ct.first;
    word &pr = pt.second = ct.second;
    pl ^= kg;
    pl -= pr;
    pl = rol(pl, ALPHA);
    pr ^= pl;
    pr = ror(pr, BETA);
    // pt.first = rol(pt.first - pt.second, ALPHA);
    // pt.second = ror(pt.second ^ pt.first, BETA);
}

// Traverse all 2^31 available plaintext pairs
// input diff is (0x40, 0)
double estimate_DLA_cor(const uint32_t& t, const uint32_t& num_rounds) {
    vector<double> cor;
    word mk[M], rk[MAX_NR];
    uint32_t tmp;
    block p0, p1;
    block c0, c1;
    word dl, dr;
    word xor_val;
    for (uint32_t i = 0; i < t; i++) {
        generate_one_user_key(mk);
        expand_key(mk, rk, num_rounds);
        uint32_t total_num = 0;
        uint32_t count = 0;
        bool stop_flag = false;
        while (total_num < (1u << 28)) {
            p0 = RAND_BLOCK;
            p1.first = p0.first ^ 0x40u;
            p1.second = p0.second ^ 0x0u;
            encrypt(p0, rk, num_rounds, c0);
            encrypt(p1, rk, num_rounds, c1);
            dl = c0.first ^ c1.first;
            dr = c0.second ^ c1.second;
            // Filter a subset
            c0.first = ror(c0.first, ALPHA) + c0.second;
            c1.first = ror(c1.first, ALPHA) + c1.second;

            // Cor = 5.6
            // if ((((c0.first ^ c0.second) >> 1) & 0xf) == 0) continue;
            // if ((((c1.first ^ c1.second) >> 1) & 0xf) == 0) continue;
            // xor_val = (dl >> 5) ^ (dl >> 12) ^ (dl >> 13) ^ (dl >> 14) ^ (dr >> 5) ^ (dr >> 13) ^ (dr >> 14);

            // Cor = 8.72
            // if ((((c0.first ^ c0.second) >> 6) & 0x1) == 0) continue;
            // if ((((c1.first ^ c1.second) >> 6) & 0x1) == 0) continue;
            xor_val = (dl >> 7) ^ (dl >> 14) ^ (dr >> 5) ^ (dr >> 7) ^ (dr >> 14);

            // if ((((c0.first ^ c0.second) >> 8) & 0x1) == 0) continue;
            // if ((((c1.first ^ c1.second) >> 8) & 0x1) == 0) continue;
            // xor_val = (dl >> 0) ^ (dr >> 0) ^ (dr >> 5) ^ (dr >> 14);

            count += (uint32_t)(xor_val & 0x1u);
            total_num++;
        }
        // for (uint32_t tmp0 = 0; tmp0 < (1u << 22); tmp0++) { // 0~21
        //     for (uint32_t tmp1 = 0; tmp1 < (1u << 9); tmp1++) { // 23~31
        //         tmp = tmp0 | (tmp1 << 23);
        //         p0.first = (word)((tmp >> 16) & 0xffffu);
        //         p0.second = (word)(tmp & 0xffffu);
        //         p1.first = p0.first ^ 0x40u;
        //         p1.second = p0.second ^ 0x0u;
        //         encrypt(p0, rk, num_rounds, c0);
        //         encrypt(p1, rk, num_rounds, c1);
        //         dl = c0.first ^ c1.first;
        //         dr = c0.second ^ c1.second;
        //         // Filter a subset
        //         c0.first = ror(c0.first, ALPHA) + c0.second;
        //         c1.first = ror(c1.first, ALPHA) + c1.second;
        //         // if ((((c0.first ^ c0.second) >> 6) & 0x1) == 0) continue;
        //         // if ((((c1.first ^ c1.second) >> 6) & 0x1) == 0) continue;
        //         // Extract dl[2] ^ dl[8] ^ dr[0] ^ dr[2] ^ dr[7] ^ dr[8]
        //         // xor_val = (dl >> 2) ^ (dl >> 8) ^ (dr) ^ (dr >> 2) ^ (dr >> 7) ^ (dr >> 8);
        //         xor_val = (dl >> 7) ^ (dl >> 14) ^ (dr >> 5) ^ (dr >> 7) ^ (dr >> 14);
        //         count += (uint32_t)(xor_val & 0x1u);
        //         total_num++;
        //         if (total_num == (1u << 24)) {
        //             stop_flag = true;
        //             break;
        //         }
        //     }
        //     if (stop_flag) {
        //         break;
        //     }
        // }
        double correlation = 1 - 2.0 * (double)count / (double)total_num;
        cor.push_back(correlation);
        printf("Test one key index: %u, correlation is %f, absolute correlation is 2^(%f)\n", i, correlation, log2(fabs(correlation)));
    }
    sort(cor.begin(), cor.end(), [](const double& a, const double& b) {
        return fabs(a) < fabs(b);
    });
    double mean, stddev;
    cal_mean_std(cor.data(), cor.size(), mean, stddev);
    printf("Correlation mean: 2^(%f), stddev: 2^(%f)\n", log2(abs(mean)), log2(stddev));
    for (uint32_t i = 0; i < t; i++) {
        printf("Correlation %u: %f, absolute correlation is 2^(%f)\n", i, cor[i], log2(fabs(cor[i])));
    }
    double median_cor = (cor[t / 2] + cor[(t - 1) / 2]) / 2.0;
    return median_cor;
}

void test_conditional_diff_pr(const block& in_diff, const block& out_diff, const uint32_t& diff_nr, const uint32_t& num_samples, vector<linear_constraint>& linear_constraints) {
    uint32_t count = 0;
    word mk[M], rk[MAX_NR];
    for (uint32_t i = 0; i < num_samples; i++) {
        generate_one_user_key(mk);
        expand_key(mk, rk, diff_nr);
        block p0 = gen_rand_block_with_linear_constraint(linear_constraints, &global_random_generator);
        block p1 = {p0.first ^ in_diff.first, p0.second ^ in_diff.second};
        block c0, c1;
        encrypt(p0, rk, diff_nr, c0);
        encrypt(p1, rk, diff_nr, c1);
        if ((c0.first ^ c1.first == out_diff.first) && (c0.second ^ c1.second == out_diff.second)) {
            count++;
        }
    }
    printf("Conditional differential probability: 2^(%f)\n", log2((double)count / num_samples));
}

void test_encryption_efficiency(const uint32_t& structure_size, const uint32_t& num_structures, const word& kg_space) {
    block *c = new block[structure_size];
    for (uint32_t i = 0; i < structure_size; i++) {
        c[i] = RAND_BLOCK;
    }
    auto start_t = chrono::system_clock::now();
    for (uint32_t i = 0; i < num_structures; i++) {
        for (word kg = 0; kg < kg_space; kg++) {
            for (uint32_t j = 0; j < structure_size; j++) {
                decrypt_one_round(c[j], kg, c[j]);
                // c[j] = decrypt_one_round({c[j].first ^ kg, c[j].second});
                // dec_one_round(c[j], kg, c[j]);
            }
        }
    }
    auto end_t = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_t - start_t);
    double test_sec = (double)duration.count() / 1000.0;
    uint64_t total_enc_num = (uint64_t)structure_size * (uint64_t)num_structures * (uint64_t)kg_space;
    printf("Time used: %.2f seconds, encryption efficiency: 2^(%.2f) one-round encryptions/second\n", test_sec, log2((double)total_enc_num / test_sec));
    // Print a random value
    word check_val = 0;
    for (uint32_t i = 0; i < structure_size; i++) {
        check_val ^= c[i].first ^ c[i].second;
    }
    printf("Check value: %x\n", check_val);
    delete[] c;
}

void test_extract_counter_index_stage1(const uint32_t& structure_size, const uint32_t& num_structures) {
    block *p0 = new block[structure_size], *p1 = new block[structure_size];
    uint32_t counter[1u << 13] = {0}; // 13-bit counter
    for (uint32_t i = 0; i < structure_size; i++) {
        p0[i] = RAND_BLOCK;
        p1[i] = RAND_BLOCK;
    }
    auto start_t = chrono::system_clock::now();
    for (uint32_t i = 0; i < num_structures; i++) {
        for (uint32_t j = 0; j < structure_size; j++) {
            uint32_t index = calculate_counter_index_stage1(p0[j], p1[j]);
            counter[index]++;
        }
    }
    auto end_t = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_t - start_t);
    double test_sec = (double)duration.count() / 1000.0;
    uint64_t total_cal_num = (uint64_t)structure_size * (uint64_t)num_structures;
    printf("Time used: %.2f seconds, counter index extraction efficiency: 2^(%.2f) extractions/second\n", test_sec, log2((double)total_cal_num / test_sec));
    delete[] p0; delete[] p1;
    // Print a random value
    uint32_t check_val = 0;
    for (uint32_t i = 0; i < (1u << 13); i++) {
        check_val += counter[i];
    }
    printf("Check value: %u\n", check_val);
}

int main() {
    // Test conditinoal differential probability
    global_random_generator.set_rand_seed(time(NULL));
    // block in_diff = {0x3800u, 0x10u}, out_diff = {0x40u, 0x0u};
    // vector<linear_constraint> linear_constraints = {
    //     {{27, 4}, 2, 0},
    //     {{28, 4}, 2, 1},
    //     {{6, 5}, 2, 0},
    // };
    // test_conditional_diff_pr(in_diff, out_diff, 1, 1u << 20, linear_constraints);
    // Estimate DLA correlation
    // double median_cor = estimate_DLA_cor(100, 8);
    // printf("Estimated median DLA correlation over 100 keys: %f, absolute correlation is 2^(%f)\n", median_cor, log2(fabs(median_cor)));
    test_encryption_efficiency(1u << 27, 1u, 1u << 6);
    // test_extract_counter_index_stage1(1u << 27, 1u << 6);
    return 0;
}