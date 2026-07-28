#include "speck32.h"
#include "util.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <vector>
#include <thread>
#include <chrono>
#include <string>
#include <random>
#include <algorithm>
using namespace std;

#define TEST_TK false
#define ATTACK_THREAD_NUM 50

block input_diff;
uint32_t attack_nr;
vector<linear_constraint> guessed_k0_bits;
vector<linear_constraint> pseudoplaintext_constraints;
uint32_t N_stage1;
uint32_t N_stage2;
uint32_t N_stage3;
uint32_t N_stage4;
uint32_t N_rest;
word k0_space_size;
vector<block> constrainted_p0;
vector<block> constrainted_p1;
vector<word> k1_space_stage1;
vector<word> k2_space_stage1;
vector<word> k1_space_stage2;
vector<word> k2_space_stage2;
vector<word> k1_space_stage3;
vector<word> k2_space_stage3;
vector<word> k3_space_stage3;
vector<word> k1_space_stage4;
vector<word> k2_space_stage4;
vector<word> k3_space_stage4;
vector<word> k1_space_rest;
vector<word> k2_space_rest;
vector<word> k3_space_rest;
uint32_t partition_pos1_stage1;
uint32_t m1_stage1;
uint32_t partition_pos1_stage2;
uint32_t m1_stage2;
uint32_t partition_pos2_1_stage3;
uint32_t m2_1_stage3;
uint32_t partition_pos2_2_stage3;
uint32_t m2_2_stage3;
uint32_t partition_pos2_stage4;
uint32_t m2_stage4;
double c_stage1;
double c_stage2;
double c_stage3;
double c_stage4;

struct surviving_k1_k2_info {
    word k1;
    word k2;
    double cor;
};

struct surviving_kg_info {
    word k1;
    word k2;
    word k3;
    double cor;
};

block decrypt_one_round(const block& ct) {
    block pt = ct;
    pt.first = rol(pt.first - pt.second, ALPHA);
    pt.second = ror(pt.second ^ pt.first, BETA);
    return pt;
}

word extract_target_val(const word& x, const uint32_t& target_pos, const word& mask_val) {
    return (x >> target_pos) & mask_val;
}

uint32_t get_final_inner_product_stage1(const uint32_t& counter_index, const word& k2) {
    word example_l0 = ((counter_index & 0x1800) >> 2) | ((counter_index & 0x400) >> 10);
    word example_y0 = ((counter_index & 0x300) << 1) | ((counter_index & 0x80) >> 7);
    word example_l1 = ((counter_index & 0x60) << 4) | ((counter_index & 0x10) >> 4);
    word example_y1 = ((counter_index & 0xc) << 7) | ((counter_index & 0x2) >> 1);
    word dx = ((example_l0 ^ k2) - example_y0) ^ ((example_l1 ^ k2) - example_y1);
    return ((dx >> 1) ^ (dx >> 11) ^ counter_index) & 0x1;
}

uint32_t get_final_inner_product_stage2(const uint32_t& counter_index, const word& k2) {
    word example_l0 = ((counter_index & 0x1e0000) >> 10) | ((counter_index & 0x10000) >> 16);
    word example_y0 = ((counter_index & 0xf000) >> 5) | ((counter_index & 0x800) >> 11);
    word example_l1 = ((counter_index & 0x780)) | ((counter_index & 0x40) >> 6);
    word example_y1 = ((counter_index & 0x3c) << 5) | ((counter_index & 0x2) >> 1);
    word dx = ((example_l0 ^ k2) - example_y0) ^ ((example_l1 ^ k2) - example_y1);
    return ((dx >> 1) ^ (dx >> 11) ^ counter_index) & 0x1;
}

uint32_t get_final_inner_product_stage3(const uint32_t& counter_index, const word& k3) {
    word example_l0 = ((counter_index & 0x180) >> 2);
    word example_y0 = ((counter_index & 0x60));
    word example_l1 = ((counter_index & 0x18) << 2);
    word example_y1 = ((counter_index & 0x6) << 4);
    word dx = ((example_l0 ^ k3) - example_y0) ^ ((example_l1 ^ k3) - example_y1);
    return ((dx >> 7) ^ counter_index) & 0x1;
}

uint32_t get_final_inner_product_stage4(const uint32_t& counter_index, const word& k3) {
    word example_l0 = ((counter_index & 0x180));
    word example_y0 = ((counter_index & 0x60) << 2);
    word example_l1 = ((counter_index & 0x18) << 4);
    word example_y1 = ((counter_index & 0x6) << 6);
    word dx = ((example_l0 ^ k3) - example_y0) ^ ((example_l1 ^ k3) - example_y1);
    return ((dx >> 9) ^ counter_index) & 0x1;
}

uint32_t get_partition_index(const block& ly0, const block& ly1, const uint32_t& partition_pos, const uint32_t& m, const word& m_mask_val) {
    uint32_t index = extract_target_val(ly0.first, partition_pos, m_mask_val) ^ extract_target_val(ly0.second, partition_pos, m_mask_val);
    index = (index << m) | (extract_target_val(ly1.first, partition_pos, m_mask_val) ^ extract_target_val(ly1.second, partition_pos, m_mask_val));
    return index;
}

uint32_t get_partition_index_two_pos(const block& ly0, const block& ly1, const uint32_t& partition_pos1, const uint32_t& m1, const word& m1_mask_val, const uint32_t& partition_pos2, const uint32_t& m2, const word& m2_mask_val) {
    uint32_t index = extract_target_val(ly0.first, partition_pos1, m1_mask_val) ^ extract_target_val(ly0.second, partition_pos1, m1_mask_val);
    index = (index << m1) | (extract_target_val(ly1.first, partition_pos1, m1_mask_val) ^ extract_target_val(ly1.second, partition_pos1, m1_mask_val));
    index = (index << m2) | (extract_target_val(ly0.first, partition_pos2, m2_mask_val) ^ extract_target_val(ly0.second, partition_pos2, m2_mask_val));
    index = (index << m2) | (extract_target_val(ly1.first, partition_pos2, m2_mask_val) ^ extract_target_val(ly1.second, partition_pos2, m2_mask_val));
    return index;
}

void partition_data(const block ly0_vec[], const block ly1_vec[], const uint32_t& N, const uint32_t& partition_pos, const uint32_t& m, block partitioned_ly0_vec[], block partitioned_ly1_vec[], vector<uint32_t>& subset_index) {
    vector<data_subset> partitioned_data;
    partitioned_data.resize(1 << (2 * m));
    const word m_mask_val = (1 << m) - 1;
    for (uint32_t i = 0; i < N; i++) {
        const block& ly0 = ly0_vec[i];
        const block& ly1 = ly1_vec[i];
        uint32_t partition_index = get_partition_index(ly0, ly1, partition_pos, m, m_mask_val);
        partitioned_data[partition_index].ly0.push_back(ly0);
        partitioned_data[partition_index].ly1.push_back(ly1);
    }
    subset_index.clear();
    subset_index.push_back(0);
    uint32_t index = 0;
    for (uint32_t partition_index = 0; partition_index < partitioned_data.size(); partition_index++) {
        for (uint32_t i = 0; i < partitioned_data[partition_index].ly0.size(); i++) {
            partitioned_ly0_vec[index] = partitioned_data[partition_index].ly0[i];
            partitioned_ly1_vec[index] = partitioned_data[partition_index].ly1[i];
            index++;
        }
        subset_index.push_back(index);
    }
    assert(index == N);
}

void partition_data_with_two_pos(const block ly0_vec[], const block ly1_vec[], const uint32_t& N, const uint32_t& partition_pos1, const uint32_t& m1, const uint32_t& partition_pos2, const uint32_t& m2, block partitioned_ly0_vec[], block partitioned_ly1_vec[], vector<uint32_t>& subset_index) {
    vector<data_subset> partitioned_data;
    partitioned_data.resize(1 << (2 * m1 + 2 * m2));
    const word m1_mask_val = (1 << m1) - 1;
    const word m2_mask_val = (1 << m2) - 1;
    for (uint32_t i = 0; i < N; i++) {
        const block& ly0 = ly0_vec[i];
        const block& ly1 = ly1_vec[i];
        uint32_t partition_index = get_partition_index_two_pos(ly0, ly1, partition_pos1, m1, m1_mask_val, partition_pos2, m2, m2_mask_val);
        partitioned_data[partition_index].ly0.push_back(ly0);
        partitioned_data[partition_index].ly1.push_back(ly1);
    }
    subset_index.clear();
    subset_index.push_back(0);
    uint32_t index = 0;
    for (uint32_t partition_index = 0; partition_index < partitioned_data.size(); partition_index++) {
        for (uint32_t i = 0; i < partitioned_data[partition_index].ly0.size(); i++) {
            partitioned_ly0_vec[index] = partitioned_data[partition_index].ly0[i];
            partitioned_ly1_vec[index] = partitioned_data[partition_index].ly1[i];
            index++;
        }
        subset_index.push_back(index);
    }
    assert(index == N);
}

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

uint32_t calculate_counter_index_stage2(const block& pt0, const block& pt1) {
    // Counter index: pt0.first[10~7,0] || pt0.second[10~7,0] || pt1.first[10~7,0] || pt1.second[10~7,0] || xor((pt0.first^pt1.first)[1,11],(pt0.second^pt1.second)[0,1,2,7,8,11])
    //                20~16             || 15~11              || 10~6              || 5~1                || 0
    uint32_t index = pt0.first;
    index = ((index & 0x780) << 10) | ((index & 0x1) << 16);
    index |= ((pt0.second & 0x780) << 5) | ((pt0.second & 0x1) << 11);
    index |= ((pt1.first & 0x780)) | ((pt1.first & 0x1) << 6);
    index |= ((pt1.second & 0x780) >> 5) | ((pt1.second & 0x1) << 1);   
    word xor_l = pt0.first ^ pt1.first, xor_y = pt0.second ^ pt1.second;
    word xor_val = (xor_l >> 1) ^ (xor_l >> 11) ^ (xor_y) ^ (xor_y >> 1) ^ (xor_y >> 2) ^ (xor_y >> 7) ^ (xor_y >> 8) ^ (xor_y >> 11);
    index |= (xor_val & 0x1);
    return index;
}

uint32_t calculate_counter_index_stage3(const block& pt0, const block& pt1) {
    // Counter index: pt0.first[6~5] || pt0.second[6~5] || pt1.first[6~5] || pt1.second[6~5] || xor((pt0.first^pt1.first)[0,7],(pt0.second^pt1.second)[0,5,14])
    //                8~7            || 6~5             || 4~3            || 2~1             || 0
    uint32_t index = ((pt0.first & 0x60) << 2);
    index |= ((pt0.second & 0x60));
    index |= ((pt1.first & 0x60) >> 2);
    index |= ((pt1.second & 0x60) >> 4);
    word xor_l = pt0.first ^ pt1.first, xor_y = pt0.second ^ pt1.second;
    word xor_val = ((xor_l) ^ (xor_l >> 7) ^ (xor_y) ^ (xor_y >> 5) ^ (xor_y >> 14));
    index |= (xor_val & 0x1);
    return index;
}

uint32_t calculate_counter_index_stage4(const block& pt0, const block& pt1) {
    // Counter index: pt0.first[8~7] || pt0.second[8~7] || pt1.first[8~7] || pt1.second[8~7] || xor((pt0.first^pt1.first)[9],(pt0.second^pt1.second)[0,5,9,14])
    //                8~7            || 6~5             || 4~3            || 2~1             || 0
    uint32_t index = ((pt0.first & 0x180));
    index |= ((pt0.second & 0x180) >> 2);
    index |= ((pt1.first & 0x180) >> 4);
    index |= ((pt1.second & 0x180) >> 6);
    word xor_l = pt0.first ^ pt1.first, xor_y = pt0.second ^ pt1.second;
    word xor_val = ((xor_l >> 9) ^ (xor_y) ^ (xor_y >> 5) ^ (xor_y >> 9) ^ (xor_y >> 14));
    index |= (xor_val & 0x1);
    return index;
}

bool check_valid_counter_stage1(const uint32_t& counter_index, const word& k2) {
    // Check (k2[10~9] ^ counter_index[12~11]) == counter_index[9~8]
    word k2_part = k2 >> 9;
    if (((k2_part ^ (counter_index >> 11) ^ (counter_index >> 8)) & 0x3) == 0) return false;
    // Check (k2[10~9] ^ counter_index[6~5]) == counter_index[3~2]
    if (((k2_part ^ (counter_index >> 5) ^ (counter_index >> 2)) & 0x3) == 0) return false;
    return true;
}

bool check_valid_counter_stage2(const uint32_t& counter_index, const word& k2) {
    // Check (k2[10~7] ^ counter_index[20~17]) == counter_index[15~12]
    word k2_part = k2 >> 7;
    if (((k2_part ^ (counter_index >> 17) ^ (counter_index >> 12)) & 0xf) == 0) return false;
    // Check (k2[10~7] ^ counter_index[10~7]) == counter_index[5~2]
    if (((k2_part ^ (counter_index >> 7) ^ (counter_index >> 2)) & 0xf) == 0) return false;
    return true;
}

bool check_valid_counter_stage3(const uint32_t& counter_index, const word& k3) {
    // Check (k3[5~6] ^ counter_index[7~8]) == counter_index[5~6]
    word k3_part = k3 >> 5;
    if (((k3_part ^ (counter_index >> 7) ^ (counter_index >> 5)) & 0x3) == 0) return false;
    // Check (k3[5~6] ^ counter_index[3~4]) == counter_index[1~2]
    if (((k3_part ^ (counter_index >> 3) ^ (counter_index >> 1)) & 0x3) == 0) return false;
    return true;
}

bool check_valid_counter_stage4(const uint32_t& counter_index, const word& k3) {
    // Check (k3[7~8] ^ counter_index[7~8]) == counter_index[5~6]
    word k3_part = k3 >> 7;
    if (((k3_part ^ (counter_index >> 7) ^ (counter_index >> 5)) & 0x3) == 0) return false;
    // Check (k3[7~8] ^ counter_index[3~4]) == counter_index[1~2]
    if (((k3_part ^ (counter_index >> 3) ^ (counter_index >> 1)) & 0x3) == 0) return false;
    return true;
}

uint32_t count_valid_num(const vector<block>& p0, const vector<block>& p1, const word rk[]) {
    block tmp_c0, tmp_c1;
    uint32_t valid_num = 0;
    for (uint32_t i = 0; i < p0.size(); i++) {
        encrypt(p0[i], rk, 2, tmp_c0);
        encrypt(p1[i], rk, 2, tmp_c1);
        if ((tmp_c0.first ^ tmp_c1.first) == 0x0040 && (tmp_c0.second ^ tmp_c1.second) == 0x0000) {
            valid_num++;
        }
    }
    return valid_num;
}

class ThreadTask {
private:
    bool attack_with_one_structure_stage4(const block ly0_vec[], const block ly1_vec[], const word& k1_stage3, const word& k2_stage3, const word& k3_stage3, surviving_kg_info& returned_kg) {
        // Backdoor: only using the correct key
        if (TEST_TK) {
            thread_k1_space_stage4.clear();
            thread_k1_space_stage4.push_back(debug_tk1_stage4);
            thread_k2_space_stage4.clear();
            thread_k2_space_stage4.push_back(debug_tk2_stage4);
            thread_k3_space_stage4.clear();
            thread_k3_space_stage4.push_back(debug_tk3_stage4);
        }

        block *pseudo_ly0_vec = new block[N_stage4], *pseudo_ly1_vec = new block[N_stage4];
        uint32_t subset_num = 1 << (2 * m2_stage4);
        vector<uint32_t> subset_index;
        vector<uint32_t> counter(1 << 9); // 9-bit counter
        vector<uint32_t> valid_counter_index;
        word m2_mask_val = (1 << m2_stage4) - 1;
        bool kg_surviving_flag = false;
        // 1. Guess k1
        for (const word& k1_part : thread_k1_space_stage4) {
            word k1 = k1_part | k1_stage3;
            // 2. Decrypt ly0_vec and ly1_vec with k1
            for (uint32_t i = 0; i < N_stage4; i++) {
                pseudo_ly0_vec[i] = decrypt_one_round({ly0_vec[i].first ^ k1, ly0_vec[i].second});
                pseudo_ly1_vec[i] = decrypt_one_round({ly1_vec[i].first ^ k1, ly1_vec[i].second});
            }
            // 3. Partition data
            partition_data(pseudo_ly0_vec, pseudo_ly1_vec, N_stage4, partition_pos2_stage4, m2_stage4, pseudo_ly0_vec, pseudo_ly1_vec, subset_index);
            assert(subset_index.size() == subset_num + 1);
            // 4. Guess k2
            for (const word& k2_part : thread_k2_space_stage4) {
                word k2 = k2_part | k2_stage3;
                // Extract k2 partition vals
                word k2_partition_val = extract_target_val(k2, partition_pos2_stage4, m2_mask_val);
                // Reset counter
                for (uint32_t i = 0; i < counter.size(); i++) {
                    counter[i] = 0;
                }
                // Analyze each data subset
                for (uint32_t partition_index = 0; partition_index < subset_num; partition_index++) {
                    // Filter invalid subset
                    if (((partition_index ^ k2_partition_val) & m2_mask_val) == 0 || (((partition_index >> m2_stage4) ^ k2_partition_val) & m2_mask_val) == 0) continue;
                    // Decrypt valid data
                    for (uint32_t i = subset_index[partition_index]; i < subset_index[partition_index + 1]; i++) {
                        block pt0 = decrypt_one_round({pseudo_ly0_vec[i].first ^ k2, pseudo_ly0_vec[i].second});
                        block pt1 = decrypt_one_round({pseudo_ly1_vec[i].first ^ k2, pseudo_ly1_vec[i].second});
                        uint32_t counter_index = calculate_counter_index_stage4(pt0, pt1);
                        counter[counter_index]++;
                    }
                }
                // 5. Guess k3
                for (const word& k3_part : thread_k3_space_stage4) {
                    word k3 = k3_part | k3_stage3;
                    uint32_t T = 0;
                    uint32_t N = 0;
                    for (uint32_t counter_index = 0; counter_index < counter.size(); counter_index++) {
                        // Filter invalid data
                        if (!check_valid_counter_stage4(counter_index, k3)) {
                            continue;
                        }
                        T += get_final_inner_product_stage4(counter_index, k3) * counter[counter_index];
                        N += counter[counter_index];
                    }
                    double cor = 2.0 * T / N - 1.0;
                    if (cor > c_stage4) {
                        kg_surviving_flag = true;
                        debug_surviving_kg_num_stage4++;
                        if (k1 == debug_tk1_stage4 && k2 == debug_tk2_stage4 && k3 == debug_tk3_stage4) {
                            fprintf(thread_output_file, "[Debug] Found the correct key guess in stage4!\n");
                            debug_success_stage4 = true;
                        }
                        if (returned_kg.cor < cor) {
                            returned_kg.k1 = k1;
                            returned_kg.k2 = k2;
                            returned_kg.k3 = k3;
                            returned_kg.cor = cor;
                        }
                    }
                }
            }
        }
        delete[] pseudo_ly0_vec; delete[] pseudo_ly1_vec;
        return kg_surviving_flag;
    }

    bool attack_with_one_structure_stage3(const block ly0_vec[], const block ly1_vec[], const word& k1_stage2, const word& k2_stage2, surviving_kg_info& returned_kg) {
        if (TEST_TK) {
            // Backdoor: only using the correct key
            thread_k1_space_stage3.clear();
            thread_k1_space_stage3.push_back(debug_tk1_stage3);
            thread_k2_space_stage3.clear();
            thread_k2_space_stage3.push_back(debug_tk2_stage3);
            thread_k3_space_stage3.clear();
            thread_k3_space_stage3.push_back(debug_tk3_stage3);
        }

        vector<surviving_kg_info> surviving_k1_k2_k3_vec;
        block *pseudo_ly0_vec = new block[N_stage3], *pseudo_ly1_vec = new block[N_stage3];
        uint32_t subset_num = 1 << (2 * m2_1_stage3 + 2 * m2_2_stage3);
        vector<uint32_t> subset_index;
        vector<uint32_t> counter(1 << 9); // 9-bit counter
        vector<uint32_t> valid_counter_index;
        word m2_1_mask_val = (1 << m2_1_stage3) - 1;
        word m2_2_mask_val = (1 << m2_2_stage3) - 1;
        bool surviving_flag = false;
        // 1. Guess k1
        for (const word& k1_part : thread_k1_space_stage3) {
            word k1 = k1_part | k1_stage2;
            // 2. Decrypt ly0_vec and ly1_vec with k1
            for (uint32_t i = 0; i < N_stage3; i++) {
                pseudo_ly0_vec[i] = decrypt_one_round({ly0_vec[i].first ^ k1, ly0_vec[i].second});
                pseudo_ly1_vec[i] = decrypt_one_round({ly1_vec[i].first ^ k1, ly1_vec[i].second});
            }
            // 3. Partition data
            partition_data_with_two_pos(pseudo_ly0_vec, pseudo_ly1_vec, N_stage3, partition_pos2_1_stage3, m2_1_stage3, partition_pos2_2_stage3, m2_2_stage3, pseudo_ly0_vec, pseudo_ly1_vec, subset_index);
            assert(subset_index.size() == subset_num + 1);
            // 4. Guess k2
            for (const word& k2_part : thread_k2_space_stage3) {
                word k2 = k2_part | k2_stage2;
                // Extract k2 partition vals
                word k2_partition_val1 = extract_target_val(k2, partition_pos2_1_stage3, m2_1_mask_val);
                word k2_partition_val2 = extract_target_val(k2, partition_pos2_2_stage3, m2_2_mask_val);
                // Reset counter
                for (uint32_t i = 0; i < counter.size(); i++) {
                    counter[i] = 0;
                }
                // Analyze each data subset
                for (uint32_t partition_index = 0; partition_index < subset_num; partition_index++) {
                    // Filter invalid subset
                    uint32_t partition_index1 = partition_index >> (2 * m2_2_stage3);
                    uint32_t partition_index2 = partition_index;
                    if (((partition_index1 ^ k2_partition_val1) & m2_1_mask_val) == 0 || (((partition_index1 >> m2_1_stage3) ^ k2_partition_val1) & m2_1_mask_val) == 0) continue;
                    if (((partition_index2 ^ k2_partition_val2) & m2_2_mask_val) == 0 || (((partition_index2 >> m2_2_stage3) ^ k2_partition_val2) & m2_2_mask_val) == 0) continue;
                    // Decrypt valid data
                    for (uint32_t i = subset_index[partition_index]; i < subset_index[partition_index + 1]; i++) {
                        block pt0 = decrypt_one_round({pseudo_ly0_vec[i].first ^ k2, pseudo_ly0_vec[i].second});
                        block pt1 = decrypt_one_round({pseudo_ly1_vec[i].first ^ k2, pseudo_ly1_vec[i].second});
                        uint32_t counter_index = calculate_counter_index_stage3(pt0, pt1);
                        counter[counter_index]++;
                    }
                }
                // 5. Guess k3
                for (const word& k3 : thread_k3_space_stage3) {
                    uint32_t T = 0;
                    uint32_t N = 0;
                    for (uint32_t counter_index = 0; counter_index < counter.size(); counter_index++) {
                        // Filter invalid data
                        if (!check_valid_counter_stage3(counter_index, k3)) {
                            continue;
                        }
                        T += get_final_inner_product_stage3(counter_index, k3) * counter[counter_index];
                        N += counter[counter_index];
                    }
                    double cor = 1.0 - 2.0 * T / N;
                    if (cor > c_stage3) {
                        surviving_k1_k2_k3_vec.push_back({k1, k2, k3, cor});
                    }
                }
            }
        }
        delete[] pseudo_ly0_vec; delete[] pseudo_ly1_vec;
        if (surviving_k1_k2_k3_vec.size() == 0) {
            return false;
        }
        // Sort surviving (k1,k2,k3) pairs according to their correlation value
        sort(surviving_k1_k2_k3_vec.begin(), surviving_k1_k2_k3_vec.end(), [](const surviving_kg_info& a, const surviving_kg_info& b) {
            return a.cor > b.cor;
        });
        // Only traverse the top-two surviving (k1,k2,k3) pairs in stage3 attack
        uint32_t top_num = min((uint32_t)2, (uint32_t)surviving_k1_k2_k3_vec.size());
        debug_surviving_kg_num_stage3 += top_num;
        // 6. For each surviving (k1,k2,k3), conduct stage4 attack
        for (uint32_t i = 0; i < top_num; i++) {
            word k1 = surviving_k1_k2_k3_vec[i].k1;
            word k2 = surviving_k1_k2_k3_vec[i].k2;
            word k3 = surviving_k1_k2_k3_vec[i].k3;
            if (k1 == debug_tk1_stage3 && k2 == debug_tk2_stage3 && k3 == debug_tk3_stage3) {
                fprintf(thread_output_file, "[Debug] Found the correct key guess in stage3!\n");
                debug_success_stage3 = true;
            } else if (TEST_TK) {
                continue;
            }
            surviving_flag = attack_with_one_structure_stage4(ly0_vec, ly1_vec, k1, k2, k3, returned_kg);
            if (surviving_flag) {
                break;
            }
        }
        return surviving_flag;
    }

    bool attack_with_one_structure_stage2(const block ly0_vec[], const block ly1_vec[], const word& k1_stage1, const word& k2_stage1, surviving_kg_info& returned_kg) {
        vector<surviving_k1_k2_info> surviving_k1_k2_vec;
        if (TEST_TK) {
            // Backdoor: only guessing the correct key
            thread_k1_space_stage2.clear();
            thread_k1_space_stage2.push_back(debug_tk1_stage2);
            thread_k2_space_stage2.clear();
            thread_k2_space_stage2.push_back(debug_tk2_stage2);
        }

        // 1. Partition data in stage2
        bool surviving_flag = false;
        uint32_t subset_num = 1 << (2 * m1_stage2);
        vector<uint32_t> subset_index;
        block *partitioned_ly0_vec = new block[N_stage2], *partitioned_ly1_vec = new block[N_stage2];
        partition_data(ly0_vec, ly1_vec, N_stage2, partition_pos1_stage2, m1_stage2, partitioned_ly0_vec, partitioned_ly1_vec, subset_index);
        assert(subset_index.size() == subset_num + 1);
        // 2. Guess k1
        vector<uint32_t> counter(1 << 21); // 21-bit counter
        word m1_mask_val = (1 << m1_stage2) - 1;
        for (const word& k1_part : thread_k1_space_stage2) {
            word k1 = k1_stage1 | k1_part;
            // Extract k1 partition val
            word k1_partition_val = extract_target_val(k1, partition_pos1_stage2, m1_mask_val);
            // Reset counter
            for (uint32_t i = 0; i < counter.size(); i++) {
                counter[i] = 0;
            }
            // Analyze each data subset
            for (uint32_t partition_index1 = 0; partition_index1 < subset_num; partition_index1++) {
                // Filter invalid subset
                if (((partition_index1 ^ k1_partition_val) & m1_mask_val) == 0 || (((partition_index1 >> m1_stage2) ^ k1_partition_val) & m1_mask_val) == 0) continue;
                // Decrypt valid data
                for (uint32_t i = subset_index[partition_index1]; i < subset_index[partition_index1 + 1]; i++) {
                    block pt0 = decrypt_one_round({partitioned_ly0_vec[i].first ^ k1, partitioned_ly0_vec[i].second});
                    block pt1 = decrypt_one_round({partitioned_ly1_vec[i].first ^ k1, partitioned_ly1_vec[i].second});
                    uint32_t counter_index = calculate_counter_index_stage2(pt0, pt1);
                    counter[counter_index]++;
                }
            }
            // 3. Guess k2
            for (const word& k2_part : thread_k2_space_stage2) {
                word k2 = k2_part | k2_stage1;
                uint32_t T = 0;
                uint32_t N = 0;
                for (uint32_t counter_index = 0; counter_index < counter.size(); counter_index++) {
                    // Filter invalid data
                    if (!check_valid_counter_stage2(counter_index, k2)) {
                        continue;
                    }
                    T += get_final_inner_product_stage2(counter_index, k2) * counter[counter_index];
                    N += counter[counter_index];
                }
                // Calculate correlation and filter wrong key guess
                double cor = 1.0 - 2.0 * T / N;
                if (cor > c_stage2) {
                    // Save a surviving key guess
                    surviving_k1_k2_vec.push_back({k1, k2, cor});
                }
            }
        }
        delete[] partitioned_ly0_vec; delete[] partitioned_ly1_vec;
        if (surviving_k1_k2_vec.size() == 0) {
            return false;
        }
        // Sort surviving (k1,k2) pairs according to their correlation value
        sort(surviving_k1_k2_vec.begin(), surviving_k1_k2_vec.end(), [](const surviving_k1_k2_info& a, const surviving_k1_k2_info& b) {
            return a.cor > b.cor;
        });
        // Only traverse the top-two surviving (k1,k2) pairs in stage2 attack
        uint32_t top_num = min((uint32_t)2, (uint32_t)surviving_k1_k2_vec.size());
        // 4. For each surviving (k1,k2), conduct stage3 attack
        debug_surviving_kg_num_stage2 += top_num;
        for (uint32_t i = 0; i < top_num; i++) {
            word k1 = surviving_k1_k2_vec[i].k1;
            word k2 = surviving_k1_k2_vec[i].k2;
            if (k1 == debug_tk1_stage2 && k2 == debug_tk2_stage2) {
                fprintf(thread_output_file, "[Debug] Found the correct key guess in stage2!\n");
                debug_success_stage2 = true;
            } else if (TEST_TK) {
                continue;
            }
            surviving_flag = attack_with_one_structure_stage3(ly0_vec, ly1_vec, k1, k2, returned_kg);
            if (surviving_flag) {
                break;
            }
        }
        return surviving_flag;
    }

    bool attack_with_one_structure_stage1(const block ly0_vec[], const block ly1_vec[], surviving_kg_info& returned_kg) {
        if (TEST_TK) {
            // Backdoor: only guessing the correct key
            thread_k1_space_stage1.clear();
            thread_k1_space_stage1.push_back(debug_tk1_stage1);
            thread_k2_space_stage1.clear();
            thread_k2_space_stage1.push_back(debug_tk2_stage1);
        }

        bool surviving_flag = false;
        vector<surviving_k1_k2_info> surviving_k1_k2_vec;
        // 1. Partition data in stage1
        uint32_t subset_num = 1 << (2 * m1_stage1);
        vector<uint32_t> subset_index;
        block *partitioned_ly0_vec = new block[N_stage1], *partitioned_ly1_vec = new block[N_stage1];
        partition_data(ly0_vec, ly1_vec, N_stage1, partition_pos1_stage1, m1_stage1, partitioned_ly0_vec, partitioned_ly1_vec, subset_index);
        assert(subset_index.size() == subset_num + 1);
        // 2. Guess k1
        vector<uint32_t> counter(1 << 13); // 13-bit counter
        word m1_mask_val = (1 << m1_stage1) - 1;
        for (const word& k1 : thread_k1_space_stage1) {
            // Extract k1 partition val
            word k1_partition_val = extract_target_val(k1, partition_pos1_stage1, m1_mask_val);
            // Reset counter
            for (uint32_t i = 0; i < counter.size(); i++) {
                counter[i] = 0;
            }
            // Analyze each data subset
            for (uint32_t partition_index1 = 0; partition_index1 < subset_num; partition_index1++) {
                // Filter invalid subset
                if (((partition_index1 ^ k1_partition_val) & m1_mask_val) == 0 || (((partition_index1 >> m1_stage1) ^ k1_partition_val) & m1_mask_val) == 0) continue;
                // Decrypt valid data
                for (uint32_t i = subset_index[partition_index1]; i < subset_index[partition_index1 + 1]; i++) {
                    block pt0 = decrypt_one_round({partitioned_ly0_vec[i].first ^ k1, partitioned_ly0_vec[i].second});
                    block pt1 = decrypt_one_round({partitioned_ly1_vec[i].first ^ k1, partitioned_ly1_vec[i].second});
                    uint32_t counter_index = calculate_counter_index_stage1(pt0, pt1);
                    counter[counter_index]++;
                }
            }
            // 3. Guess k2
            for (const word& k2 : thread_k2_space_stage1) {
                uint32_t T = 0;
                uint32_t N = 0;
                for (uint32_t counter_index = 0; counter_index < counter.size(); counter_index++) {
                    // Filter invalid data
                    if (!check_valid_counter_stage1(counter_index, k2)) {
                        continue;
                    }
                    T += get_final_inner_product_stage1(counter_index, k2) * counter[counter_index];
                    N += counter[counter_index];
                }
                // Calculate correlation and filter wrong key guess
                double cor = 1.0 - 2.0 * T / N;
                if (cor > c_stage1) {
                    // Save a surviving key guess
                    surviving_k1_k2_vec.push_back({k1, k2, cor});
                    debug_surviving_kg_num_stage1++;
                }
            }
        }
        delete[] partitioned_ly0_vec; delete[] partitioned_ly1_vec;
        if (surviving_k1_k2_vec.size() == 0) {
            return false;
        }
        // Sort surviving (k1,k2) pairs according to their correlation value
        sort(surviving_k1_k2_vec.begin(), surviving_k1_k2_vec.end(), [](const surviving_k1_k2_info& a, const surviving_k1_k2_info& b) {
            return a.cor > b.cor;
        });
        // Only traverse the top-64 surviving (k1,k2) pairs in stage1 attack
        uint32_t top_num = min((uint32_t)64, (uint32_t)surviving_k1_k2_vec.size());
        // 4. For each surviving (k1,k2), conduct stage2 attack
        fprintf(thread_output_file, "[Debug] Stage1 found %lu surviving (k1,k2) pairs\n", surviving_k1_k2_vec.size());
        fflush(thread_output_file);
        for (uint32_t i = 0; i < top_num; i++) {
            word k1 = surviving_k1_k2_vec[i].k1;
            word k2 = surviving_k1_k2_vec[i].k2;
            if (k1 == debug_tk1_stage1 && k2 == debug_tk2_stage1) {
                fprintf(thread_output_file, "[Debug] Found the correct key guess in stage1!\n");
                debug_success_stage1 = true;
            } else if (TEST_TK) {
                continue;
            }
            surviving_flag = attack_with_one_structure_stage2(ly0_vec, ly1_vec, k1, k2, returned_kg);
            if (surviving_flag) {
                break;
            }
        }
        return surviving_flag;
    }
    void recover_complete_rk(word rk[], word l[]) {
        for (int i = attack_nr - 2; i >= attack_nr - 4; i--) {
            l[i] = rol((rk[i+1] ^ rol(rk[i], BETA) ^ i) - rk[i], ALPHA);
        }
        for (int i = attack_nr - 5; i >= 0; i--) {
            dec_one_round(l[i+3], rk[i+1], (word)i, l[i], rk[i]);
        }
    }
    bool attack_key_exhaustion(const block check_p[], const block check_c[], const word& guessed_k1, const word& guessed_k2, const word& guessed_k3, word mk[]) {
        block tmp_c;
        word rk[MAX_NR], l[MAX_NR];
        for (word k1_part : thread_k1_space_rest) {
            word k1 = guessed_k1 | k1_part;
            for (word k2_part : thread_k2_space_rest) {
                word k2 = guessed_k2 | k2_part;
                for (word k3_part : thread_k3_space_rest) {
                    word k3 = (guessed_k3 | k3_part) ^ ((k2 & 0x8000) >> 9);
                    for (uint32_t k4_part = 0; k4_part < (1u << 16); k4_part++) {
                        word k4 = k4_part;
                        rk[attack_nr - 1] = k1;
                        rk[attack_nr - 2] = k2;
                        rk[attack_nr - 3] = k3;
                        rk[attack_nr - 4] = k4;
                        recover_complete_rk(rk, l);
                        encrypt(check_p[0], rk, attack_nr, tmp_c);
                        if (tmp_c.first != check_c[0].first || tmp_c.second != check_c[0].second) {
                            continue;
                        }
                        bool check_flag = true;
                        for (uint32_t i = 1; i < N_rest; i++) {
                            encrypt(check_p[i], rk, attack_nr, tmp_c);
                            if (tmp_c.first != check_c[i].first || tmp_c.second != check_c[i].second) {
                                check_flag = false;
                                break;
                            }
                        }
                        if (check_flag) {
                            // Recover the complete master key
                            mk[0] = l[2];
                            mk[1] = l[1];
                            mk[2] = l[0];
                            mk[3] = rk[0];
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }
public:
    uint32_t attack_num;
    RandomGenerator *thread_random_generator = nullptr;
    FILE *thread_output_file;
    // Thread own settings
    vector<word> thread_k1_space_stage1;
    vector<word> thread_k2_space_stage1;
    vector<word> thread_k1_space_stage2;
    vector<word> thread_k2_space_stage2;
    vector<word> thread_k1_space_stage3;
    vector<word> thread_k2_space_stage3;
    vector<word> thread_k3_space_stage3;
    vector<word> thread_k1_space_stage4;
    vector<word> thread_k2_space_stage4;
    vector<word> thread_k3_space_stage4;
    vector<word> thread_k1_space_rest;
    vector<word> thread_k2_space_rest;
    vector<word> thread_k3_space_rest;
    word debug_tk0_value;
    word debug_tk1_stage1;
    word debug_tk2_stage1;
    word debug_tk1_stage2;
    word debug_tk2_stage2;
    word debug_tk1_stage3;
    word debug_tk2_stage3;
    word debug_tk3_stage3;
    word debug_tk1_stage4;
    word debug_tk2_stage4;
    word debug_tk3_stage4;
    bool debug_success_stage1;
    bool debug_success_stage2;
    bool debug_success_stage3;
    bool debug_success_stage4;
    uint64_t debug_surviving_kg_num_stage1;
    uint64_t debug_surviving_kg_num_stage2;
    uint64_t debug_surviving_kg_num_stage3;
    uint64_t debug_surviving_kg_num_stage4;
    uint64_t total_surviving_kg_num_stage1 = 0;
    uint64_t total_surviving_kg_num_stage2 = 0;
    uint64_t total_surviving_kg_num_stage3 = 0;
    uint64_t total_surviving_kg_num_stage4 = 0;
    uint64_t tk_survive_num_stage1 = 0;
    uint64_t tk_survive_num_stage2 = 0;
    uint64_t tk_survive_num_stage3 = 0;
    uint64_t tk_survive_num_stage4 = 0;
    uint64_t structure_consumption_sum = 0;
    uint64_t attack_time_ms_sum = 0;
    uint64_t kg_surviving_time = 0;
    uint64_t attack_success_time = 0;

    void key_recovery_attack() {
        word user_key[M], rk[MAX_NR], recovered_mk[M];
        block *c0 = new block[N_stage1], *c1 = new block[N_stage1];
        surviving_kg_info surviving_kg;
        bool kg_surviving;
        bool attack_res;
        vector<linear_constraint> kg0_constraints = guessed_k0_bits;
        for (uint32_t attack_index = 0; attack_index < attack_num; attack_index++) {
            fprintf(thread_output_file, "Attack index %u:\n", attack_index);
            generate_one_user_key(user_key, thread_random_generator);
            expand_key(user_key, rk, attack_nr);
            debug_tk0_value = extract_target_linear_value(rk[0], guessed_k0_bits);
            debug_tk1_stage1 = rk[attack_nr - 1] & 0x1f3f;
            debug_tk2_stage1 = rk[attack_nr - 2] & 0x601;
            debug_tk1_stage2 = rk[attack_nr - 1] & 0x1fbf;
            debug_tk2_stage2 = rk[attack_nr - 2] & 0x781;
            debug_tk1_stage3 = rk[attack_nr - 1] & 0x1fff;
            debug_tk2_stage3 = rk[attack_nr - 2] & 0x7783;
            debug_tk3_stage3 = (rk[attack_nr - 3] & 0x60) ^ ((rk[attack_nr - 2] & 0x8000) >> 9);
            debug_tk1_stage4 = rk[attack_nr - 1] & 0x7fff;
            debug_tk2_stage4 = rk[attack_nr - 2] & 0x778f;
            debug_tk3_stage4 = (rk[attack_nr - 3] & 0x1e0) ^ ((rk[attack_nr - 2] & 0x8000) >> 9);

            debug_surviving_kg_num_stage1 = 0;
            debug_surviving_kg_num_stage2 = 0;
            debug_surviving_kg_num_stage3 = 0;
            debug_surviving_kg_num_stage4 = 0;
            debug_success_stage1 = false;
            debug_success_stage2 = false;
            debug_success_stage3 = false;
            debug_success_stage4 = false;
            kg_surviving = false;
            surviving_kg.cor = 0;
            auto start_t = chrono::system_clock::now();
            fprintf(thread_output_file, "Tk0 value is %x\n", debug_tk0_value);
            // 1. Guess k0
            for (word kg0_value = 0; kg0_value < k0_space_size; kg0_value++) {
                fprintf(thread_output_file, "Trying kg0 value %x\n", kg0_value);
                // Generate corresponding guess of k0
                for (uint32_t i = 0; i < guessed_k0_bits.size(); i++) {
                    kg0_constraints[i].xor_value = (kg0_value >> (guessed_k0_bits.size() - 1 - i)) & 0x1;
                }
                word kg0 = gen_rand_word_with_linear_constraint(kg0_constraints, thread_random_generator);
                // 2. Generate pseudoplaintext pairs conforming to the linear constraints, and collect ciphertext pairs
                for (uint32_t i = 0; i < N_stage1; i++) {
                    dec_one_round(constrainted_p0[i], kg0, c0[i]);
                    dec_one_round(constrainted_p1[i], kg0, c1[i]);
                    encrypt(c0[i], rk, attack_nr, c0[i]);
                    encrypt(c1[i], rk, attack_nr, c1[i]);
                    c0[i].second = ror(c0[i].second ^ c0[i].first, BETA);
                    c1[i].second = ror(c1[i].second ^ c1[i].first, BETA);
                }
                // 3. Attack with one structure in stage1
                kg_surviving = attack_with_one_structure_stage1(c0, c1, surviving_kg);
                fflush(thread_output_file);
                if (kg_surviving) {
                    break;
                }
            }
            if (kg_surviving) {
                kg_surviving_time++;
                fprintf(thread_output_file, "A key guess is returned. Kg difference is (%x,%x,%x).\n", surviving_kg.k1 ^ debug_tk1_stage4, surviving_kg.k2 ^ debug_tk2_stage4, surviving_kg.k3 ^ debug_tk3_stage4);
                // Conduct key exhaustion to recover the complete key
                // First prepare some plaintext-ciphertext pairs to check the correctness of a key guess
                block *check_p = new block[N_rest], *check_c = new block[N_rest];
                for (uint32_t i = 0; i < N_rest; i++) {
                    check_p[i] = RAND_BLOCK_X(thread_random_generator);
                    encrypt(check_p[i], rk, attack_nr, check_c[i]);
                }
                // Then conduct key exhaustion
                if (TEST_TK) {
                    kg_surviving = false;
                } else {
                    kg_surviving = kg_surviving && attack_key_exhaustion(check_p, check_c, surviving_kg.k1, surviving_kg.k2, surviving_kg.k3, recovered_mk);
                }
                delete[] check_p; delete[] check_c;
                if (kg_surviving && recovered_mk[0] == user_key[0] && recovered_mk[1] == user_key[1] && recovered_mk[2] == user_key[2] && recovered_mk[3] == user_key[3]) {
                    attack_success_time++;
                    fprintf(thread_output_file, "Attack succeeded!\n");
                } else {
                    fprintf(thread_output_file, "Attack failed in key exhaustion stage!\n");
                }
            } else {
                fprintf(thread_output_file, "Attack failed because no key guess survived!\n");
            }
            auto end_t = chrono::system_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(end_t - start_t);
            attack_time_ms_sum += duration.count();
            tk_survive_num_stage1 += debug_success_stage1;
            tk_survive_num_stage2 += debug_success_stage2;
            tk_survive_num_stage3 += debug_success_stage3;
            tk_survive_num_stage4 += debug_success_stage4;
            total_surviving_kg_num_stage1 += debug_surviving_kg_num_stage1;
            total_surviving_kg_num_stage2 += debug_surviving_kg_num_stage2;
            total_surviving_kg_num_stage3 += debug_surviving_kg_num_stage3;
            total_surviving_kg_num_stage4 += debug_surviving_kg_num_stage4;
            fprintf(thread_output_file, "Kg surviving num in stage1: %u\n", debug_surviving_kg_num_stage1);
            fprintf(thread_output_file, "Kg surviving num in stage2: %u\n", debug_surviving_kg_num_stage2);
            fprintf(thread_output_file, "Kg surviving num in stage3: %u\n", debug_surviving_kg_num_stage3);
            fprintf(thread_output_file, "Kg surviving num in stage4: %u\n", debug_surviving_kg_num_stage4);
            fprintf(thread_output_file, "Time used: %.2f seconds\n", (double)duration.count() / 1000.0);
            fflush(thread_output_file);
        }
        fprintf(thread_output_file, "Kg surviving rate: %.2f. Attack success rate: %.2f\n", (double)kg_surviving_time / attack_num, (double)attack_success_time / attack_num);
        fprintf(thread_output_file, "Average surviving key guess number in stage1: %.2f\n", (double)total_surviving_kg_num_stage1 / attack_num);
        fprintf(thread_output_file, "Average surviving key guess number in stage2: %.2f\n", (double)total_surviving_kg_num_stage2 / attack_num);
        fprintf(thread_output_file, "Average surviving key guess number in stage3: %.2f\n", (double)total_surviving_kg_num_stage3 / attack_num);
        fprintf(thread_output_file, "Average surviving key guess number in stage4: %.2f\n", (double)total_surviving_kg_num_stage4 / attack_num);
        fprintf(thread_output_file, "Tk surviving rate in stage1: %.2f\n", (double)tk_survive_num_stage1 / attack_num);
        fprintf(thread_output_file, "Tk surviving rate in stage2: %.2f\n", (double)tk_survive_num_stage2 / attack_num);
        fprintf(thread_output_file, "Tk surviving rate in stage3: %.2f\n", (double)tk_survive_num_stage3 / attack_num);
        fprintf(thread_output_file, "Tk surviving rate in stage4: %.2f\n", (double)tk_survive_num_stage4 / attack_num);
        fprintf(thread_output_file, "Average structure consumption: %.2f\n", (double)structure_consumption_sum / attack_num);
        fprintf(thread_output_file, "Average time used: %.2f seconds\n", (double)attack_time_ms_sum / attack_num / 1000.0);
        fflush(thread_output_file);
        delete[] c0; delete[] c1;
    }
};

void set_attack_parameters() {
    // Basic attack parameters
    input_diff = {0x3800u, 0x10u};
    attack_nr = 13;

    // Linear constraints for pseudoplaintexts
    guessed_k0_bits = {
        {{11, 4}, 2, 0},
        {{12, 4}, 2, 0},
        {{5, 6}, 2, 0},
    };
    pseudoplaintext_constraints = {
        {{27, 4}, 2, 0},
        {{28, 4}, 2, 1},
        {{6, 5}, 2, 0},
    };

    // Data size
    N_stage1 = 145121238u;
    N_stage2 = 144510653u;
    N_stage3 = 78020165u;
    N_stage4 = 88045857u;
    N_rest = 10u;

    // Used data
    constrainted_p0.resize(N_stage1);
    constrainted_p1.resize(N_stage1);
    generate_plaintext_pairs_for_13r_attack(constrainted_p0.data(), constrainted_p1.data(), N_stage1);

    k0_space_size = 1u << guessed_k0_bits.size();

    // Guessed k1 in stage1: k1[0~5] || k1[8~12]
    k1_space_stage1.clear();
    for (word k1_part1 = 0; k1_part1 < (1 << 6); k1_part1++) {
        for (word k1_part2 = 0; k1_part2 < (1 << 5); k1_part2++) {
            word k1 = (k1_part2 << 8) | k1_part1;
            k1_space_stage1.push_back(k1);
        }
    }
    // Guessed k2 in stage1: k2[0] || k2[9~10]
    k2_space_stage1.clear();
    for (word k2_part1 = 0; k2_part1 < (1 << 1); k2_part1++) {
        for (word k2_part2 = 0; k2_part2 < (1 << 2); k2_part2++) {
            word k2 = (k2_part2 << 9) | k2_part1;
            k2_space_stage1.push_back(k2);
        }
    }

    // Guessed k1 in stage2: k1[7]
    k1_space_stage2.clear();
    for (word k1 = 0; k1 < (1 << 1); k1++) {
        k1_space_stage2.push_back(k1 << 7);
    }

    // Guessed k2 in stage2: k2[7~8]
    k2_space_stage2.clear();
    for (word k2 = 0; k2 < (1 << 2); k2++) {
        k2_space_stage2.push_back(k2 << 7);
    }

    // Guessed k1 in stage3: k1[6]
    k1_space_stage3.clear();
    for (word k1_part = 0; k1_part < (1 << 1); k1_part++) {
        word k1 = k1_part << 6;
        k1_space_stage3.push_back(k1);
    }

    // Guessed k2 in stage3: k2[1,12~14]
    k2_space_stage3.clear();
    for (word k2_part1 = 0; k2_part1 < (1 << 1); k2_part1++) {
        for (word k2_part2 = 0; k2_part2 < (1 << 3); k2_part2++) {
            word k2 = (k2_part1 << 1) | (k2_part2 << 12);
            k2_space_stage3.push_back(k2);
        }
    }

    // Guessed k3 in stage3: k3[5~6]
    k3_space_stage3.clear();
    for (word k3_part = 0; k3_part < (1 << 2); k3_part++) {
        word k3 = k3_part << 5;
        k3_space_stage3.push_back(k3);
    }

    // Guessed k1 in stage4: k1[13~14]
    k1_space_stage4.clear();
    for (word k1_part = 0; k1_part < (1 << 2); k1_part++) {
        word k1 = k1_part << 13;
        k1_space_stage4.push_back(k1);
    }

    // Guessed k2 in stage4: k2[2~3]
    k2_space_stage4.clear();
    for (word k2_part = 0; k2_part < (1 << 2); k2_part++) {
        word k2 = k2_part << 2;
        k2_space_stage4.push_back(k2);
    }

    // Guessed k3 in stage4: k3[7~8]
    k3_space_stage4.clear();
    for (word k3_part = 0; k3_part < (1 << 2); k3_part++) {
        word k3 = k3_part << 7;
        k3_space_stage4.push_back(k3);
    }

    // Guessed k1 in key exhaustion: k1[15]
    k1_space_rest.clear();
    for (word k1_part = 0; k1_part < (1 << 1); k1_part++) {
        word k1 = k1_part << 15;
        k1_space_rest.push_back(k1);
    }

    // Guessed k2 in key exhaustion: k2[4~6,11,15]
    k2_space_rest.clear();
    for (word k2_part1 = 0; k2_part1 < (1 << 3); k2_part1++) {
        for (word k2_part2 = 0; k2_part2 < (1 << 1); k2_part2++) {
            for (word k2_part3 = 0; k2_part3 < (1 << 1); k2_part3++) {
                word k2 = (k2_part1 << 4) | (k2_part2 << 11) | (k2_part3 << 15);
                k2_space_rest.push_back(k2);
            }
        }
    }

    // Guessed k3 in key exhaustion: k3[0~4,9~15]
    k3_space_rest.clear();
    for (word k3_part1 = 0; k3_part1 < (1 << 5); k3_part1++) {
        for (word k3_part2 = 0; k3_part2 < (1 << 7); k3_part2++) {
            word k3 = k3_part1 | (k3_part2 << 9);
            k3_space_rest.push_back(k3);
        }
    }

    // Set data partition parameters
    partition_pos1_stage1 = 8;
    m1_stage1 = 1;
    partition_pos1_stage2 = 7;
    m1_stage2 = 2;
    partition_pos2_1_stage3 = 7;
    m2_1_stage3 = 4;
    partition_pos2_2_stage3 = 12;
    m2_2_stage3 = 2;
    partition_pos2_stage4 = 7;
    m2_stage4 = 2;

    // Set correlation threshold
    c_stage1 = pow(2, -11.0344);
    c_stage2 = pow(2, -10.9853);
    c_stage3 = pow(2, -9.1934);
    c_stage4 = pow(2, -9.7786);
}

int main(int argc, char* argv[]) {
    assert(argc >= 2);
    uint32_t n = atoi(argv[1]);
    uint32_t output_init_index = 0;
    if (argc >= 3) {
        output_init_index = atoi(argv[2]);
    }
    random_device rd;
    global_random_generator.set_rand_seed(rd());
    printf("Preparing attack...\n");
    set_attack_parameters();
    uint32_t n_per_thread = (n + ATTACK_THREAD_NUM - 1) / ATTACK_THREAD_NUM;
    thread* thread_pool[ATTACK_THREAD_NUM];
    ThreadTask task_pool[ATTACK_THREAD_NUM];
    uint32_t tmp_n = n;
    for (int i = 0; i < ATTACK_THREAD_NUM; i++) {
        if (n_per_thread > tmp_n)
            task_pool[i].attack_num = tmp_n;
        else
            task_pool[i].attack_num = n_per_thread;
        tmp_n -= task_pool[i].attack_num;
        string thread_output_file_path = "./attack_records/key_recovery_attack_thread_" + to_string(i + output_init_index) + ".txt";
        task_pool[i].thread_output_file = fopen(thread_output_file_path.c_str(), "w");
        task_pool[i].thread_random_generator = new RandomGenerator(rd());
        task_pool[i].thread_k1_space_stage1 = k1_space_stage1;
        task_pool[i].thread_k2_space_stage1 = k2_space_stage1;
        task_pool[i].thread_k1_space_stage2 = k1_space_stage2;
        task_pool[i].thread_k2_space_stage2 = k2_space_stage2;
        task_pool[i].thread_k1_space_stage3 = k1_space_stage3;
        task_pool[i].thread_k2_space_stage3 = k2_space_stage3;
        task_pool[i].thread_k3_space_stage3 = k3_space_stage3;
        task_pool[i].thread_k1_space_stage4 = k1_space_stage4;
        task_pool[i].thread_k2_space_stage4 = k2_space_stage4;
        task_pool[i].thread_k3_space_stage4 = k3_space_stage4;
        task_pool[i].thread_k1_space_rest = k1_space_rest;
        task_pool[i].thread_k2_space_rest = k2_space_rest;
        task_pool[i].thread_k3_space_rest = k3_space_rest;
    }
    // generate threads
    for (int i = 0; i < ATTACK_THREAD_NUM; i++) {
        thread_pool[i] = new thread(&ThreadTask::key_recovery_attack, &task_pool[i]);
    }
    printf("Begin attack. Total attack number: %u. Total thread number: %d.\n", n, ATTACK_THREAD_NUM);
    // join threads
    for (int i = 0; i < ATTACK_THREAD_NUM; i++) {
        thread_pool[i]->join();
    }
    // summarize results
    printf("Attack Summary:\n");
    uint64_t total_attack_time_ms = 0;
    uint64_t total_structure_consumption = 0;
    uint64_t total_kg_surviving_time = 0;
    uint64_t total_attack_success_time = 0;
    uint64_t total_surviving_kg_num_stage1 = 0;
    uint64_t total_surviving_kg_num_stage2 = 0;
    uint64_t total_surviving_kg_num_stage3 = 0;
    uint64_t total_surviving_kg_num_stage4 = 0;
    uint64_t total_tk_survive_num_stage1 = 0;
    uint64_t total_tk_survive_num_stage2 = 0;
    uint64_t total_tk_survive_num_stage3 = 0;
    uint64_t total_tk_survive_num_stage4 = 0;
    for (int i = 0; i < ATTACK_THREAD_NUM; i++) {
        total_attack_time_ms += task_pool[i].attack_time_ms_sum;
        total_structure_consumption += task_pool[i].structure_consumption_sum;
        total_kg_surviving_time += task_pool[i].kg_surviving_time;
        total_attack_success_time += task_pool[i].attack_success_time;
        total_surviving_kg_num_stage1 += task_pool[i].total_surviving_kg_num_stage1;
        total_surviving_kg_num_stage2 += task_pool[i].total_surviving_kg_num_stage2;
        total_surviving_kg_num_stage3 += task_pool[i].total_surviving_kg_num_stage3;
        total_surviving_kg_num_stage4 += task_pool[i].total_surviving_kg_num_stage4;
        total_tk_survive_num_stage1 += task_pool[i].tk_survive_num_stage1;
        total_tk_survive_num_stage2 += task_pool[i].tk_survive_num_stage2;
        total_tk_survive_num_stage3 += task_pool[i].tk_survive_num_stage3;
        total_tk_survive_num_stage4 += task_pool[i].tk_survive_num_stage4;
        fclose(task_pool[i].thread_output_file);
        delete thread_pool[i];
        delete task_pool[i].thread_random_generator;
    }
    printf("Kg surviving rate: %.2f. Attack success rate: %.2f\n", (double)total_kg_surviving_time / n, (double)total_attack_success_time / n);
    printf("Average structure consumption: %.2f\n", (double)total_structure_consumption / n);
    printf("Average time used: %.2f seconds\n", (double)total_attack_time_ms / n / 1000.0);
    printf("Average surviving key guess number in stage1: %.2f\n", (double)total_surviving_kg_num_stage1 / n);
    printf("Average surviving key guess number in stage2: %.2f\n", (double)total_surviving_kg_num_stage2 / n);
    printf("Average surviving key guess number in stage3: %.2f\n", (double)total_surviving_kg_num_stage3 / n);
    printf("Average surviving key guess number in stage4: %.2f\n", (double)total_surviving_kg_num_stage4 / n);
    printf("Tk surviving rate in stage1: %.2f\n", (double)total_tk_survive_num_stage1 / n);
    printf("Tk surviving rate in stage2: %.2f\n", (double)total_tk_survive_num_stage2 / n);
    printf("Tk surviving rate in stage3: %.2f\n", (double)total_tk_survive_num_stage3 / n);
    printf("Tk surviving rate in stage4: %.2f\n", (double)total_tk_survive_num_stage4 / n);
    return 0;
}