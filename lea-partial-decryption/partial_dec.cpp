#include "partial_dec.h"

#include <initializer_list>

// bits a..b inclusive
static word range_mask(int a, int b) {
    const word hi = (b == 31) ? 0xFFFFFFFFu : ((word{1} << (b + 1)) - 1u);
    return hi & ~((word{1} << a) - 1u);
}

static word bits_mask(std::initializer_list<int> bits) {
    word m = 0;
    for (int b : bits) m |= word{1} << b;
    return m;
}

const word GAMMA_MASK[4] = {bits_mask({5, 6, 9, 18}), bits_mask({3, 4, 24, 27}),
                            bits_mask({6, 29}), bits_mask({0, 29})};

const word M_VY1_18 = range_mask(0, 8) | range_mask(14, 28);
const word M_VX1_18 = bits_mask({0, 21, 22, 23});
const word M_VY2_18 = range_mask(0, 8) | range_mask(15, 17) | range_mask(20, 28);
const word M_VX2_18 = bits_mask({24, 25});
const word M_VY3_18 = range_mask(0, 5) | bits_mask({8, 17}) | range_mask(21, 28);
const word M_VY1_17 = range_mask(1, 5) | range_mask(22, 28);
const word M_VY2_17 = bits_mask({5}) | range_mask(26, 28);
const word M_VY3_17 = range_mask(27, 28);

KeyBits85 extract_keybits85(const word Ka[6], const word Kb[6]) {
    KeyBits85 kb;
    kb.rk0_18 = Kb[0] & M_VY1_18;
    kb.rk1_18 = Kb[1] & M_VX1_18;
    kb.rk1x2_18 = (Kb[1] ^ Kb[2]) & M_VY2_18;
    kb.rk3_18 = Kb[3] & M_VX2_18;
    kb.rk3x4_18 = (Kb[3] ^ Kb[4]) & M_VY3_18;
    kb.rk5x0 = (Kb[5] ^ Ka[0]) & M_VY1_17;
    kb.rk1x2_17 = (Ka[1] ^ Ka[2]) & M_VY2_17;
    kb.rk3x4_17 = (Ka[3] ^ Ka[4]) & M_VY3_17;
    return kb;
}

void kb85_randomize_unguessed(KeyBits85& kb, const word rv[8]) {
    kb.rk0_18 = (kb.rk0_18 & M_VY1_18) | (rv[0] & ~M_VY1_18);
    kb.rk1_18 = (kb.rk1_18 & M_VX1_18) | (rv[1] & ~M_VX1_18);
    kb.rk1x2_18 = (kb.rk1x2_18 & M_VY2_18) | (rv[2] & ~M_VY2_18);
    kb.rk3_18 = (kb.rk3_18 & M_VX2_18) | (rv[3] & ~M_VX2_18);
    kb.rk3x4_18 = (kb.rk3x4_18 & M_VY3_18) | (rv[4] & ~M_VY3_18);
    kb.rk5x0 = (kb.rk5x0 & M_VY1_17) | (rv[5] & ~M_VY1_17);
    kb.rk1x2_17 = (kb.rk1x2_17 & M_VY2_17) | (rv[6] & ~M_VY2_17);
    kb.rk3x4_17 = (kb.rk3x4_17 & M_VY3_17) | (rv[7] & ~M_VY3_17);
}

int gamma_parity(const block& d) {
    const word t = (d.x0 & GAMMA_MASK[0]) ^ (d.x1 & GAMMA_MASK[1]) ^
                   (d.x2 & GAMMA_MASK[2]) ^ (d.x3 & GAMMA_MASK[3]);
    return __builtin_popcount(t) & 1;
}

// ---- canonical 85-bit indexing ----

static const word FIELD_MASKS[8] = {M_VY1_18, M_VX1_18, M_VY2_18, M_VX2_18,
                                    M_VY3_18, M_VY1_17, M_VY2_17, M_VY3_17};

static word field_value(const KeyBits85& kb, int f) {
    switch (f) {
        case 0: return kb.rk0_18;
        case 1: return kb.rk1_18;
        case 2: return kb.rk1x2_18;
        case 3: return kb.rk3_18;
        case 4: return kb.rk3x4_18;
        case 5: return kb.rk5x0;
        case 6: return kb.rk1x2_17;
        default: return kb.rk3x4_17;
    }
}

static void set_field_value(KeyBits85& kb, int f, word v) {
    switch (f) {
        case 0: kb.rk0_18 = v; break;
        case 1: kb.rk1_18 = v; break;
        case 2: kb.rk1x2_18 = v; break;
        case 3: kb.rk3_18 = v; break;
        case 4: kb.rk3x4_18 = v; break;
        case 5: kb.rk5x0 = v; break;
        case 6: kb.rk1x2_17 = v; break;
        default: kb.rk3x4_17 = v; break;
    }
}

int kb85_get_bit(const KeyBits85& kb, int idx) {
    for (int f = 0; f < 8; ++f) {
        const word m = FIELD_MASKS[f];
        const int nb = __builtin_popcount(m);
        if (idx < nb) {
            for (int b = 0; b < 32; ++b) {
                if ((m >> b) & 1u) {
                    if (idx == 0) return (field_value(kb, f) >> b) & 1u;
                    --idx;
                }
            }
        }
        idx -= nb;
    }
    return 0;  // idx out of range
}

void kb85_set_bit(KeyBits85& kb, int idx, int v) {
    for (int f = 0; f < 8; ++f) {
        const word m = FIELD_MASKS[f];
        const int nb = __builtin_popcount(m);
        if (idx < nb) {
            for (int b = 0; b < 32; ++b) {
                if ((m >> b) & 1u) {
                    if (idx == 0) {
                        word fv = field_value(kb, f);
                        if (v) fv |= word{1} << b;
                        else fv &= ~(word{1} << b);
                        set_field_value(kb, f, fv);
                        return;
                    }
                    --idx;
                }
            }
        }
        idx -= nb;
    }
}

static int borrow_known(const SubResult& s, int bit) {
    const word m = word{1} << bit;
    return ((s.bk0 & m) && (s.bk1 & m)) ? 1 : 0;
}

DecResult partial_dec_2r(const KeyBits85& kb, const block& C, const block& Cp,
                         bool early_term) {
    DecResult res;
    res.valid = false;
    res.stat_computable = false;
    res.stat = 0;
    res.borrow_ok = 0;

    // Free values from the ciphertext (Figure 6 caption).
    const TWord x0_18 = tw_known(C.x3, Cp.x3);
    const TWord z1_18 = tw_known((C.x0 >> 9) | (C.x0 << 23), (Cp.x0 >> 9) | (Cp.x0 << 23));
    const TWord z2_18 = tw_known((C.x1 << 5) | (C.x1 >> 27), (Cp.x1 << 5) | (Cp.x1 >> 27));
    const TWord z3_18 = tw_known((C.x2 << 3) | (C.x2 >> 29), (Cp.x2 << 3) | (Cp.x2 >> 29));

    // ---- round 18: x^19 -> x^18 ----
    const TWord y1_18 = tw_xor_key(x0_18, kb.rk0_18, M_VY1_18);
    const SubResult s1 = tw_sub(z1_18, y1_18);  // v1_18 = z1_18 - y1_18
    const TWord x1_18 = tw_xor_key(s1.v, kb.rk1_18, M_VX1_18);
    const TWord y2_18 = tw_xor_key(s1.v, kb.rk1x2_18, M_VY2_18);
    const SubResult s2 = tw_sub(z2_18, y2_18);  // v2_18
    res.borrow_ok |= static_cast<uint8_t>(borrow_known(s2, 17) << B2_18_17);
    if (early_term && !(res.borrow_ok & (1u << B2_18_17))) return res;

    const TWord x2_18 = tw_xor_key(s2.v, kb.rk3_18, M_VX2_18);
    const TWord y3_18 = tw_xor_key(s2.v, kb.rk3x4_18, M_VY3_18);
    const SubResult s3 = tw_sub(z3_18, y3_18);  // v3_18
    res.borrow_ok |= static_cast<uint8_t>(borrow_known(s3, 9) << B3_18_9);
    res.borrow_ok |= static_cast<uint8_t>(borrow_known(s3, 18) << B3_18_18);
    if (early_term &&
        (res.borrow_ok & ((1u << B3_18_9) | (1u << B3_18_18))) !=
            ((1u << B3_18_9) | (1u << B3_18_18)))
        return res;

    // x0^17 = x3^18 = v3_18 ^ rk5^18; rk5^18 alone is not guessed, so all
    // values are unknown but the difference is preserved.
    const TWord x0_17 = tw_xor_key(s3.v, 0, 0);
    const TWord z1_17 = tw_ror(x0_18, 9);
    const TWord z2_17 = tw_rol(x1_18, 5);
    const TWord z3_17 = tw_rol(x2_18, 3);

    // ---- round 17: x^18 -> x^17 ----
    const TWord y1_17 = tw_xor_key(s3.v, kb.rk5x0, M_VY1_17);
    const SubResult s4 = tw_sub(z1_17, y1_17);  // v1_17
    const TWord y2_17 = tw_xor_key(s4.v, kb.rk1x2_17, M_VY2_17);
    const SubResult s5 = tw_sub(z2_17, y2_17);  // v2_17
    const TWord y3_17 = tw_xor_key(s5.v, kb.rk3x4_17, M_VY3_17);
    const SubResult s6 = tw_sub(z3_17, y3_17);  // v3_17

    // Output words (blue: only the difference is needed). The remaining
    // un-guessed keys rk1^17, rk3^17, rk5^17 do not affect the difference.
    const TWord X[4] = {x0_17, s4.v, s5.v, s6.v};

    res.borrow_ok |= static_cast<uint8_t>(borrow_known(s5, 6) << B2_17_6);
    res.borrow_ok |= static_cast<uint8_t>(borrow_known(s4, 3) << B1_17_3);
    res.borrow_ok |= static_cast<uint8_t>(borrow_known(s6, 29) << B3_17_29);
    res.borrow_ok |= static_cast<uint8_t>(borrow_known(s4, 24) << B1_17_24);
    res.valid = (res.borrow_ok == 0x7F);
    if (!res.valid) return res;

    res.stat_computable = true;
    word diff_masked = 0;
    for (int w = 0; w < 4; ++w) {
        if ((X[w].dk & GAMMA_MASK[w]) != GAMMA_MASK[w]) res.stat_computable = false;
        diff_masked ^= X[w].dv & GAMMA_MASK[w];
    }
    res.stat = __builtin_popcount(diff_masked) & 1;
    return res;
}
