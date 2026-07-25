#pragma once
#include "lea.h"
#include "tristate.h"

/**
 * Two-round partial decryption of LEA (paper Section 6.1, Figure 6):
 * from the ciphertext pair (x^19, x'^19) back to x^17, tracking known /
 * unknown bits under a partial key guess, and deciding sample validity
 * via the 7 target borrow bits (paper Table 5).
 *
 * Round keys are named as in the paper: rk^18 is the LAST round key,
 * rk^17 the second-to-last. For the 2-round encryption experiment with
 * round keys (Ka, Kb), set rk^17 := Ka and rk^18 := Kb.
 */

// Gamma_out = ([5,6,9,18], [3,4,24,27], [6,29], [0,29])
extern const word GAMMA_MASK[4];

// Bit masks of the guessed (combinations of) round-key bits, 85 bits total.
extern const word M_VY1_18;  // rk0^18          {0:8, 14:28}        24 bits
extern const word M_VX1_18;  // rk1^18          {0, 21, 22, 23}      4 bits
extern const word M_VY2_18;  // rk1^18 ^ rk2^18 {0:8, 15:17, 20:28} 21 bits
extern const word M_VX2_18;  // rk3^18          {24, 25}             2 bits
extern const word M_VY3_18;  // rk3^18 ^ rk4^18 {0:5, 8, 17, 21:28} 16 bits
extern const word M_VY1_17;  // rk5^18 ^ rk0^17 {1:5, 22:28}        12 bits
extern const word M_VY2_17;  // rk1^17 ^ rk2^17 {5, 26:28}           4 bits
extern const word M_VY3_17;  // rk3^17 ^ rk4^17 {27:28}              2 bits

// 85-bit partial-decryption key guess. Only the bits inside the
// corresponding mask of each field are meaningful.
struct KeyBits85 {
    word rk0_18;
    word rk1_18;
    word rk1x2_18;
    word rk3_18;
    word rk3x4_18;
    word rk5x0;
    word rk1x2_17;
    word rk3x4_17;
};

// Extract the correct 85 key bits from the two full round keys
// (Ka = rk^17 = first of the two rounds, Kb = rk^18 = second).
// Bits outside the guessed masks are zeroed: the struct carries exactly
// the 85 guessed bits and nothing more.
KeyBits85 extract_keybits85(const word Ka[6], const word Kb[6]);

// Overwrite the un-guessed bits (outside the masks) of every field with
// the given random words. The 85 guessed bits are preserved. This models
// the attack setting: only the 85 bits are (correctly) guessed, all other
// round-key bits are arbitrary. Note: partial_dec_2r never reads bits
// outside the masks, so this fill cannot change any result.
void kb85_randomize_unguessed(KeyBits85& kb, const word rv[8]);

// Canonical 85-bit indexing of a KeyBits85: fields in declaration order
// (rk0_18, rk1_18, rk1x2_18, rk3_18, rk3x4_18, rk5x0, rk1x2_17, rk3x4_17),
// each field's guessed bits enumerated LSB to MSB. Total 85 bits, idx in [0, 85).
int kb85_get_bit(const KeyBits85& kb, int idx);
void kb85_set_bit(KeyBits85& kb, int idx, int v);

// The 7 target borrow bits of paper Table 5 (fixed order).
enum BorrowTarget {
    B2_17_6 = 0,  // b_2^17[6],  single,        m=1, obs [5]
    B3_18_9,      // b_3^18[9],  single,        m=1, obs [8]
    B3_18_18,     // b_3^18[18], single,        m=1, obs [17]
    B1_17_3,      // b_1^17[3],  single,        m=2, obs [1:2]
    B3_17_29,     // b_3^17[29], two simultaneous,   t=1, m=2
    B2_18_17,     // b_2^18[17], two simultaneous,   t=1, m=2
    B1_17_24,     // b_1^17[24], three simultaneous, s=t=1, m=2
    N_BORROW_TARGETS
};

struct DecResult {
    bool valid;             // all 7 borrows determined for BOTH elements
    bool stat_computable;   // every Gamma_out bit of Delta x^17 computable
    int stat;               // <Gamma_out, Delta x^17> in {0,1}, if computable
    uint8_t borrow_ok;      // bit i: target borrow i determined for both elements
};

// Two-round partial decryption of a ciphertext pair.
// With early_term=true, returns as soon as a target borrow is found
// undeterminable (valid=false); borrow_ok then only contains the targets
// checked so far. Pass early_term=false to always evaluate all 7 targets
// (e.g. when measuring per-borrow ratios).
DecResult partial_dec_2r(const KeyBits85& kb, const block& C, const block& Cp,
                         bool early_term = true);

// <Gamma_out, d> for a fully known difference d.
int gamma_parity(const block& d);
