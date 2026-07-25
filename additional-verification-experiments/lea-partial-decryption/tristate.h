#pragma once
#include "lea.h"

/**
 * Three-value (0/1/unknown) propagation of a word PAIR through the LEA
 * two-round partial decryption (paper Section 6.1, Figure 6).
 *
 * For each element e in {0, 1} of the pair we track which bit values are
 * known under the current partial key guess. For the pair we additionally
 * track which DIFFERENCE bits are computable: a difference bit may be
 * computable even when neither element's value is known (blue bits,
 * D_x \ V_x), via Delta v[i] = Delta y[i] ^ Delta z[i] ^ Delta b[i] with a
 * recovered borrow difference.
 *
 * Invariant: dk is a superset of (k0 & k1), and where k0 & k1 has a bit,
 * dv equals v0 ^ v1 there.
 */
struct TWord {
    word v0, v1;  // bit values of the two elements (meaningful where k0/k1 set)
    word k0, k1;  // per-element known masks
    word dk, dv;  // dk: difference-known mask; dv: difference bits (valid where dk set)
};

TWord tw_known(word a, word b);  // fully known pair (e.g. values from ciphertext)
TWord tw_unknown();              // nothing known

// XOR with partially-known key bits: value known only where both the input
// bit and the key bit are known; the difference is unaffected by the key.
TWord tw_xor_key(const TWord& x, word key_val, word key_mask);

TWord tw_rol(const TWord& x, int r);
TWord tw_ror(const TWord& x, int r);

struct SubResult {
    TWord v;        // v = z (-) y, three-value
    word bk0, bv0;  // borrow known-mask / borrow bits of element 0 (bit j = b[j])
    word bk1, bv1;  // element 1
};

// Modular subtraction z (-) y with three-value borrow propagation:
//   b[0] = 0; if y[j], z[j] both known and differ -> b[j+1] = y[j] (extinguish);
//   if both known and equal -> propagate b[j]; otherwise b[j+1] unknown
//   (until the next detectable extinguish bit).
SubResult tw_sub(const TWord& z, const TWord& y);
