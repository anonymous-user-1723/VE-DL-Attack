#include "tristate.h"

static word rol32(word x, int r) { return (x << r) | (x >> (32 - r)); }
static word ror32(word x, int r) { return (x >> r) | (x << (32 - r)); }

TWord tw_known(word a, word b) {
    TWord w;
    w.v0 = a;
    w.v1 = b;
    w.k0 = w.k1 = 0xFFFFFFFFu;
    w.dk = 0xFFFFFFFFu;
    w.dv = a ^ b;
    return w;
}

TWord tw_unknown() {
    TWord w;
    w.v0 = w.v1 = 0;
    w.k0 = w.k1 = 0;
    w.dk = w.dv = 0;
    return w;
}

TWord tw_xor_key(const TWord& x, word key_val, word key_mask) {
    TWord r;
    r.v0 = x.v0 ^ key_val;
    r.v1 = x.v1 ^ key_val;
    r.k0 = x.k0 & key_mask;
    r.k1 = x.k1 & key_mask;
    r.dk = x.dk;  // the key cancels in the difference
    r.dv = x.dv;
    return r;
}

TWord tw_rol(const TWord& x, int r) {
    TWord w;
    w.v0 = rol32(x.v0, r);
    w.v1 = rol32(x.v1, r);
    w.k0 = rol32(x.k0, r);
    w.k1 = rol32(x.k1, r);
    w.dk = rol32(x.dk, r);
    w.dv = rol32(x.dv, r);
    return w;
}

TWord tw_ror(const TWord& x, int r) {
    TWord w;
    w.v0 = ror32(x.v0, r);
    w.v1 = ror32(x.v1, r);
    w.k0 = ror32(x.k0, r);
    w.k1 = ror32(x.k1, r);
    w.dk = ror32(x.dk, r);
    w.dv = ror32(x.dv, r);
    return w;
}

// One element of v = z (-) y with three-value borrow propagation,
// word-parallel (parallel-prefix / Kogge-Stone style).
//
//   e = positions where y, z are both known and differ (borrow extinguished:
//       b[j+1] = y[j]);
//   p = positions where both known and equal (borrow propagates).
// Borrow b[j+1] is known iff there is an extinguish position k <= j with all
// positions k+1..j propagating, or bits 0..j all propagate (then b[j+1] =
// b[0] = 0). The borrow value is y at the last extinguish position.
static void sub_one(word zv, word zk, word yv, word yk, word& vv, word& vk, word& bk,
                    word& bv) {
    const word both = yk & zk;
    const word diff = yv ^ zv;
    const word e = diff & both;   // extinguish
    const word p = ~diff & both;  // propagate

    // Parallel prefix over (generate=e, propagate=p), simultaneously
    // spreading the extinguish value B = y at the last extinguish position.
    word G = e, B = yv & e, P = p;
    G |= P & (G << 1);
    B |= (B << 1) & P;
    P &= P << 1;
    G |= P & (G << 2);
    B |= (B << 2) & P;
    P &= P << 2;
    G |= P & (G << 4);
    B |= (B << 4) & P;
    P &= P << 4;
    G |= P & (G << 8);
    B |= (B << 8) & P;
    P &= P << 8;
    G |= P & (G << 16);
    B |= (B << 16) & P;
    // G[j]: borrow into bit j+1 determined by an extinguish at <= j.
    // Base case (no extinguish): borrow into j+1 = b[0] = 0 iff bits 0..j all
    // propagate, i.e. the lowest non-propagate bit is above j:
    const word np = ~p;
    const word AP = (np & (0u - np)) - 1u;  // AND p[0..j] (np==0 -> all ones)
    bk = 1u | ((G | AP) << 1);
    bv = B << 1;  // bv[0] = 0

    vv = yv ^ zv ^ bv;
    vk = both & bk;
}

SubResult tw_sub(const TWord& z, const TWord& y) {
    SubResult r;
    r.v = tw_unknown();
    sub_one(z.v0, z.k0, y.v0, y.k0, r.v.v0, r.v.k0, r.bk0, r.bv0);
    sub_one(z.v1, z.k1, y.v1, y.k1, r.v.v1, r.v.k1, r.bk1, r.bv1);
    // Delta v[j] = Delta y[j] ^ Delta z[j] ^ Delta b[j]; computable iff all
    // three differences are known (borrow difference needs both elements' borrows).
    const word dbk = r.bk0 & r.bk1;
    const word dbv = r.bv0 ^ r.bv1;
    r.v.dk = y.dk & z.dk & dbk;
    r.v.dv = y.dv ^ z.dv ^ dbv;  // meaningful where dk set
    return r;
}
