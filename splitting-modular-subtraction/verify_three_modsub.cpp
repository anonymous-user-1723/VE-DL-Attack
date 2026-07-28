/**
 * Verification of simultaneous splitting for THREE modular subtractions.
 * Spec: ../verification-plan.md (Lemma 5 / Table 3)
 *
 * Default: constructive (attack-style chained CDB) + T-trial IT sanity check.
 * Optional --mode it: key-perturbation invariance only.
 *
 * Build:
 *   g++ -O2 -std=c++17 -o verify_three_modsub verify_three_modsub.cpp
 *
 * Examples:
 *   ./verify_three_modsub
 *   ./verify_three_modsub --s 1 --t 1 --m 3
 *   ./verify_three_modsub --mode it --s 1 --t 1 --m 2 --N 20000 --T 32
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

enum class Mode { Constructive, IT };

struct Config {
    int n = 32;
    int i = 16;
    int s = -1;  // <0 means run default grid
    int t = -1;
    int m = -1;
    uint64_t N = 100000;
    int T = 32;
    uint64_t seed = 1;
    bool grid = false;
    Mode mode = Mode::Constructive;
};

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n"
        << "\n"
        << "Options (defaults from verification-plan.md):\n"
        << "  --n <int>       word size (default: 32)\n"
        << "  --i <int>       target borrow bit index (default: 16)\n"
        << "  --s <int>       splitting parameter s\n"
        << "  --t <int>       splitting parameter t\n"
        << "  --m <int>       splitting parameter m\n"
        << "  --N <uint>      number of public samples (default: 100000)\n"
        << "  --T <int>       key perturbations per sample (default: 32)\n"
        << "  --seed <uint>   RNG seed (default: 1)\n"
        << "  --mode <name>   constructive (default) | it\n"
        << "  --grid          run default (s,t,m) grid\n"
        << "  -h, --help      show this help\n"
        << "\n"
        << "If none of --s/--t/--m is set, the default parameter grid is run:\n"
        << "  (s,t,m) = (1,1,2),(1,1,3),(1,1,4),(1,1,5),(1,2,2),(1,2,3),(1,2,4)\n"
        << "\n"
        << "Mode notes:\n"
        << "  constructive : Lemma 5 attack-style chained CDB + T-trial IT sanity\n"
        << "  it           : unknown-key perturbation invariance only\n"
        << "Pass: constructive uses |p_hat-p| threshold; it requires p_hat >= p.\n";
}

static bool parse_args(int argc, char** argv, Config& cfg) {
    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        auto need = [&](const char* name) -> const char* {
            if (a + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                return nullptr;
            }
            return argv[++a];
        };

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            std::exit(0);
        } else if (arg == "--n") {
            const char* v = need("--n");
            if (!v) return false;
            cfg.n = std::atoi(v);
        } else if (arg == "--i") {
            const char* v = need("--i");
            if (!v) return false;
            cfg.i = std::atoi(v);
        } else if (arg == "--s") {
            const char* v = need("--s");
            if (!v) return false;
            cfg.s = std::atoi(v);
        } else if (arg == "--t") {
            const char* v = need("--t");
            if (!v) return false;
            cfg.t = std::atoi(v);
        } else if (arg == "--m") {
            const char* v = need("--m");
            if (!v) return false;
            cfg.m = std::atoi(v);
        } else if (arg == "--N") {
            const char* v = need("--N");
            if (!v) return false;
            cfg.N = static_cast<uint64_t>(std::strtoull(v, nullptr, 10));
        } else if (arg == "--T") {
            const char* v = need("--T");
            if (!v) return false;
            cfg.T = std::atoi(v);
        } else if (arg == "--seed") {
            const char* v = need("--seed");
            if (!v) return false;
            cfg.seed = static_cast<uint64_t>(std::strtoull(v, nullptr, 10));
        } else if (arg == "--grid") {
            cfg.grid = true;
        } else if (arg == "--mode") {
            const char* v = need("--mode");
            if (!v) return false;
            std::string m = v;
            if (m == "constructive" || m == "lemma" || m == "c") {
                cfg.mode = Mode::Constructive;
            } else if (m == "it" || m == "invariance" || m == "perturb") {
                cfg.mode = Mode::IT;
            } else {
                std::cerr << "Unknown --mode: " << m << "\n";
                return false;
            }
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return false;
        }
    }
    return true;
}

static uint64_t word_mask(int n) {
    if (n <= 0 || n > 64) return 0;
    if (n == 64) return ~uint64_t{0};
    return (uint64_t{1} << n) - 1;
}

static uint64_t make_known_mask(int lo, int hi, int n) {
    if (lo > hi) return 0;
    uint64_t m = 0;
    for (int b = lo; b <= hi; ++b) {
        if (b >= 0 && b < n) m |= (uint64_t{1} << b);
    }
    return m;
}

static int get_bit(uint64_t x, int bit) {
    return static_cast<int>((x >> bit) & 1ULL);
}

static uint64_t bit_slice(uint64_t x, int lo, int hi) {
    if (lo > hi) return 0;
    const int width = hi - lo + 1;
    if (width >= 64) return x >> lo;
    return (x >> lo) & ((uint64_t{1} << width) - 1);
}

/** First bit in [lo:hi] where y and z differ; -1 if none or lo>hi. */
static int find_first_diff(uint64_t y, uint64_t z, int lo, int hi) {
    if (lo > hi) return -1;
    for (int b = lo; b <= hi; ++b) {
        if (get_bit(y, b) != get_bit(z, b)) return b;
    }
    return -1;
}

static uint64_t random_word(std::mt19937_64& rng, uint64_t wmask) {
    return rng() & wmask;
}

static uint64_t perturb_unknown(uint64_t k, uint64_t known_mask, uint64_t wmask,
                                std::mt19937_64& rng) {
    const uint64_t rnd = random_word(rng, wmask);
    const uint64_t perturb_mask = (~known_mask) & wmask;
    return (k & known_mask) | (rnd & perturb_mask);
}

/**
 * Three modular subtractions (Fig. 4):
 *   y3 = x4 ^ k3;  x3 = z3 boxminus y3
 *   y2 = x3 ^ k2;  x2 = z2 boxminus y2
 *   y1 = x2 ^ k1;  x1 = z1 boxminus y1
 * return b1[i]
 */
static int borrow_bit_3(uint64_t z1, uint64_t z2, uint64_t z3, uint64_t x4,
                        uint64_t k1, uint64_t k2, uint64_t k3,
                        int i, uint64_t wmask) {
    const uint64_t y3 = (x4 ^ k3) & wmask;
    const uint64_t x3 = (z3 - y3) & wmask;
    const uint64_t y2 = (x3 ^ k2) & wmask;
    const uint64_t x2 = (z2 - y2) & wmask;
    const uint64_t y1 = (x2 ^ k1) & wmask;
    const uint64_t x1 = (z1 - y1) & wmask;
    const uint64_t b1 = x1 ^ y1 ^ z1;
    return get_bit(b1, i);
}

/**
 * Lemma 5 constructive validity (attack-style chained CDB).
 *
 * 1) y3=x4XORk3; j3 = first diff of (y3,z3) from lo3=i-m-t-s.
 *    If no j3 or j3+1 > i-1 => invalid (x3[j3+1:i-1] empty).
 * 2) i2=max(j3+1,i-m-t); y2 known on [i2:i-1]; j2 = first diff of (y2,z2)
 *    from i2. If no j2 or j2+1 > i-1 => invalid.
 * 3) i1=max(j2+1,i-m); valid iff (y1XORz1)[i1:i-1] != 0.
 */
static bool is_valid_constructive_3(uint64_t z1, uint64_t z2, uint64_t z3,
                                    uint64_t x4, uint64_t k1, uint64_t k2,
                                    uint64_t k3, int i, int s, int t, int m,
                                    uint64_t wmask) {
    const int hi = i - 1;
    const int lo3 = i - m - t - s;
    const int k2_lo = i - m - t;
    const int k1_lo = i - m;

    // --- step 1: CDB on third subtraction ---
    const uint64_t y3 = (x4 ^ k3) & wmask;
    const int j3 = find_first_diff(y3, z3, lo3, hi);
    if (j3 < 0 || j3 + 1 > hi) return false;

    const uint64_t x3 = (z3 - y3) & wmask;
    const uint64_t y2 = (x3 ^ k2) & wmask;

    // --- step 2: CDB on second subtraction from i2 ---
    const int i2 = std::max(j3 + 1, k2_lo);
    if (i2 > hi) return false;

    const int j2 = find_first_diff(y2, z2, i2, hi);
    if (j2 < 0 || j2 + 1 > hi) return false;

    const uint64_t x2 = (z2 - y2) & wmask;
    const uint64_t y1 = (x2 ^ k1) & wmask;

    // --- step 3: first-subtraction window for b1[i] ---
    const int i1 = std::max(j2 + 1, k1_lo);
    if (i1 > hi) return false;

    return bit_slice(y1 ^ z1, i1, hi) != 0;
}

/**
 * Sanity check for a constructive-valid sample: fix known key windows, randomize
 * unknown bits of (k1,k2,k3) T times; b1[i] must match the true value.
 */
static bool constructive_sanity_it_check(
    uint64_t z1, uint64_t z2, uint64_t z3, uint64_t x4, uint64_t k1, uint64_t k2,
    uint64_t k3, int i, uint64_t known_k1, uint64_t known_k2, uint64_t known_k3,
    uint64_t wmask, int T, std::mt19937_64& rng) {
    const int b1_star = borrow_bit_3(z1, z2, z3, x4, k1, k2, k3, i, wmask);
    for (int tau = 0; tau < T; ++tau) {
        const uint64_t k1p = perturb_unknown(k1, known_k1, wmask, rng);
        const uint64_t k2p = perturb_unknown(k2, known_k2, wmask, rng);
        const uint64_t k3p = perturb_unknown(k3, known_k3, wmask, rng);
        const int b1 = borrow_bit_3(z1, z2, z3, x4, k1p, k2p, k3p, i, wmask);
        if (b1 != b1_star) {
            std::cerr << "ERROR: constructive-valid sample failed IT sanity check\n"
                      << "  i=" << i << " tau=" << tau << "/" << T << "\n"
                      << "  b1_star=" << b1_star << " b1_perturbed=" << b1 << "\n"
                      << "  z1=0x" << std::hex << z1 << " z2=0x" << z2
                      << " z3=0x" << z3 << " x4=0x" << x4 << "\n"
                      << "  k1=0x" << k1 << " k2=0x" << k2 << " k3=0x" << k3
                      << "\n"
                      << "  k1p=0x" << k1p << " k2p=0x" << k2p << " k3p=0x"
                      << k3p << std::dec << "\n";
            return false;
        }
    }
    return true;
}

static bool is_valid_it_3(uint64_t z1, uint64_t z2, uint64_t z3, uint64_t x4,
                          uint64_t k1, uint64_t k2, uint64_t k3, int i,
                          uint64_t known_k1, uint64_t known_k2, uint64_t known_k3,
                          uint64_t wmask, int T, std::mt19937_64& rng) {
    const int b1_star =
        borrow_bit_3(z1, z2, z3, x4, k1, k2, k3, i, wmask);
    for (int tau = 0; tau < T; ++tau) {
        const uint64_t k1p = perturb_unknown(k1, known_k1, wmask, rng);
        const uint64_t k2p = perturb_unknown(k2, known_k2, wmask, rng);
        const uint64_t k3p = perturb_unknown(k3, known_k3, wmask, rng);
        if (borrow_bit_3(z1, z2, z3, x4, k1p, k2p, k3p, i, wmask) != b1_star) {
            return false;
        }
    }
    return true;
}

/** Lemma 5 simultaneous splitting proportion */
static double theory_p3(int s, int t, int m) {
    const double mt = static_cast<double>(m) * std::ldexp(1.0, -t);
    const double inner =
        1.0 + static_cast<double>(m + 2 * t - 1) * std::ldexp(1.0, -(s + 1));
    return 1.0 - (1.0 + mt * inner) * std::ldexp(1.0, -m);
}

static double theory_p3_obo(int s, int t, int m) {
    return (1.0 - std::ldexp(1.0, -s)) * (1.0 - std::ldexp(1.0, -t)) *
           (1.0 - std::ldexp(1.0, -m));
}

struct RunResult {
    int n = 0;
    int i = 0;
    int s = 0;
    int t = 0;
    int m = 0;
    uint64_t N = 0;
    int T = 0;
    uint64_t seed = 0;
    uint64_t count_valid = 0;
    double p_hat = 0;
    double p_theory = 0;
    double p_obo = 0;
    double abs_err = 0;
    double threshold = 0;
    bool pass = false;
    std::string skip_reason;
};

static RunResult run_one(const Config& base, int s, int t, int m,
                         uint64_t seed_offset) {
    RunResult r;
    r.n = base.n;
    r.s = s;
    r.t = t;
    r.m = m;
    r.N = base.N;
    r.T = base.T;
    r.seed = base.seed + seed_offset;
    r.p_theory = theory_p3(s, t, m);
    r.p_obo = theory_p3_obo(s, t, m);

    if (base.n < 2 || base.n > 64) {
        r.skip_reason = "n must be in [2,64]";
        return r;
    }
    if (s < 1 || t < 1 || m < 1) {
        r.skip_reason = "s,t,m must be >= 1";
        return r;
    }
    if (base.T < 1) {
        r.skip_reason = "T must be >= 1";
        return r;
    }
    if (base.N < 1) {
        r.skip_reason = "N must be >= 1";
        return r;
    }

    int i = base.i;
    if (i <= s + t + m) {
        i = s + t + m + 1;
    }
    if (i >= base.n) {
        r.skip_reason = "i >= n after adjustment (need i > s+t+m and i < n)";
        return r;
    }
    r.i = i;

    const uint64_t wmask = word_mask(base.n);
    const uint64_t known_k1 = make_known_mask(i - m, i - 1, base.n);
    const uint64_t known_k2 = make_known_mask(i - m - t, i - 1, base.n);
    const uint64_t known_k3 = make_known_mask(i - m - t - s, i - 1, base.n);

    std::mt19937_64 rng(r.seed);

    uint64_t count_valid = 0;
    for (uint64_t sample = 0; sample < base.N; ++sample) {
        const uint64_t k1 = random_word(rng, wmask);
        const uint64_t k2 = random_word(rng, wmask);
        const uint64_t k3 = random_word(rng, wmask);
        const uint64_t z1 = random_word(rng, wmask);
        const uint64_t z2 = random_word(rng, wmask);
        const uint64_t z3 = random_word(rng, wmask);
        const uint64_t x4 = random_word(rng, wmask);

        bool is_valid = false;
        if (base.mode == Mode::Constructive) {
            if (is_valid_constructive_3(z1, z2, z3, x4, k1, k2, k3, i, s, t, m,
                                        wmask)) {
                if (!constructive_sanity_it_check(
                        z1, z2, z3, x4, k1, k2, k3, i, known_k1, known_k2,
                        known_k3, wmask, base.T, rng)) {
                    r.skip_reason =
                        "constructive-valid sample failed IT sanity check "
                        "(see stderr)";
                    return r;
                }
                is_valid = true;
            }
        } else {
            is_valid = is_valid_it_3(z1, z2, z3, x4, k1, k2, k3, i, known_k1,
                                     known_k2, known_k3, wmask, base.T, rng);
        }
        if (is_valid) ++count_valid;
    }

    r.count_valid = count_valid;
    r.p_hat = static_cast<double>(count_valid) / static_cast<double>(base.N);
    r.abs_err = std::fabs(r.p_hat - r.p_theory);
    const double p = r.p_theory;
    r.threshold = 5.0 * std::sqrt(p * (1.0 - p) / static_cast<double>(base.N));
    if (base.mode == Mode::IT) {
        r.pass = (r.p_hat >= r.p_theory);
    } else {
        r.pass = (r.abs_err <= r.threshold);
    }
    return r;
}

static void print_header() {
    std::cout << std::left
              << std::setw(4) << "s"
              << std::setw(4) << "t"
              << std::setw(4) << "m"
              << std::setw(4) << "i"
              << std::setw(10) << "N"
              << std::setw(6) << "T"
              << std::setw(12) << "p_hat"
              << std::setw(12) << "p_theory"
              << std::setw(12) << "p_obo"
              << std::setw(12) << "abs_err"
              << std::setw(12) << "threshold"
              << std::setw(6) << "pass"
              << "\n";
}

static void print_row(const RunResult& r) {
    if (!r.skip_reason.empty()) {
        std::cout << std::left
                  << std::setw(4) << r.s
                  << std::setw(4) << r.t
                  << std::setw(4) << r.m
                  << " SKIP: " << r.skip_reason << "\n";
        return;
    }
    std::cout << std::left
              << std::setw(4) << r.s
              << std::setw(4) << r.t
              << std::setw(4) << r.m
              << std::setw(4) << r.i
              << std::setw(10) << r.N
              << std::setw(6) << r.T
              << std::setw(12) << std::setprecision(8) << std::fixed << r.p_hat
              << std::setw(12) << r.p_theory
              << std::setw(12) << r.p_obo
              << std::setw(12) << r.abs_err
              << std::setw(12) << r.threshold
              << std::setw(6) << (r.pass ? "Y" : "N")
              << "\n";
}

int main(int argc, char** argv) {
    Config cfg;
    if (!parse_args(argc, argv, cfg)) return 1;

    std::vector<std::tuple<int, int, int>> jobs;
    const bool any_stm = (cfg.s >= 0 || cfg.t >= 0 || cfg.m >= 0);
    const bool want_grid = cfg.grid || !any_stm;
    if (want_grid) {
        const int grid_stm[][3] = {
            {1, 1, 2}, {1, 1, 3}, {1, 1, 4}, {1, 1, 5},
            {1, 2, 2}, {1, 2, 3}, {1, 2, 4},
        };
        for (const auto& stm : grid_stm) {
            jobs.emplace_back(stm[0], stm[1], stm[2]);
        }
    } else if (cfg.s >= 0 && cfg.t >= 0 && cfg.m >= 0) {
        jobs.emplace_back(cfg.s, cfg.t, cfg.m);
    } else {
        std::cerr << "All of --s, --t, and --m must be set for a single run "
                     "(or omit all / use --grid to run the default grid).\n";
        return 1;
    }

    std::cout << "# Three modular subtractions (simultaneous splitting)\n"
              << "# mode="
              << (cfg.mode == Mode::Constructive ? "constructive" : "it")
              << " n=" << cfg.n
              << " i_requested=" << cfg.i
              << " N=" << cfg.N
              << " T=" << cfg.T
              << " seed=" << cfg.seed
              << "\n";
    print_header();

    int fails = 0;
    int skipped = 0;
    for (size_t idx = 0; idx < jobs.size(); ++idx) {
        const int s = std::get<0>(jobs[idx]);
        const int t = std::get<1>(jobs[idx]);
        const int m = std::get<2>(jobs[idx]);
        RunResult r = run_one(cfg, s, t, m, static_cast<uint64_t>(idx) * 10007ULL);
        print_row(r);
        if (!r.skip_reason.empty()) {
            ++skipped;
        } else if (!r.pass) {
            ++fails;
        }
    }

    std::cout << "# done: jobs=" << jobs.size()
              << " fail=" << fails
              << " skip=" << skipped
              << "\n";
    return (fails > 0) ? 2 : 0;
}
