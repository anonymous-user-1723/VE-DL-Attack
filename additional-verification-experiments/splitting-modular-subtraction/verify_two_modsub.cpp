/**
 * End-to-end verification of simultaneous splitting for TWO modular subtractions.
 * Spec: ../verification-plan.md  (Lemma 4 / Table 3)
 *
 * Default validity test: constructive determination (Lemma 4, attack-style
 * unified CDB path), with a T-trial IT sanity check on each constructive-valid
 * sample. Optional --mode it uses key-perturbation invariance only.
 *
 * Build:
 *   g++ -O2 -std=c++17 -o verify_two_modsub verify_two_modsub.cpp
 *
 * Examples:
 *   ./verify_two_modsub
 *   ./verify_two_modsub --t 1 --m 3
 *   ./verify_two_modsub --n 32 --i 16 --t 1 --m 2 --N 20000 --T 32 --seed 1
 *   ./verify_two_modsub --mode it --t 1 --m 3 --N 20000
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

enum class Mode { Constructive, IT };

struct Config {
    int n = 32;
    int i = 16;
    int t = -1;  // <0 means run default grid
    int m = -1;
    uint64_t N = 100000;
    int T = 32;  // used only in --mode it
    uint64_t seed = 1;
    bool grid = false;
    Mode mode = Mode::Constructive;
    bool print_cdb_breakdown = true;  // IT CDB table; disabled for --grid
};

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n"
        << "\n"
        << "Options (defaults from verification-plan.md):\n"
        << "  --n <int>       word size (default: 32)\n"
        << "  --i <int>       target borrow bit index (default: 16)\n"
        << "  --t <int>       splitting parameter t\n"
        << "  --m <int>       splitting parameter m\n"
        << "  --N <uint>      number of public samples (default: 100000)\n"
        << "  --T <int>       key perturbations per sample for --mode it (default: 32)\n"
        << "  --seed <uint>   RNG seed (default: 1)\n"
        << "  --mode <name>   constructive (default) | it\n"
        << "  --grid          run default (t,m) grid\n"
        << "  -h, --help      show this help\n"
        << "\n"
        << "If neither --t nor --m is set, the default parameter grid is run:\n"
        << "  (t,m) = (1,2),(1,3),(1,4),(1,5),(2,2),(2,3),(2,4)\n"
        << "\n"
        << "Mode notes:\n"
        << "  constructive : Lemma 4 attack-style CDB path + T-trial IT sanity check\n"
        << "  it           : invariance under unknown-key perturbation only\n"
        << "                 (overcounts vs Lemma 4)\n";
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

/** bits x[lo:hi] as integer; empty range -> 0 */
static uint64_t bit_slice(uint64_t x, int lo, int hi) {
    if (lo > hi) return 0;
    const int width = hi - lo + 1;
    if (width >= 64) return x >> lo;
    return (x >> lo) & ((uint64_t{1} << width) - 1);
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
 * Two modular subtractions:
 *   y2 = x3 ^ k2;  x2 = z2 boxminus y2
 *   y1 = x2 ^ k1;  x1 = z1 boxminus y1
 * return b1[i]
 */
static int borrow_bit_2(uint64_t z1, uint64_t z2, uint64_t x3,
                        uint64_t k1, uint64_t k2, int i, uint64_t wmask) {
    const uint64_t y2 = (x3 ^ k2) & wmask;
    const uint64_t x2 = (z2 - y2) & wmask;
    const uint64_t y1 = (x2 ^ k1) & wmask;
    const uint64_t x1 = (z1 - y1) & wmask;
    const uint64_t b1 = x1 ^ y1 ^ z1;
    return get_bit(b1, i);
}

/**
 * Lemma 4 constructive validity (attack-style, unified CDB path).
 *
 * 1) y2 = x3 XOR k2; find first j >= lo2 = i-m-t with CDB^{(y2,z2)}_{lo2,j}
 *    (first bit where y2 and z2 differ). Then b2[j+1] is known => x2[j+1:i-1]
 *    known. If no such j in [lo2:i-1], or j+1 > i-1, not valid.
 * 2) With k1[i-m:i-1] known, y1 is known on [j':i-1] where
 *    j' = max(j+1, i-m). Valid iff (y1 XOR z1)[j':i-1] != 0.
 */
static bool is_valid_constructive_2(uint64_t z1, uint64_t z2, uint64_t x3,
                                    uint64_t k1, uint64_t k2,
                                    int i, int t, int m, uint64_t wmask) {
    const uint64_t y2 = (x3 ^ k2) & wmask;
    const int lo2 = i - m - t;
    const int mid = i - m;
    const int hi = i - 1;

    // First differing bit j of (y2, z2) on [lo2:hi] = CDB_{lo2,j}
    int j = -1;
    for (int b = lo2; b <= hi; ++b) {
        if (get_bit(y2, b) != get_bit(z2, b)) {
            j = b;
            break;
        }
    }
    if (j < 0) {
        return false;  // y2[lo2:hi] == z2[lo2:hi]
    }
    if (j + 1 > hi) {
        return false;  // j == i-1: no x2 bit in [j+1:hi]
    }

    // Full chain to read true y1 on the attacker-known segment
    // (equals what the attack computes on [j':hi]).
    const uint64_t x2 = (z2 - y2) & wmask;
    const uint64_t y1 = (x2 ^ k1) & wmask;

    const int jp = std::max(j + 1, mid);  // j' = max(j+1, i-m)
    return bit_slice(y1 ^ z1, jp, hi) != 0;
}

/**
 * Sanity check for a constructive-valid sample: with known key windows fixed,
 * randomize unknown bits of (k1,k2) T times; b1[i] must always equal the true
 * value. Returns true if all T trials match. On mismatch, prints an error and
 * returns false (caller should abort).
 */
static bool constructive_sanity_it_check(uint64_t z1, uint64_t z2, uint64_t x3,
                                         uint64_t k1, uint64_t k2, int i,
                                         uint64_t known_k1, uint64_t known_k2,
                                         uint64_t wmask, int T,
                                         std::mt19937_64& rng) {
    const int b1_star = borrow_bit_2(z1, z2, x3, k1, k2, i, wmask);
    for (int tau = 0; tau < T; ++tau) {
        const uint64_t k1p = perturb_unknown(k1, known_k1, wmask, rng);
        const uint64_t k2p = perturb_unknown(k2, known_k2, wmask, rng);
        const int b1 = borrow_bit_2(z1, z2, x3, k1p, k2p, i, wmask);
        if (b1 != b1_star) {
            std::cerr << "ERROR: constructive-valid sample failed IT sanity check\n"
                      << "  i=" << i << " tau=" << tau << "/" << T << "\n"
                      << "  b1_star=" << b1_star << " b1_perturbed=" << b1 << "\n"
                      << "  z1=0x" << std::hex << z1 << " z2=0x" << z2
                      << " x3=0x" << x3 << "\n"
                      << "  k1=0x" << k1 << " k2=0x" << k2 << "\n"
                      << "  k1p=0x" << k1p << " k2p=0x" << k2p << std::dec
                      << "\n";
            return false;
        }
    }
    return true;
}

static bool is_valid_it_2(uint64_t z1, uint64_t z2, uint64_t x3,
                          uint64_t k1, uint64_t k2,
                          int i, int /*t*/, int /*m*/,
                          uint64_t known_k1, uint64_t known_k2,
                          uint64_t wmask, int T, std::mt19937_64& rng) {
    const int b1_star = borrow_bit_2(z1, z2, x3, k1, k2, i, wmask);
    for (int tau = 0; tau < T; ++tau) {
        const uint64_t k1p = perturb_unknown(k1, known_k1, wmask, rng);
        const uint64_t k2p = perturb_unknown(k2, known_k2, wmask, rng);
        if (borrow_bit_2(z1, z2, x3, k1p, k2p, i, wmask) != b1_star) {
            return false;
        }
    }
    return true;
}

static double theory_p2(int t, int m) {
    return 1.0 - (1.0 + static_cast<double>(m) * std::ldexp(1.0, -t)) *
                     std::ldexp(1.0, -m);
}

static double theory_p2_obo(int t, int m) {
    return (1.0 - std::ldexp(1.0, -t)) * (1.0 - std::ldexp(1.0, -m));
}

/** Classify (y2,z2) on [i1:i2]:
 *  return j if CDB_{i1,j}^{(y2,z2)} for some j in [i1,i2];
 *  return -1 if y2[i1:i2] == z2[i1:i2] (entire window equal).
 */
static int classify_cdb_or_equal(uint64_t y2, uint64_t z2, int i1, int i2) {
    for (int j = i1; j <= i2; ++j) {
        const bool prefix_zero =
            (j == i1) ? true : (bit_slice(y2 ^ z2, i1, j - 1) == 0);
        if (prefix_zero && get_bit(y2, j) != get_bit(z2, j)) {
            return j;
        }
    }
    return -1;
}

struct RunResult {
    int n = 0;
    int i = 0;
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

    // IT-mode stratification by CDB^{(y2,z2)}_{i1,j} / window-equal
    // bins[0 .. width-1] <-> CDB at j = i1 + offset
    // bins[width]        <-> y2[i1:i2] == z2[i1:i2]
    bool has_it_breakdown = false;
    int cdb_i1 = 0;  // i-m-t
    int cdb_i2 = 0;  // i-1
    int cdb_width = 0;
    std::vector<uint64_t> bin_total;
    std::vector<uint64_t> bin_valid;
};

static RunResult run_one(const Config& base, int t, int m, uint64_t seed_offset) {
    RunResult r;
    r.n = base.n;
    r.t = t;
    r.m = m;
    r.N = base.N;
    r.T = base.T;
    r.seed = base.seed + seed_offset;
    r.p_theory = theory_p2(t, m);
    r.p_obo = theory_p2_obo(t, m);

    if (base.n < 2 || base.n > 64) {
        r.skip_reason = "n must be in [2,64]";
        return r;
    }
    if (t < 1 || m < 1) {
        r.skip_reason = "t,m must be >= 1";
        return r;
    }
    if (base.T < 1) {
        r.skip_reason = "T must be >= 1 (IT mode / constructive sanity check)";
        return r;
    }
    if (base.N < 1) {
        r.skip_reason = "N must be >= 1";
        return r;
    }

    int i = base.i;
    if (i <= t + m) {
        i = t + m + 1;
    }
    if (i >= base.n) {
        r.skip_reason = "i >= n after adjustment (need i > t+m and i < n)";
        return r;
    }
    r.i = i;

    const uint64_t wmask = word_mask(base.n);
    const uint64_t known_k1 = make_known_mask(i - m, i - 1, base.n);
    const uint64_t known_k2 = make_known_mask(i - m - t, i - 1, base.n);

    int cdb_i1 = 0;
    int cdb_i2 = 0;
    int cdb_width = 0;
    if (base.mode == Mode::IT && base.print_cdb_breakdown) {
        cdb_i1 = i - m - t;
        cdb_i2 = i - 1;
        cdb_width = m + t;
        const size_t n_bins = size_t(cdb_width) + 1;
        r.has_it_breakdown = true;
        r.cdb_i1 = cdb_i1;
        r.cdb_i2 = cdb_i2;
        r.cdb_width = cdb_width;
        r.bin_total.assign(n_bins, 0);
        r.bin_valid.assign(n_bins, 0);
    }

    std::mt19937_64 rng(r.seed);

    uint64_t count_valid = 0;
    for (uint64_t s = 0; s < base.N; ++s) {
        // Fixed true keys + random public words (Lemma 4 sampling model).
        // Keys are re-drawn each sample; by symmetry this matches fixing keys.
        const uint64_t k1 = random_word(rng, wmask);
        const uint64_t k2 = random_word(rng, wmask);
        const uint64_t z1 = random_word(rng, wmask);
        const uint64_t z2 = random_word(rng, wmask);
        const uint64_t x3 = random_word(rng, wmask);

        bool is_valid = false;
        if (base.mode == Mode::Constructive) {
            if (is_valid_constructive_2(z1, z2, x3, k1, k2, i, t, m, wmask)) {
                // Constructive-valid => b1[i] must be invariant under unknown-key
                // randomization; verify with T trials before counting.
                if (!constructive_sanity_it_check(z1, z2, x3, k1, k2, i,
                                                  known_k1, known_k2, wmask,
                                                  base.T, rng)) {
                    r.skip_reason =
                        "constructive-valid sample failed IT sanity check "
                        "(see stderr)";
                    return r;
                }
                is_valid = true;
            }
        } else {
            is_valid = is_valid_it_2(z1, z2, x3, k1, k2, i, t, m,
                                     known_k1, known_k2, wmask, base.T, rng);
            if (r.has_it_breakdown) {
                // Stratify by CDB^{(y2,z2)}_{i-m-t,j} / window equal
                const uint64_t y2 = (x3 ^ k2) & wmask;
                const int cdb_j = classify_cdb_or_equal(y2, z2, cdb_i1, cdb_i2);
                const size_t bin =
                    (cdb_j < 0) ? size_t(cdb_width)  // equal case
                                : size_t(cdb_j - cdb_i1);
                ++r.bin_total[bin];
                if (is_valid) ++r.bin_valid[bin];
            }
        }
        if (is_valid) ++count_valid;
    }

    r.count_valid = count_valid;
    r.p_hat = static_cast<double>(count_valid) / static_cast<double>(base.N);
    r.abs_err = std::fabs(r.p_hat - r.p_theory);
    const double p = r.p_theory;
    r.threshold = 5.0 * std::sqrt(p * (1.0 - p) / static_cast<double>(base.N));
    // constructive: |p_hat - p_theory| within statistical threshold
    // it: only require p_hat >= p_theory (IT rate is expected to be >= lemma rate)
    if (base.mode == Mode::IT) {
        r.pass = (r.p_hat >= r.p_theory);
    } else {
        r.pass = (r.abs_err <= r.threshold);
    }
    return r;
}

static void print_header() {
    std::cout << std::left
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

static void print_it_breakdown(const RunResult& r) {
    if (!r.has_it_breakdown) return;

    std::cout << "# IT breakdown by CDB^{(y2,z2)}_{" << r.cdb_i1 << ",j} "
              << "for j in {" << r.cdb_i1 << ",...," << r.cdb_i2 << "} "
              << "and case y2[" << r.cdb_i1 << ":" << r.cdb_i2
              << "]=z2[" << r.cdb_i1 << ":" << r.cdb_i2 << "]\n";
    std::cout << std::left
              << std::setw(28) << "case"
              << std::setw(8) << "j"
              << std::setw(12) << "count"
              << std::setw(12) << "valid"
              << std::setw(14) << "frac_in_bin"
              << std::setw(14) << "p_valid|bin"
              << std::setw(14) << "contrib"
              << "\n";

    uint64_t sum_count = 0;
    uint64_t sum_valid = 0;
    for (size_t b = 0; b < r.bin_total.size(); ++b) {
        const uint64_t c = r.bin_total[b];
        const uint64_t vv = r.bin_valid[b];
        sum_count += c;
        sum_valid += vv;
        const double frac_in_bin =
            (r.N > 0) ? static_cast<double>(c) / static_cast<double>(r.N) : 0.0;
        const double p_cond =
            (c > 0) ? static_cast<double>(vv) / static_cast<double>(c) : 0.0;
        const double contrib =
            (r.N > 0) ? static_cast<double>(vv) / static_cast<double>(r.N) : 0.0;

        std::string case_name;
        int j_print = -1;
        if (static_cast<int>(b) < r.cdb_width) {
            j_print = r.cdb_i1 + static_cast<int>(b);
            case_name = "CDB_{" + std::to_string(r.cdb_i1) + "," +
                        std::to_string(j_print) + "}";
        } else {
            case_name = "y2==z2 on window";
        }

        std::cout << std::left
                  << std::setw(28) << case_name
                  << std::setw(8) << (j_print < 0 ? "-" : std::to_string(j_print))
                  << std::setw(12) << c
                  << std::setw(12) << vv
                  << std::setw(14) << std::setprecision(6) << std::fixed << frac_in_bin
                  << std::setw(14) << p_cond
                  << std::setw(14) << contrib
                  << "\n";
    }

    std::cout << "# check: sum_count=" << sum_count << " sum_valid=" << sum_valid
              << " p_hat=" << (r.N ? double(sum_valid) / double(r.N) : 0.0) << "\n";
}

static void print_row(const RunResult& r) {
    if (!r.skip_reason.empty()) {
        std::cout << std::left
                  << std::setw(4) << r.t
                  << std::setw(4) << r.m
                  << " SKIP: " << r.skip_reason << "\n";
        return;
    }
    std::cout << std::left
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
    print_it_breakdown(r);
}

int main(int argc, char** argv) {
    Config cfg;
    if (!parse_args(argc, argv, cfg)) return 1;

    std::vector<std::pair<int, int>> jobs;
    const bool want_grid = cfg.grid || (cfg.t < 0 && cfg.m < 0);
    if (want_grid) {
        cfg.print_cdb_breakdown = false;  // grid: summary table only
        const int grid_tm[][2] = {
            {1, 2}, {1, 3}, {1, 4}, {1, 5},
            {2, 2}, {2, 3}, {2, 4},
        };
        for (const auto& tm : grid_tm) jobs.emplace_back(tm[0], tm[1]);
    } else if (cfg.t >= 0 && cfg.m >= 0) {
        jobs.emplace_back(cfg.t, cfg.m);
    } else {
        std::cerr << "Both --t and --m must be set for a single run "
                     "(or omit both / use --grid to run the default grid).\n";
        return 1;
    }

    std::cout << "# Two modular subtractions (simultaneous splitting)\n"
              << "# mode=" << (cfg.mode == Mode::Constructive ? "constructive" : "it")
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
        const int t = jobs[idx].first;
        const int m = jobs[idx].second;
        RunResult r = run_one(cfg, t, m, static_cast<uint64_t>(idx) * 10007ULL);
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
