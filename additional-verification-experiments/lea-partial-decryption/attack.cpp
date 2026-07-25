/**
 * 7-round reduced key-recovery attack on LEA, verifying the paper's
 * two-round partial decryption (Section 6.1, Figure 6) in an attack
 * setting.
 *
 *   7 rounds = 5-round DL distinguisher + 2-round partial decryption D.
 *     Delta_in = (0,[30],0,0) = (0, 40000000, 0, 0)
 *     Gamma_out = ([5,6,9,18],[3,4,24,27],[6,29],[0,29])
 *     distinguisher correlation ~ +2^-2.72
 *
 * The partial decryption covers x^7 -> x^5, isomorphic to the paper's
 * x^19 -> x^17: rk^17 := rk[5], rk^18 := rk[6] of the 7-round schedule.
 *
 * Key guessing: the 85-bit partial key (the combinations of round-key
 * bits listed in the paper, Section 6.1 / Table 6) with 6 bits as the
 * guess space (64 candidates); the other 79 bits are fixed to their
 * correct values.
 *
 * Per candidate k: T(k) = sum over valid pairs of (-1)^<Gamma_out, Delta x^5>,
 * normalized to Q(k) = T(k) / N_v(k). Default setting: a = 30, P_S = 0.9,
 * N_v = 2308, N = 43130000, threshold t = 0.1251 (paper's formula (1)).
 *
 * Valid-sample cap: with N total pairs, at most L = int(N * 2^-14.19) valid
 * samples are USED per candidate (the first L ones); if fewer pairs are
 * valid, all of them are used. This aligns the effective sample size across
 * candidates.
 *
 * Embedded sanity check (--check, default on): for every valid pair of the
 * CORRECT candidate, the statistic is compared with a full two-round
 * decryption under the real round keys (must be 0 mismatches).
 *
 * Build:
 *   g++ -O2 -std=c++17 -o attack attack.cpp partial_dec.cpp tristate.cpp lea.cpp -pthread
 *
 * Example:
 *   ./attack --bits 18,57,8,64,47,23 --trials 100 --seed 100 --threads 8
 */

#include "lea.h"
#include "partial_dec.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr int N_CAND = 64;  // 2^6 candidates

// Default Delta_in = (0,[30],0,0)
const block DELTA_IN_DEFAULT = {0, 0x40000000u, 0, 0};

struct Config {
    uint64_t N = 43130000;   // pairs per attack run (~2^25.4)
    int trials = 1;
    uint64_t seed = 1;
    int threads = 0;           // 0 = hardware concurrency
    int rounds = 7;            // total rounds (distinguisher rounds + 2)
    bool check = true;         // embedded full-decryption sanity check
    bool dump = false;         // print Q of all 64 candidates per trial
    bool random_mode = false;  // candidates are fully random 85-bit guesses
    double threshold = 0.1251; // threshold on Q(k) = T(k)/N_v(k)
    uint64_t cap = 0;          // valid-sample cap L = int(N * 2^-14.19) (set in main)
    std::string bits_str;      // optional fixed guess-bit positions "i0,...,i5"
    bool use_delta = false;    // full Delta_in override via --delta
    word delta_hex[4] = {0, 0, 0, 0};
};

struct CandResult {
    uint64_t Nv = 0;      // valid samples USED (capped at L)
    int64_t T = 0;        // statistic over the used samples
    uint64_t Nv_raw = 0;  // total valid samples seen (correct candidate only)
    uint64_t sanity_mismatch = 0;  // correct candidate only
};

bool parse_args(int argc, char** argv, Config& cfg) {
    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        auto need = [&](const char* name) -> const char* {
            if (a + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                std::exit(1);
            }
            return argv[++a];
        };
        if (arg == "-h" || arg == "--help") {
            std::cout
                << "Usage: " << argv[0] << " [options]\n"
                << "  --N <pairs>      pairs per run (default 43130000)\n"
                << "  --trials <n>     number of attack runs (default 1)\n"
                << "  --seed <s>       RNG seed (default 1)\n"
                << "  --threads <t>    worker threads (default: hardware concurrency)\n"
                << "  --rounds <r>     total rounds (default 7; distinguisher = rounds-2)\n"
                << "  --threshold <f>  Q(k) threshold (default 0.1251)\n"
                << "  --delta <w0,w1,w2,w3>  Delta_in as 4 hex words (default: 0,40000000,0,0)\n"
                << "  --bits <i0,..,i5>  fix the 6 guess-bit positions in [0,85)\n"
                << "  --no-check       disable the full-decryption sanity check\n";
            std::exit(0);
        } else if (arg == "--N") {
            cfg.N = std::strtoull(need("--N"), nullptr, 10);
        } else if (arg == "--trials") {
            cfg.trials = std::atoi(need("--trials"));
        } else if (arg == "--seed") {
            cfg.seed = std::strtoull(need("--seed"), nullptr, 10);
        } else if (arg == "--threads") {
            cfg.threads = std::atoi(need("--threads"));
        } else if (arg == "--rounds") {
            cfg.rounds = std::atoi(need("--rounds"));
        } else if (arg == "--threshold") {
            cfg.threshold = std::atof(need("--threshold"));
        } else if (arg == "--delta") {
            const char* v = need("--delta");
            unsigned long w[4];
            if (std::sscanf(v, "%lx,%lx,%lx,%lx", &w[0], &w[1], &w[2], &w[3]) != 4) {
                std::cerr << "Bad --delta format, expect w0,w1,w2,w3 hex\n";
                return false;
            }
            for (int i = 0; i < 4; ++i) cfg.delta_hex[i] = static_cast<word>(w[i]);
            cfg.use_delta = true;
        } else if (arg == "--bits") {
            cfg.bits_str = need("--bits");
        } else if (arg == "--no-check") {
            cfg.check = false;
        } else if (arg == "--dump") {
            cfg.dump = true;
        } else if (arg == "--random-guesses") {
            cfg.random_mode = true;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return false;
        }
    }
    return true;
}

std::vector<int> parse_bits(const std::string& s) {
    std::vector<int> bits;
    size_t pos = 0;
    while (pos < s.size()) {
        const size_t comma = s.find(',', pos);
        bits.push_back(std::atoi(s.substr(pos, comma - pos).c_str()));
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    return bits;
}

word rand_word(std::mt19937_64& rng) { return static_cast<word>(rng()); }

block rand_block(std::mt19937_64& rng) {
    return {rand_word(rng), rand_word(rng), rand_word(rng), rand_word(rng)};
}

// One worker: regenerate the identical pair stream and process candidates
// [c_lo, c_hi). Results are accumulated into res[c].
void worker(const Config& cfg, uint64_t pair_seed, const word rk[][6],
            const block& delta, const std::vector<KeyBits85>& cand_keys,
            int correct_cand, int c_lo, int c_hi, std::vector<CandResult>& res) {
    std::mt19937_64 rng(pair_seed);
    for (uint64_t i = 0; i < cfg.N; ++i) {
        const block P = rand_block(rng);
        const block Pp = P ^ delta;
        const block C = encrypt(P, rk, cfg.rounds);
        const block Cp = encrypt(Pp, rk, cfg.rounds);
        for (int c = c_lo; c < c_hi; ++c) {
            CandResult& R = res[c];
            // capped candidates are skipped entirely (except that the raw
            // valid count of the correct candidate is still tracked)
            if (R.Nv >= cfg.cap && c != correct_cand) continue;
            const DecResult r = partial_dec_2r(cand_keys[c], C, Cp);
            if (!r.valid) continue;
            if (c == correct_cand) ++R.Nv_raw;
            if (R.Nv >= cfg.cap) continue;  // count raw only, do not use
            ++R.Nv;
            R.T += r.stat ? -1 : 1;
            if (cfg.check && c == correct_cand) {
                // full two-round decryption with the real round keys
                const block X6 = decrypt(C, rk + cfg.rounds - 2, 2);
                const block X6p = decrypt(Cp, rk + cfg.rounds - 2, 2);
                if (r.stat != gamma_parity(X6 ^ X6p)) ++R.sanity_mismatch;
            }
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    Config cfg;
    if (!parse_args(argc, argv, cfg)) return 1;
    if (cfg.threads <= 0) {
        cfg.threads = static_cast<int>(std::thread::hardware_concurrency());
        if (cfg.threads <= 0) cfg.threads = 1;
    }
    cfg.threads = std::min(cfg.threads, N_CAND);
    cfg.cap = static_cast<uint64_t>(static_cast<double>(cfg.N) * std::pow(2.0, -14.19));

    std::vector<int> fixed_bits;
    if (!cfg.bits_str.empty()) {
        fixed_bits = parse_bits(cfg.bits_str);
        if (fixed_bits.size() != 6) {
            std::cerr << "--bits needs exactly 6 positions\n";
            return 1;
        }
    }

    // Delta_in (default: (0,[30],0,0); override with --delta)
    block delta = DELTA_IN_DEFAULT;
    if (cfg.use_delta)
        delta = {cfg.delta_hex[0], cfg.delta_hex[1], cfg.delta_hex[2],
                 cfg.delta_hex[3]};

    std::cout << "# reduced attack on LEA (" << (cfg.rounds - 2)
              << "r distinguisher + 2r partial dec, total " << cfg.rounds
              << " rounds)\n"
              << "# N=" << cfg.N << " trials=" << cfg.trials << " seed=" << cfg.seed
              << " threads=" << cfg.threads << " check=" << (cfg.check ? "on" : "off")
              << "\n# Delta_in=(" << std::hex << std::setfill('0') << std::setw(8)
              << delta.x0 << "," << std::setw(8) << delta.x1 << "," << std::setw(8)
              << delta.x2 << "," << std::setw(8) << delta.x3 << ")\n"
              << std::dec << std::setfill(' ')
              << "# statistic Q(k)=T(k)/N_v(k), threshold t=" << cfg.threshold
              << ", valid-sample cap L=" << cfg.cap << " (= int(N*2^-14.19))\n";

    std::mt19937_64 rng(cfg.seed);
    int n_success = 0;

    for (int t = 0; t < cfg.trials; ++t) {
        const auto t0 = std::chrono::steady_clock::now();

        // master key and round schedule
        word mk[8] = {0};
        for (int i = 0; i < 4; ++i) mk[i] = rand_word(rng);
        word rk[MAX_NR][6];
        expand_key(mk, cfg.rounds, rk, 128);

        // correct 85-bit key and the 64 candidates
        const KeyBits85 correct = extract_keybits85(rk[cfg.rounds - 2], rk[cfg.rounds - 1]);
        int correct_cand = 0;
        std::vector<int> pos;
        std::vector<KeyBits85> cand_keys(N_CAND, correct);
        if (cfg.random_mode) {
            // fully random (wrong) 85-bit guesses; no correct candidate
            correct_cand = -1;
            for (int c = 0; c < N_CAND; ++c) {
                KeyBits85 k;
                k.rk0_18 = rand_word(rng);
                k.rk1_18 = rand_word(rng);
                k.rk1x2_18 = rand_word(rng);
                k.rk3_18 = rand_word(rng);
                k.rk3x4_18 = rand_word(rng);
                k.rk5x0 = rand_word(rng);
                k.rk1x2_17 = rand_word(rng);
                k.rk3x4_17 = rand_word(rng);
                cand_keys[c] = k;
            }
        } else {
            // 6 guess-bit positions
            pos = fixed_bits;
            if (pos.empty()) {
                std::vector<int> all(85);
                for (int i = 0; i < 85; ++i) all[i] = i;
                std::shuffle(all.begin(), all.end(), rng);
                pos.assign(all.begin(), all.begin() + 6);
                std::sort(pos.begin(), pos.end());
            }
            for (int j = 0; j < 6; ++j)
                correct_cand |= kb85_get_bit(correct, pos[j]) << j;
            for (int c = 0; c < N_CAND; ++c)
                for (int j = 0; j < 6; ++j)
                    kb85_set_bit(cand_keys[c], pos[j], (c >> j) & 1);
        }

        const uint64_t pair_seed = rng();

        // run candidates in parallel
        std::vector<CandResult> res(N_CAND);
        {
            std::vector<std::thread> pool;
            std::vector<std::vector<CandResult>> partial(cfg.threads,
                                                         std::vector<CandResult>(N_CAND));
            const int per = (N_CAND + cfg.threads - 1) / cfg.threads;
            for (int th = 0; th < cfg.threads; ++th) {
                const int lo = th * per;
                const int hi = std::min(N_CAND, lo + per);
                if (lo >= hi) break;
                pool.emplace_back(worker, std::cref(cfg), pair_seed, std::cref(rk),
                                  std::cref(delta), std::cref(cand_keys), correct_cand,
                                  lo, hi, std::ref(partial[th]));
            }
            for (auto& th : pool) th.join();
            for (const auto& part : partial)
                for (int c = 0; c < N_CAND; ++c) {
                    res[c].Nv += part[c].Nv;
                    res[c].T += part[c].T;
                    res[c].Nv_raw += part[c].Nv_raw;
                    res[c].sanity_mismatch += part[c].sanity_mismatch;
                }
        }

        // evaluate with the normalized statistic Q(k) = T(k)/N_v(k)
        auto Q = [&](int c) -> double {
            return res[c].Nv ? static_cast<double>(res[c].T) / res[c].Nv
                             : -1.0;
        };
        const auto t1 = std::chrono::steady_clock::now();
        const double secs = std::chrono::duration<double>(t1 - t0).count();

        if (cfg.random_mode) {
            // no correct candidate: report the distribution of Q over the
            // fully random guesses
            double mean = 0, qmax = -1.0, qmin = 1.0;
            int npass = 0;
            for (int c = 0; c < N_CAND; ++c) {
                mean += Q(c);
                qmax = std::max(qmax, Q(c));
                qmin = std::min(qmin, Q(c));
                if (Q(c) > cfg.threshold) ++npass;
            }
            mean /= N_CAND;
            std::cout << "trial " << t << " (random guesses): meanQ="
                      << std::setprecision(4) << mean << " min=" << qmin
                      << " max=" << qmax << " pass(t=" << cfg.threshold
                      << ")=" << npass << "/" << N_CAND << "  (" << std::fixed
                      << std::setprecision(1) << secs << "s)\n";
            if (cfg.dump) {
                std::cout << "  Q values:";
                for (int c = 0; c < N_CAND; ++c) {
                    if (c % 12 == 0) std::cout << "\n   ";
                    std::cout << " " << std::setprecision(3) << Q(c);
                }
                std::cout << "\n";
            }
            continue;
        }

        const double Q_correct = Q(correct_cand);
        double max_wrong_Q = -1.0;
        int wrong_pass = 0;
        int rank = 1;  // rank of correct candidate by Q (1 = best)
        for (int c = 0; c < N_CAND; ++c) {
            if (c == correct_cand) continue;
            max_wrong_Q = std::max(max_wrong_Q, Q(c));
            if (Q(c) > cfg.threshold) ++wrong_pass;
            if (Q(c) >= Q_correct) ++rank;
        }
        const bool success = Q_correct > cfg.threshold;
        n_success += success;

        // wrong candidates passing the threshold, with Hamming distance to correct
        std::string passers;
        for (int c = 0; c < N_CAND; ++c) {
            if (c == correct_cand) continue;
            if (Q(c) > cfg.threshold) {
                char buf[96];
                std::snprintf(buf, sizeof(buf), " (c=%d,hd=%d,Q=%.3f,Nv=%llu)", c,
                              __builtin_popcount(static_cast<unsigned>(c ^ correct_cand)),
                              Q(c), static_cast<unsigned long long>(res[c].Nv));
                passers += buf;
            }
        }

        std::cout << "trial " << t << ": mk=" << std::hex << std::setfill('0')
                  << std::setw(8) << mk[0] << std::setw(8) << mk[1] << std::setw(8)
                  << mk[2] << std::setw(8) << mk[3] << std::dec << std::setfill(' ')
                  << " bits=";
        for (int j = 0; j < 6; ++j) std::cout << pos[j] << (j < 5 ? "," : "");
        std::cout << " correct_cand=" << correct_cand << "\n"
                  << "  valid: raw=" << res[correct_cand].Nv_raw << " (ratio "
                  << std::setprecision(6)
                  << static_cast<double>(res[correct_cand].Nv_raw) / cfg.N
                  << ")  used N_v=" << res[correct_cand].Nv << " (cap L=" << cfg.cap
                  << ")\n"
                  << "  T(correct)=" << res[correct_cand].T
                  << "  Q(correct)=" << std::setprecision(4) << Q_correct
                  << "  t=" << cfg.threshold << "  rank=" << rank << "/64"
                  << "  max Q(wrong)=" << std::setprecision(4) << max_wrong_Q
                  << "  wrong_pass=" << wrong_pass
                  << "  sanity_mismatch=" << res[correct_cand].sanity_mismatch << "\n"
                  << "  wrong passers:" << (wrong_pass ? passers : " none") << "\n"
                  << "  " << (success ? "SUCCESS" : "FAIL") << "  (" << std::fixed
                  << std::setprecision(1) << secs << "s)\n";

        if (cfg.dump) {
            std::cout << "  all candidates [c:Q(hd)] (* = correct):\n";
            for (int c = 0; c < N_CAND; ++c) {
                if (c % 8 == 0) std::cout << "   ";
                char buf[64];
                std::snprintf(buf, sizeof(buf), " %d:%.3f(%d)%s", c, Q(c),
                              __builtin_popcount(static_cast<unsigned>(c ^ correct_cand)),
                              c == correct_cand ? "*" : "");
                std::cout << std::setw(14) << std::left << buf;
                if (c % 8 == 7) std::cout << "\n";
            }
        }
    }

    std::cout << "summary: " << n_success << "/" << cfg.trials << " trials succeeded"
              << "\n";
    return 0;
}
