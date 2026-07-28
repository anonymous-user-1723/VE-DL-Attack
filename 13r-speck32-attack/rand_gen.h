#pragma once
#include <stdint.h>
#include <random>

class RandomGenerator {
    private:
        std::mt19937 engine;
        std::uniform_int_distribution<uint16_t> dist_uint16;
    public:
        void set_rand_seed(uint32_t);
        uint16_t random_16bits();
        RandomGenerator(uint32_t);
};

extern RandomGenerator global_random_generator;