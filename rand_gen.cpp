#include "rand_gen.h"
#include <string.h>
#include <stdio.h>

RandomGenerator::RandomGenerator(uint32_t seed) {
	engine.seed(seed);
	dist_uint16 = std::uniform_int_distribution<uint16_t>(0, 0xFFFFU);
}

void RandomGenerator::set_rand_seed(uint32_t seed) {
	engine.seed(seed);
}

uint16_t RandomGenerator::random_16bits() {
	return dist_uint16(engine);
}

RandomGenerator global_random_generator(2333);