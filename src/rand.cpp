#include <random>
#include "rand.h"

// Set up our generator once, and keep calling from it.

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis(0.0, 1.0);

double get_rand() {
	return dis(gen);
}

// This works better than rand() / RAND_MAX