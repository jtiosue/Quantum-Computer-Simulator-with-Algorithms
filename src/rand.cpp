#include <ctime>
#include <cstdlib>
#include "rand.h"

void set_srand() {
	srand(time(NULL));
	rand();
	// The first one is always basically the same for some reason.
	// So get that out of the way.
}

double get_rand() {
	return double(rand()) / double(RAND_MAX);
}