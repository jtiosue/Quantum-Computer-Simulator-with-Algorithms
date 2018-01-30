#ifndef QALGORITHMS_INCLUDE
#define QALGORITHMS_INCLUDE

#include "quantum.h"

// Useful functions

unsigned int char_to_int(char c);

unsigned int gcd(unsigned int a, unsigned int b);

unsigned int binary_to_base10(string s);

string base10_to_binary(unsigned int x);



// Qalgorithms

void QFT(Register *reg, unsigned int start=0, unsigned int end=0);

void IQFT(Register *reg, unsigned int start=0, unsigned int end=0);

void Add(Register *reg);

unsigned int Add(unsigned int a, unsigned int b);

void ModAdd(Register *reg);

unsigned int ModAdd(unsigned int a, unsigned int b, unsigned int N);

unsigned int Grover(unsigned int omega, unsigned int num_bits, bool verbose=true);

#endif