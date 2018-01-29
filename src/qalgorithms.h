#ifndef QALGORITHMS_INCLUDE
#define QALGORITHMS_INCLUDE

#include "quantum.h"

// Useful functions

int char_to_int(char c);

int gcd(unsigned int a, unsigned int b);

int binary_to_base10(string s);



// Qalgorithms

void QFT(Register *reg, int start=0, int end=0);

void IQFT(Register *reg, int start=0, int end=0);

int Shor(unsigned int N, unsigned int depth_limit=10);

void Add(Register *reg);

int Add(unsigned int a, unsigned int b);

void ModAdd(Register *reg);

int ModAdd(unsigned int a, unsigned int b, unsigned int N);

int Grover(unsigned int omega, unsigned int num_bits, bool verbose=true);

int function(int x);

#endif