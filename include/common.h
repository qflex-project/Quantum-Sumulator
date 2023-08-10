#ifndef _COMMON_H_
#define _COMMON_H_

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define complex _Complex

const int PT_TAM = 1;
const int QB_LIMIT = 30;

const int TAM_ARG = 5;

const int SHIFT = 0;
const int CTRL_MASK = 1;
const int CTRL_VALUE = 2;
const int CTRL_REG_MASK = 3;
const int CTRL_REG_VALUE = 4;

enum MatrixType { DENSE, DIAG_PRI, DIAG_SEC };

void printMem(const float complex *mem, const int qubits);
void printMemExp(const float complex *mem, const int qubits, const int reg1,
                 const int reg2, const long n);
void printMemCheckExp(const float complex *mem, const int qubits, const long n,
                      const long a, const long N);

long modular_pow(long base, long exponent, long modulus);
long long_pow(long base, long exponent);

void swap_value(int *v1, int *v2);
void swap_ptr(float **ptr1, float **ptr2);
void swap_ptr(float complex **ptr1, float complex **ptr2);

int timeval_subtract(struct timeval *result, const struct timeval *t2,
                     const struct timeval *t1);

#endif
