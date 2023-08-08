#include "common.h"

#include <math.h>

const float round_precision_float = 0.00000000001f;

void printMem(const float complex *mem, const int qubits) {
  long size = long_pow(2, qubits);

  float real = 0.0;
  float imag = 0.0;
  for (long i = 0; i < size; i++) {
    real = crealf(mem[i]);
    imag = cimagf(mem[i]);
    printf("%ld:\t%.6f %.6f\n", i, real, imag);
  }
}

void printMemExp(const float complex *mem, const int qubits, const int reg1,
                 const int reg2, const long n) {
  long size = long_pow(2, qubits);
  long mask = long_pow(2, n) - 1;

  float real = 0.0;
  float imag = 0.0;
  float f = 0.0;
  long last_X = 0;
  for (long i = 0; i < size; i++) {
    real = imag = 0;
    f = fabs(crealf(mem[i]));
    if (f > round_precision_float) {
      real = crealf(mem[i]);
    }

    f = fabs(cimagf(mem[i]));
    if (f > round_precision_float) {
      imag = cimagf(mem[i]);
    }

    if ((imag != 0) || (real != 0)) {
      long X = (i >> (qubits - reg1 - n)) & mask;
      long Exp = (i >> (qubits - reg2 - n)) & mask;
      printf("%ld\t>>  X: %ld\tExp: %ld\tDif: %ld\t\t\tV: %f %f\n", i, X, Exp,
             X - last_X, crealf(mem[i]), cimagf(mem[i]));
      last_X = X;
    }
  }
}

void printMemCheckExp(const float complex *mem, const int qubits,
                      const long width, const long a, const long N) {
  long size = long_pow(2, qubits);

  long mask = long_pow(2, width) - 1;

  float real = 0.0;
  float imag = 0.0;
  float f = 0.0;
  long last_X = 0;
  for (long i = 0; i < size; i++) {
    real = imag = 0;
    f = fabs(crealf(mem[i]));
    if (f > round_precision_float) {
      real = crealf(mem[i]);
    }

    f = fabs(cimagf(mem[i]));
    if (f > round_precision_float) {
      imag = cimagf(mem[i]);
    }

    if ((imag != 0) || (real != 0)) {
      long X = (i >> (2 * width + 2));
      long Exp = (i >> (width + 2)) & mask;

      printf("%ld\t>>  X: %ld\tExp: %ld\tDif: %ld\t\t", i, X, Exp, X - last_X);
      last_X = X;
      if (modular_pow(a, X, N) != Exp) {
        printf("Errado\n");
      } else {
        printf("\n");
      }
    }
  }
}

////////////////////////////////////////////////////////////
void swap_value(int *v1, int *v2) {
  int aux = *v1;
  *v1 = *v2;
  *v2 = aux;
}

void swap_ptr(float **ptr1, float **ptr2) {
  float *aux = *ptr1;
  *ptr1 = *ptr2;
  *ptr2 = aux;
}

void swap_ptr(float complex **ptr1, float complex **ptr2) {
  float complex *aux = *ptr1;
  *ptr1 = *ptr2;
  *ptr2 = aux;
}
//////////////////////////////////////////////////////////

int timeval_subtract(struct timeval *result, const struct timeval *t2,
                     const struct timeval *t1) {
  long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) -
                  (t1->tv_usec + 1000000 * t1->tv_sec);
  result->tv_sec = diff / 1000000;
  result->tv_usec = diff % 1000000;

  return (diff < 0);
}

long modular_pow(long base, long exponent, long modulus) {
  long result = 1;
  base = base % modulus;
  while (exponent) {
    if (exponent & 1) result = (result * base) % modulus;
    exponent = exponent >> 1;
    base = (base * base) % modulus;
  }
  return result;
}

long long_pow(long base, long exponent) { return (long)pow(base, exponent); }