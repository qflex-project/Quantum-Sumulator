#ifndef _COMMON_H_
#define _COMMON_H_

#define complex _Complex

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sys/time.h>

const int PT_TAM = 1;
const int QB_LIMIT = 30;

const int TAM_ARG = 5;

const int SHIFT = 0;
const int CTRL_MASK = 1;
const int CTRL_VALUE = 2;
const int CTRL_REG_MASK = 3;
const int CTRL_REG_VALUE = 4;

enum MatrixType { DENSE, DIAG_PRI, DIAG_SEC };

enum ExecutionType { t_CPU, t_PAR_CPU, t_GPU, t_HYBRID, t_SPEC };

class Group {
 public:
  std::vector<std::string> ops;
  std::vector<long> pos_ops;
  std::vector<bool> ctrl;
  std::vector<long> pos_ctrl;

  Group(){};
  bool isAfected(int pos, int afect);
};

struct OPSCounter {
  int total_op;
  int dense;
  int main_diag;
  int sec_diag;
  int c_dense;
  int c_main_diag;
  int c_sec_diag;
};

struct CPUParams {
  int n_threads;
  int cpu_coales;
  int cpu_region;
};

struct GPUParams {
  int multi_gpu;
  int gpu_coales;
  int gpu_region;
  int tam_block;
  int rept;
};

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
