#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <complex>

#define CHUNCK_SIZE 262144

#define PT_TAM 1
#define QB_LIMIT 30

#define DIM_BLOCK 2048

#define TAM_ARG 5

#define SHIFT 0
#define CTRL_MASK 1
#define CTRL_VALUE 2
#define CTRL_REG_MASK 3
#define CTRL_REG_VALUE 4

#define ACUMM 0
#define SHIFT_READ 0
#define SHIFT_WRITE 0
#define MAT_START 0
#define MAT_SIZE 0
#define MAT_END 0

static const std::complex <float> COMPLEX_ZERO = std::complex <float>(0.0, 0.0);
static const std::complex <float> COMPLEX_ONE = std::complex <float>(1.0, 0.0);
static const std::complex <float> COMPLEX_I = std::complex <float>(0.0, 1.0); // Imaginary unit

static const float M_PI = std::acos(-1.0); // Calculate Pi
static const float M_E = std::exp(1.0); // Calculate Euler's number (e)

enum {
	DENSE,
	DIAG_PRI,
	DIAG_SEC
};

struct PT{
	int qubits;
	std::complex <float> *matrix;
	int mat_size;
	int start, end;
	bool affected;
	
	long ctrl_value, ctrl_mask;
	long *ctrl_pos, ctrl_count;
	long *ctrl_rest, ctrl_rest_count;

	PT();

	void destructor();
	long ctrlAffect(long qubit);
	long matrixType();
	void setArgs(long *arg, long affect);
	void setArgs_soft(long *arg, long affect);
	void setArgsGPU(long *arg, int region_start, int region_size, int coalesc);
	void print();
	void printMatrix();

};

void printMem(std::complex <float>* mem, int qubits);
void printMemExp(std::complex <float>* mem, int qubits, int reg1, int reg2, long n);
void printMemCheckExp(std::complex <float>* mem, int qubits, long n, long a, long N);

long modular_pow(long base, long exponent, long modulus);

bool increasing(const PT *pt1, const PT *pt2);
bool decreasing(const PT *pt1, const PT *pt2);

void swap_value(int *v1, int* v2);
void swap_ptr(float **ptr1, float **ptr2);
void swap_ptr(std::complex <float> **ptr1, std::complex <float> **ptr2);

int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1);

bool is_valid_quantum_state(std::complex <float>* state, int qubits);

#endif

