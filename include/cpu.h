#ifndef _CPU_H_
#define _CPU_H_

#include "common.h"
#include "pt.h"

void CpuExecution1(float complex *state, PT **pts, int qubits, int it);
void CpuExecution1_1(float complex *state, const PT *pt, long mem_size);
void CpuExecution1_2(float complex *state, const PT *pt, long mem_size);
void CpuExecution1_3(float complex *state, const PT *pt, long mem_size);

#endif