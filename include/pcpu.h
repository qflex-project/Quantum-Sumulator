#ifndef _PCPU_H_
#define _PCPU_H_

#include "common.h"
#include "pt.h"

void PCpuExecution1(float complex *state, PT **pts, int qubits, long n_threads,
                    int coales, int region, int it);
void PCpuExecution1_0(float complex *state, PT **pts, int qubits, int start,
                      int end, int pos_count, int reg_id, int reg_mask);

#endif