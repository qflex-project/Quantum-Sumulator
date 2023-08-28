#ifndef _PCPU_H_
#define _PCPU_H_

#include "common.h"
#include "pt.h"

void PCpuExecution1(float complex *state, PT **pts, e_size qubits, int n_threads,
                    e_size coales, e_size region);
void PCpuExecution1_0(float complex *state, PT **pts, e_size qubits, e_size start,
                      e_size end, e_size pos_count, e_size reg_id, e_size reg_mask);

#endif