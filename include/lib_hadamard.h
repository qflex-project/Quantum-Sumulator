#ifndef _LIBGENERAL_H_
#define _LIBGENERAL_H_

#include "dgm.h"

#define complex _Complex


int HadamardNQubits(int qubits, int num_of_it, int type, const CPUParams& cpu,
                    const GPUParams& gpu);

#endif