#ifndef _HYBRID_H_
#define _HYBRID_H_

#include "common.h"
#include "pt.h"


void HybridExecution(float complex *state, PT **pts, int qubits,
                     const CPUParams& cpu_params, const GPUParams& gpu_params);


#endif