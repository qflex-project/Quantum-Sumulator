#ifndef _GPU_H_
#define _GPU_H_

#include "common.h"
#include "pt.h"

extern "C" bool setDevice(int num);

extern "C" float complex *GpuExecutionWrapper(float complex *state, PT **pts,
                                              int qubits, int multi_gpu,
                                              int coalesc, int qbs_region,
                                              int tam_block, int rept,
                                              int num_it);
extern "C" bool ProjectState(float complex *state, int qubits, int region_size,
                             long reg_id, long reg_mask, int multi_gpu);
extern "C" bool GetState(float complex *state, int qubits, int region_size,
                         long reg_id, long reg_mask, int multi_gpu);

#endif