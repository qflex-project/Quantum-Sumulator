#ifndef _DCPU_H_
#define _DCPU_H_

#include "common.h"
#include "pt.h"

typedef struct CoalescResult {
  float complex *new_state;
  long *chunk_sizes;
  long *displ;
} CoalescResult;

void DCpuExecution1(float complex *state, PT **pts, int qubits, int procs,
                    int it);
CoalescResult projectState(const float complex *state, int qubits,
                           int proj_qubits, long reg_id, long reg_mask,
                           int world_size);
bool collectState(float complex *state, CoalescResult &r, int qubits,
                  int proj_qubits, long reg_id, long reg_mask, int world_size);
#endif