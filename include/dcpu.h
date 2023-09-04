#ifndef _DCPU_H_
#define _DCPU_H_

#include "common.h"
#include "pt.h"

typedef struct CoalescResult {
  float complex *new_state;
  int *chunk_sizes;
  int *displ;
} CoalescResult;

CoalescResult projectState(const float complex *state, int qubits,
                           int proj_qubits, long reg_id, long reg_mask,
                           int world_size);
bool collectState(float complex *state, CoalescResult &r, int qubits,
                  int proj_qubits, long reg_id, long reg_mask, int world_size);
#endif