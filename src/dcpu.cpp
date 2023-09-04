#include "dcpu.h"

#include <math.h>
#include <mpi.h>

// Coalescimento
CoalescResult projectState(const float complex *state, int qubits,
                           int proj_qubits, long reg_id, long reg_mask,
                           int world_size) {
  CoalescResult result;
  int qbs_coales = 0;
  for (int i = 0; i < qubits; i++) {
    if ((reg_mask >> i) & 1) {
      qbs_coales++;
    } else {
      break;
    }
  }

  long mem_portions = long_pow(2.0, proj_qubits - qbs_coales);
  int portion_size = 1 << qbs_coales;

  result.new_state =
      (float complex *)(malloc(sizeof(float complex) * pow(2, qubits)));
  result.chunk_sizes = (int *)(malloc(sizeof(int) * world_size));
  result.displ = (int *)(malloc(sizeof(int) * world_size));

  long inc = ~(reg_mask >> qbs_coales);

  long dest_pos = 0;
  long src_pos = 0;
  long base = 0;
  printf("projectState: %d\n", reg_id);
  for (long d = 0; d < world_size; d++) {
    result.displ[d] = dest_pos;
    printf("displ[%d]: %d\n", d, result.displ[d]);
    for (long b = mem_portions / world_size * d;
         b < mem_portions / world_size * (d + 1); b++) {
      src_pos = (base << qbs_coales) | reg_id;

      memcpy(result.new_state + dest_pos, state + src_pos,
             portion_size * sizeof(float complex));

      base = (base + inc + 1) & ~inc;
      dest_pos += portion_size;
    }
    result.chunk_sizes[d] = dest_pos - result.displ[d];
    printf("chunk_sizes[%d]: %d\n", d, result.chunk_sizes[d]);

  }
  return result;
}

bool collectState(float complex *state, CoalescResult &r, int qubits,
                  int proj_qubits, long reg_id, long reg_mask, int world_size) {
  int qbs_coales = 0;
  for (int i = 0; i < qubits; i++) {
    if ((reg_mask >> i) & 1)
      qbs_coales++;
    else
      break;
  }

  int mem_portions = pow(2.0, proj_qubits - qbs_coales);
  int portion_size = 1 << qbs_coales;

  long inc = ~(reg_mask >> qbs_coales);

  long dev_pos, pos, base = 0;
  for (int d = 0; d < world_size; d++) {
    dev_pos = 0;
    for (int b = mem_portions / world_size * d;
         b < mem_portions / world_size * (d + 1); b++) {
      pos = (base << qbs_coales) | reg_id;

      memcpy(state + pos, r.new_state + dev_pos,
             portion_size * sizeof(float complex));

      base = (base + inc + 1) & ~inc;
      dev_pos += portion_size;
    }
  }

  free(r.new_state);

  return true;
}