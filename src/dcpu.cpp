#include "dcpu.h"

#include <math.h>
#include <mpi.h>

void DCpuExecution1(float complex *state, PT **pts, int qubits, int procs,
                    int it) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // We are assuming at least 2 processes for this task
  if (world_size < 2) {
    fprintf(stderr, "World size must be greater than 1 for this to work!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  } else {
    fprintf(stderr, "Its all right!\n");
  }

  MPI_Finalize();
}
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
  result.chunk_sizes = (long *)(malloc(sizeof(long) * world_size));
  result.displ = (long *)(malloc(sizeof(long) * world_size));

  long inc = ~(reg_mask >> qbs_coales);

  long dest_pos = 0;
  long src_pos = 0;
  long base = 0;
  for (long d = 0; d < world_size; d++) {
    result.displ[d] = dest_pos;
    for (long b = mem_portions / world_size * d;
         b < mem_portions / world_size * (d + 1); b++) {
      src_pos = (base << qbs_coales) | reg_id;

      memcpy(result.new_state + dest_pos, state + src_pos,
             portion_size * sizeof(float complex));

      base = (base + inc + 1) & ~inc;
      dest_pos += portion_size;
    }
    result.chunk_sizes[d] = dest_pos - result.displ[d];
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