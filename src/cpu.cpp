#include "cpu.h"
#include <math.h>

void CpuExecution1(float complex *state, PT **pts, int qubits, int it) {
  long mem_size = long_pow(2.0, qubits);

  for (int x = 0; x < it; x++) {
    long i = 0;
    while (pts[i] != NULL) {
      long mt = pts[i]->matrixType();

      switch (mt) {
        case DENSE:
          CpuExecution1_1(state, pts[i], mem_size);
          break;
        case DIAG_PRI:
          CpuExecution1_2(state, pts[i], mem_size);
          break;
        case DIAG_SEC:
          CpuExecution1_3(state, pts[i], mem_size);
          break;
        default:
          exit(1);
      }
      i++;
    }
  }
}

void CpuExecution1_1(float complex *state, const PT *pt, long mem_size) {  // Denso
  long pos0;
  long pos1;
  long shift;

  shift = 1 << pt->end;

  float complex tmp;

  if (!pt->ctrl_count) {  // operador não controlado
    mem_size /= 2;
    for (long pos = 0; pos < mem_size; pos++) {
      pos0 = (pos * 2) - (pos & (shift - 1));
      pos1 = pos0 | shift;

      tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
      state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];
      state[pos0] = tmp;
    }
  } else {  // operador controlado
    long mask = ~(pt->ctrl_mask | shift);
    long inc = (~mask) + 1;

    for (long pos = 0; pos < mem_size; pos = (pos + inc) & mask) {
      pos0 = pos | pt->ctrl_value;
      pos1 = pos0 | shift;

      tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
      state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];
      state[pos0] = tmp;
    }
  }
}

void CpuExecution1_2(float complex *state, const PT *pt, long mem_size) {  // Diagonal Principal
  long pos0;
  long shift = pt->end;

  if (!pt->ctrl_count)  // operador não controlado
    for (long pos = 0; pos < mem_size; pos++)
      state[pos] = pt->matrix[((pos >> shift) & 1) * 3] * state[pos];
  else {  // operador controlado
    long mask = ~(pt->ctrl_mask);
    long inc = (~mask) + 1;

    for (long pos = 0; pos < mem_size; pos = (pos + inc) & mask) {
      pos0 = pos | pt->ctrl_value;

      state[pos0] = pt->matrix[((pos0 >> shift) & 1) * 3] * state[pos0];
    }
  }
}

void CpuExecution1_3(float complex *state, const PT *pt, long mem_size) {  // Diagonal Secundária
  long pos0;
  long pos1;
  long shift;

  shift = 1 << pt->end;

  float complex tmp;

  if (!pt->ctrl_count) {  // operador não controlado
    mem_size /= 2;
    for (long pos = 0; pos < mem_size; pos++) {
      pos0 = (pos * 2) - (pos & (shift - 1));
      pos1 = pos0 | shift;

      tmp = pt->matrix[1] * state[pos1];
      state[pos1] = pt->matrix[2] * state[pos0];
      state[pos0] = tmp;
    }
  } else {  // operador controlado
    long mask = ~(pt->ctrl_mask | shift);
    long inc = (~mask) + 1;

    for (long pos = 0; pos < mem_size; pos = (pos + inc) & mask) {
      pos0 = pos | pt->ctrl_value;
      pos1 = pos0 | shift;

      tmp = pt->matrix[1] * state[pos1];
      state[pos1] = pt->matrix[2] * state[pos0];
      state[pos0] = tmp;
    }
  }
}