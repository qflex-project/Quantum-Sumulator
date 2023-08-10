#ifndef _PT_H_
#define _PT_H_

#include "common.h"

#define complex _Complex

struct PT {
  int qubits;

  float complex *matrix;
  int mat_size;

  int start;
  int end;

  bool affected;

  long ctrl_value;
  long ctrl_mask;

  long *ctrl_pos;
  long ctrl_count;

  long *ctrl_rest;
  long ctrl_rest_count;

  PT() : matrix(NULL), ctrl_pos(NULL), ctrl_rest(NULL){};

  void destructor();
  long ctrlAffect(long qubit) const;
  long matrixType() const;
  void setArgs(long *arg, long affect) const;
  void setArgs_soft(long *arg, long affect) const;
  void setArgsGPU(long *arg, int region_start, int region_size,
                  int coalesc) const;
  void print() const;
  void printMatrix() const;
};

bool increasing(const PT *pt1, const PT *pt2);
bool decreasing(const PT *pt1, const PT *pt2);

#endif
