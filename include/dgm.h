#ifndef _DGM_H_
#define _DGM_H_

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "common.h"
#include "gates.h"
#include "pt.h"

float complex *GenericExecute(float complex *state, std::string function,
                              int qubits, int type, int threads, int factor);
float complex *GenericExecute(float complex *state,
                              std::vector<std::string> function, int qubits,
                              int type, int threads, int factor);

inline long LINE(long pos, long shift) { return ((pos >> shift) & 1) * 2; }
inline long BASE(long pos, long shift) { return pos & (~(1 << shift)); }

enum ExecutionType { t_CPU, t_PAR_CPU, t_GPU, t_HYBRID, t_SPEC };

class Group {
 public:
  std::vector<std::string> ops;
  std::vector<long> pos_ops;
  std::vector<bool> ctrl;
  std::vector<long> pos_ctrl;

  Group(){};
  bool isAfected(int pos, int afect);
};

struct OPSCounter {
  long total_op = 0;
  long dense = 0;
  long main_diag = 0;
  long sec_diag = 0;
  long c_dense = 0;
  long c_main_diag = 0;
  long c_sec_diag = 0;
};

struct CPUParams {
  long n_threads;
  long cpu_coales;
  long cpu_region;
};

struct GPUParams {
  int multi_gpu;
  long gpu_coales;
  long gpu_region;
  int tam_block;
  int rept;
};

class DGM {
 public:
  std::vector<std::string> diag;
  long MAX_QB;
  long MAX_PT;
  long qb_afected;

  long factor;

  int exec_type;

  CPUParams cpu_params;
  GPUParams gpu_params;

  std::vector<PT *> vec_pts;
  PT **pts;
  long qubits;

  float measure_value;

  float elapsed_time;
  struct timeval timev;

  float complex *state;

  DGM();
  ~DGM();

  bool en_print;

  void printPTs();
  void erase();
  void setExecType(int type);

  void setCpuStructure(long cpu_region, long cpu_coales);
  void setGpuStructure(long gpu_coales, long gpu_region, int rept = 1);

  void allocateMemory();
  void setMemory(float complex *mem);
  void freeMemory();
  void setMemoryValue(int pos);

  int measure(int q_pos);
  std::map<long, float> measure(std::vector<int> q_pos);
  void colapse(int q_pos, int value);

  void setFunction(std::string function, int it = 1, bool er = true);
  void setFunction(std::vector<std::string> steps, int it = 1, bool er = true);
  std::map<long, Group> genGroups(std::string step);
  void genPTs(std::map<long, Group> &gps, std::vector<PT *> &step_pts);
  void genMatrix(float complex *matrix, std::vector<float complex *> &matrices,
                 long tam, long current, long line, long column,
                 float complex cmplx);

  OPSCounter CountOps(int it = 1);

  void executeFunction(std::string function, int it = 1);
  void executeFunction(std::vector<std::string> steps, int it = 1);
  float complex *execute(int it);

  void HybridExecution(PT **pts);

  void CpuExecution1(int it);
  void CpuExecution1_1(PT *pt, long mem_size);
  void CpuExecution1_2(PT *pt, long mem_size);
  void CpuExecution1_3(PT *pt, long mem_size);
};

#endif
