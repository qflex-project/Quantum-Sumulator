#include <omp.h>
#include <sys/time.h>

#include <iostream>
#include <vector>

#include "dgm.h"
#include "lib_hadamard.h"

int main(int argc, char **argv) {
  int n_threads = 1;
  int cpu_region = 14;
  int cpu_coalesc = 11;
  int multi_gpu = 1;
  int gpu_region = 8;
  int gpu_coalesc = 4;
  int tam_block = 64;
  int rept = 2;
  int num_of_it = 3;

  if (argc < 3) {
    std::cout << "You need to define the execution parameters" << std::endl;
    return 0;
  }

  int qubits = atoi(argv[1]);

  int execType = atoi(argv[2]);
  if (execType < t_CPU || execType > t_HYBRID) {
    std::cout << "Invalid execution type: " << execType << std::endl;
    return 0;
  }

  num_of_it = atoi(argv[3]);

  if (argc > 4) {
    cpu_region = atoi(argv[4]);
  }

  if (argc > 5) {
    cpu_coalesc = atoi(argv[5]);
  }

  if (execType < t_CPU || execType > t_HYBRID) {
    std::cout << "Invalid execution type: " << execType << std::endl;
    return 0;
  }

  if (execType == t_PAR_CPU || execType == t_HYBRID) {
    n_threads = omp_get_max_threads();
  } else if (execType == t_GPU && argc > 6) {
    multi_gpu = atoi(argv[6]);
  }

  std::vector<float> amostras;

  struct timeval timev;
  struct timeval tvBegin;
  struct timeval tvEnd;
  float t;

  CPUParams cpu;
  cpu.n_threads = n_threads;
  cpu.cpu_region = cpu_region;
  cpu.cpu_coales = cpu_coalesc;

  GPUParams gpu;
  gpu.multi_gpu = multi_gpu;
  gpu.gpu_region = gpu_region;
  gpu.gpu_coales = gpu_coalesc;
  gpu.tam_block = tam_block;
  gpu.rept = rept;

  gettimeofday(&tvBegin, NULL);
  HadamardNQubits(qubits, num_of_it, execType, cpu, gpu);
  gettimeofday(&tvEnd, NULL);

  timeval_subtract(&timev, &tvEnd, &tvBegin);
  t = ((float) timev.tv_sec) + ( ((float) timev.tv_usec) / 1000000.0f);

  long reg_size_per_thread = (1 << (qubits - cpu_region));
  std::cout << "reg_size_per_thread: " << reg_size_per_thread << std::endl;

  std::cout << "Time: " << t << std::endl;

  return 0;
}