#include <omp.h>
#include <sys/time.h>

#include <iostream>
#include <vector>

#include "dgm.h"
#include "lib_grover.h"

int main(int argc, char **argv) {
  struct timeval timev;
  struct timeval tvBegin;
  struct timeval tvEnd;
  float t;
  int execType = t_CPU;
  int n_threads = 1;
  int cpu_region = 14;
  int cpu_coalesc = 11;
  int multi_gpu = 1;
  int gpu_region = 8;
  int gpu_coalesc = 4;
  int tam_block = 64;
  int rept = 2;
  int seed = 0;

  if (argc < 2) {
    std::cout << "You need to define the execution parameters" << std::endl;
    return 0;
  }

  int qubits = atoi(argv[1]);

  if (argc > 2) {
    execType = atoi(argv[2]);
  }

  if (argc > 3) {
    seed = atoi(argv[3]);
  }
  srand(seed);

  if (argc > 4) {
    cpu_region = atoi(argv[4]);
  }

  if (argc > 5) {
    cpu_coalesc = atoi(argv[5]);
  }

  int value = 10;
  if (argc > 6) {
    value = atoi(argv[6]);
  }

  if (execType < t_CPU || execType > t_HYBRID) {
    std::cout << "Invalid execution type: " << execType << std::endl;
    return 0;
  }

  if (execType == t_PAR_CPU || execType == t_HYBRID) {
    n_threads = omp_get_max_threads();
  } else if (execType == t_GPU && argc > 7) {
    multi_gpu = atoi(argv[7]);
  }

  std::cout << "Executing Grover: " << qubits << " qubits" << std::endl;
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
  int result = Grover(qubits, value, execType, cpu, gpu);
  gettimeofday(&tvEnd, NULL);
  timeval_subtract(&timev, &tvEnd, &tvBegin);
  t = ((float) timev.tv_sec) + ( ((float) timev.tv_usec) / 1000000.0f);

  long reg_size_per_thread = (1 << (qubits - cpu_region));
  std::cout << "reg_size_per_thread: " << reg_size_per_thread << std::endl;

  std::cout << "Time: " << t << std::endl;

  if (result != -1) {
    std::cout << "Result found successfully: " << result << std::endl;
  } else {
    std::cout << "Failed to find the value" << std::endl;
  }
}
