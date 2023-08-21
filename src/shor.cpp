#include <omp.h>
#include <sys/time.h>

#include <iostream>
#include <map>
#include <vector>

#include "common.h"
#include "dgm.h"
#include "lib_shor.h"

int main(int argc, char** argv) {
  std::map<int, int> qubitsMap;
  qubitsMap[15] = 57;
  qubitsMap[17] = 119;
  qubitsMap[19] = 253;
  qubitsMap[21] = 485;
  qubitsMap[23] = 1017;
  qubitsMap[25] = 2045;
  qubitsMap[27] = 2863;

  struct timeval timev;
  struct timeval tvBegin;
  struct timeval tvEnd;
  float t;
  std::vector<float> amostras;

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
  if (qubitsMap.count(qubits) == 0) {
    std::cout << "The amount of qubits does not map to a valid number to be "
                 "factored: "
              << qubits << std::endl;
    return 0;
  }

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

  if (execType < t_CPU || execType > t_HYBRID) {
    std::cout << "Invalid execution type: " << execType << std::endl;
    return 0;
  }

  if (execType == t_PAR_CPU || execType == t_HYBRID) {
    n_threads = omp_get_max_threads();
  } else if (execType == t_GPU && argc > 6) {
    multi_gpu = atoi(argv[6]);
  }

  std::vector<int> factors;

  std::cout << "Executing Shor: " << qubits << " qubits" << std::endl;

  gettimeofday(&tvBegin, NULL);
  factors =
      Shor(qubitsMap[qubits], execType, n_threads, cpu_region, cpu_coalesc,
           multi_gpu, gpu_region, gpu_coalesc, tam_block, rept);
  gettimeofday(&tvEnd, NULL);
  timeval_subtract(&timev, &tvEnd, &tvBegin);
  t = ((float) timev.tv_sec) + ( ((float) timev.tv_usec) / 1000000.0f);

  long reg_size_per_thread = (1 << (qubits - cpu_region));
  std::cout << "reg_size_per_thread: " << reg_size_per_thread << std::endl;

  std::cout << "Time: " << t << std::endl;

  if (factors.size() == 2) {
    std::cout << "Found factors: " << factors[0] << " -- " << factors[1]
              << std::endl;
  } else {
    std::cout << "Failed to find factors" << std::endl;
  }

  return 0;
}
