#include "lib_hadamard.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "common.h"
#include "dgm.h"
#include "gates.h"

int HadamardNQubits(int qubits, int num_of_it, int type, const CPUParams& cpu,
                    const GPUParams& gpu) {
  DGM dgm;
  dgm.qubits = qubits;
  dgm.exec_type = type;

  dgm.cpu_params = cpu;
  dgm.gpu_params = gpu;

  dgm.allocateMemory();
  dgm.setMemoryValue(0);

  std::string hadamardN = Hadamard(qubits, 0, qubits);
  dgm.setFunction(hadamardN, num_of_it);

  dgm.execute(1);

  printMem(dgm.state, 4);

  dgm.freeMemory();

  return 0;
}