#include "lib_grover.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "common.h"
#include "dgm.h"
#include "gates.h"

int Grover(int qubits, int value, int type, const CPUParams& cpu,
           const GPUParams& gpu) {
  DGM dgm;
  dgm.qubits = qubits;
  dgm.exec_type = type;

  dgm.cpu_params = cpu;
  dgm.gpu_params = gpu;

  dgm.allocateMemory();
  dgm.setMemoryValue(1 << (qubits - 1));

  std::string H = Hadamard(qubits, 0, qubits);
  std::string orcl = Oracle1(qubits, value);
  std::string CZ = ControledZ(qubits);

  std::vector<std::string> grover_step;

  grover_step.push_back(orcl);

  for (int i = 1; i < qubits; i++) {
    grover_step.push_back(Hadamard(qubits, i, 1));
    grover_step.push_back(Pauli_X(qubits, i, 1));
  }

  grover_step.push_back(CZ);

  for (int i = qubits - 1; i >= 1; i--) {
    grover_step.push_back(Pauli_X(qubits, i, 1));
    grover_step.push_back(Hadamard(qubits, i, 1));
  }

  int num_of_it = (int)(M_PI / 4.0 * sqrt(1 << (qubits - 1)));
  int result = 0;

  dgm.setFunction(H);
  dgm.setFunction(grover_step, num_of_it, false);
  dgm.execute(1);

  for (int i = 1; i < qubits; i++) {
    result = (result << 1) | dgm.measure(i);
  }

  dgm.freeMemory();

  return result;
}

std::string Oracle1(int qubits, int value) {
  std::vector<std::string> t(qubits);
  int ctrl_v;

  for (int i = qubits - 1; i >= 1; i--) {
    ctrl_v = value & 1;
    value = value >> 1;
    if (ctrl_v)
      t[i] = "Control1(1)";
    else
      t[i] = "Control1(0)";
  }
  t[0] = "Target1(X)";

  return concatena(t, qubits);
}

std::string ControledZ(int qubits) {
  std::vector<std::string> cz;
  cz.push_back("ID");
  for (int i = 0; i < qubits - 2; i++) cz.push_back("Control1(1)");
  cz.push_back("Target1(Z)");

  return concatena(cz, qubits);
}
