#ifndef _GATES_H_
#define _GATES_H_

#include <complex.h>
#include <math.h>

#include <map>
#include <string>
#include <vector>

#include "common.h"

std::string CNot(int qubits, int ctrl, int target, int cv = 1);
std::string Toffoli(int qubits, int ctrl1, int ctrl2, int target, int cv = 3);
std::string Controlled1(int qubits, int ctrl, int target, const std::string& op,
                        int cv = 1);
std::string Controlled2(int qubits, int ctrl1, int ctrl2, int target,
                        const std::string& op, int cv = 3);
std::string ControlledN(int qubits, std::vector<int> ctrls, int target,
                        const std::string& op, int cv = -1);
std::string Pauli_X(int qubits, int reg, int width = 1);
std::string Pauli_Z(int qubits, int reg, int width = 1);
std::string Hadamard(int qubits, int reg, int width = 1);

std::string concatena(std::vector<std::string> vec, int size, bool rev = false);

class Gates {
 public:
  // possivel vazamento de memória nessa lista estatica
  static std::map<std::string, float complex*> list;
  Gates();
  ~Gates();
  void init() const;
  float complex* getMatrix(const std::string& gateName);
  bool addGate(const std::string& name, float complex* matrix) const;
  bool addGate(const std::string& name, float complex a0, float complex a1,
               float complex a2, float complex a3) const;
};

#endif
