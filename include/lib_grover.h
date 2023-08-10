#ifndef _LIBGROVER_H_
#define _LIBGROVER_H_

#include <string>
#include <vector>

#define complex _Complex

int Grover(int qubits, int value, int type, const CPUParams& cpu,
           const GPUParams& gpu);

std::string ControledZ(int qubits);
std::string Oracle1(int qubits, int value);

#endif
