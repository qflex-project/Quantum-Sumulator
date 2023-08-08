#ifndef _LIBGROVER_H_
#define _LIBGROVER_H_

#include <string>
#include <vector>

#define complex _Complex

float Grover(long qubits, long value, int type, int n_threads, int cpu_region,
             int cpu_coales, int multi_gpu, int gpu_region, int gpu_coales,
             int tam_block, int rept);

std::string ControledZ(int qubits);
std::string Oracle1(long qubits, long int value);

#endif
