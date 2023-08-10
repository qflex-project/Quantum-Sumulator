#ifndef _LIBSHOR_H_
#define _LIBSHOR_H_

#include <string>
#include <vector>

#define complex _Complex

// N - Number to ne factored
// type - Execution Type
// threads - Number of threads to be used in case of a parallel execution on CPU
std::vector<int> Shor(long N, int type, int n_threads, int cpu_region,
                      int cpu_coalesc, int multi_gpu, int gpu_region,
                      int gpu_coalesc, int tam_block, int rept);

void ApplyQFT(int qubits, int type, int multi_gpu, int qbs_region, int coalesc,
              int tam_block, int rept);

//////////////////////////////////////////////////////////////////////////

std::vector<std::string> QFT(int qubits, int reg, int over, int width);
std::vector<std::string> QFT2(int qubits, int reg, int width);
std::vector<std::string> RQFT(int qubits, int reg, int over, int width);
std::vector<std::string> CSwapR(int qubits, int ctrl, int reg1, int reg2,
                                int width);
std::vector<std::string> SwapOver(int qubits, int reg, int width);
std::string genRot(int qubits, int reg, int value);

std::vector<std::string> CMultMod(int qubits, int ctrl, int reg1, int reg2,
                                  int over, int over_bool, int width, long a,
                                  long N);
std::vector<std::string> CRMultMod(int qubits, int ctrl, int reg1, int reg2,
                                   int over, int over_bool, int width, long a,
                                   long N);
std::vector<std::string> C2AddMod(int qubits, int ctrl1, int ctrl2, int reg,
                                  int over, int over_bool, int width, long a,
                                  long N);
std::vector<std::string> C2SubMod(int qubits, int ctrl1, int ctrl2, int reg,
                                  int over, int over_bool, int width, long a,
                                  long N);

std::string C2AddF(int qubits, int ctrl1, int ctrl2, int reg, int over,
                   long num, int width);
std::string CAddF(int qubits, int ctrl1, int reg, int over, long num,
                  int width);
std::string AddF(int qubits, int reg, int over, long num, int width);
std::vector<std::string> AddF(int qubits, int reg, int over, long num,
                              int width, bool controlled);

std::string C2SubF(int qubits, int ctrl1, int ctrl2, int reg, int over,
                   long num, int width);
std::string CSubF(int qubits, int ctrl1, int reg, int over, long num,
                  int width);
std::string SubF(int qubits, int reg, int over, long num, int width);
std::vector<std::string> SubF(int qubits, int reg, int over, long num,
                              int width, bool controlled);

//////////////////////////////////////////////////////////////////////////

std::string int2str(int number);

long mul_inv(long a, long b);
int revert_bits(int res, int n);
int quantum_ipow(int a, int b);

/* Calculate the greatest common divisor with Euclid's algorithm */
int quantum_gcd(int u, int v);
void quantum_frac_approx(int *a, int *b, int width);

//////////////////////////////////////////////////////////////////////////

#endif
