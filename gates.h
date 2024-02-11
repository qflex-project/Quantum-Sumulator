#ifndef _GATES_H_
#define _GATES_H_

#include "common.h"
#include <vector>
#include <map>
#include <string>
#include <complex.h>

#define _USE_MATH_DEFINES

#define ccomplex _Complex

using namespace std;

string CNot(int qubits, int ctrl, int target, int cv = 1);
string Toffoli(int qubits, int ctrl1, int ctrl2, int target, int cv = 3);
string Controlled1(int qubits, int ctrl, int target, string op, int cv = 1);
string Controlled2(int qubits, int ctrl1, int ctrl2, int target, string op, int cv = 3);
string ControlledN(int qubits, vector<int> ctrls, int target, string op, int cv = -1);
string Pauli_X(int qubits, int reg, int width = 1);
string Pauli_Z(int qubits, int reg, int width = 1);
string Hadamard(int qubits, int reg, int width = 1);

string concatena(vector <string> vec, int size, bool rev = false);


class Gates{
public:
	static map <string, std::complex <float>*> list;
	Gates();
	~Gates();
	void init();
	std::complex <float>* getMatrix(string gateName);
	bool addGate(string name, std::complex <float>* matrix);
	bool addGate(string name, std::complex <float> a0, std::complex <float> a1, std::complex <float> a2, std::complex <float> a3);
	void printGates();
};

#endif
