#include "gates.h"
#include "common.h"
#include <iostream>

using namespace std;

map <string, std::complex <float>*> Gates::list;

Gates::Gates(){
	init();
}

Gates::~Gates(){}

void Gates::init(){
	if (Gates::list.size() == 0){
		addGate("H", 1.0/sqrt(2), 1.0/sqrt(2), 1.0/sqrt(2), -1.0/sqrt(2));
		addGate("X", 0.0, 1.0, 1.0, 0.0);
		addGate("Z", 1.0, 0.0, 0.0, -1.0);
		addGate("Y", 0.0, std::complex <float>(0.0, 1.0), std::complex <float>(0.0, -1.0), 0.0);

		for (int k = 1; k < 3; k++){
			// Calculate the phase shift
			float exponent = 2 * M_PI / pow(2, k);
			std::complex <float> phaseShift = exp(std::complex <float>(0, exponent));

			string name = "R" + to_string(k);
			addGate(name, 1.0, 0, 0.0, phaseShift);
		}
	}
}

std::complex <float>* Gates::getMatrix(string gateName){
	if (Gates::list.find(gateName) == Gates::list.end()) return 0;

	return Gates::list[gateName];
}

bool Gates::addGate(string name, std::complex <float>* matrix){
	if (Gates::list.find(name) != Gates::list.end()) return false;

	Gates::list[name] = matrix;
	return true;
}

bool Gates::addGate(string name, std::complex <float> a0, std::complex <float> a1, std::complex <float> a2, std::complex <float> a3){
	if (Gates::list.find(name) != Gates::list.end()) return false;

	std::complex <float>* matrix = new std::complex <float>[4];
	matrix[0] = a0;
	matrix[1] = a1;
	matrix[2] = a2;
	matrix[3] = a3;

	Gates::list[name] = matrix;
	return true;
}

void Gates::printGates(){
	map <string, std::complex <float>*>::iterator it;
	for (it = Gates::list.begin(); it != Gates::list.end(); it++){
		cout << it->first << endl;
		auto matrix = it->second;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++)
				printf("%d: %.4f, %.4f\t", i*2+j, matrix[i*2+j].real(), matrix[i*2+j].imag());
			printf("\n");
		}
		if (isUnitary(it->first))
			printf("Unitary\n");
		else
			printf("Not unitary\n");

		printf ("\n");
	}
}

bool Gates::isUnitary(string gateName){
	if (Gates::list.find(gateName) == Gates::list.end()) return false;

	auto matrix = Gates::list[gateName];

	// Create a 2x2 matrix from the array
	std::complex <float> U[2][2] = {{matrix[0], matrix[1]}, {matrix[2], matrix[3]}};

	// Calculate the product U^{\dagger}U
	std::complex <float> product[2][2] = {{conj(U[0][0]) * U[0][0] + conj(U[0][1]) * U[0][1],
							conj(U[0][0]) * U[1][0] + conj(U[0][1]) * U[1][1]},
							{conj(U[1][0]) * U[0][0] + conj(U[1][1]) * U[0][1],
							conj(U[1][0]) * U[1][0] + conj(U[1][1]) * U[1][1]}};

	// Check if the product is close enough to the identity matrix
	float epsilon = 1e-6; // Adjust as needed
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			if (i == j) {
				if (abs(product[i][j] - COMPLEX_ONE) > epsilon) {
					return false;
				}
			} else {
				if (abs(product[i][j]) > epsilon) {
					return false;
				}
			}
		}
	}

	return true;
}

/////////////////////////////////////////////////////////////////////

string CNot(int qubits, int ctrl, int target, int cv){
	vector <string> cn (qubits, "ID");
	cn[ctrl] = "Control1(0)";
	if (cv) cn[ctrl] = "Control1(1)";
	cn[target] = "Target1(X)";

	return concatena(cn, qubits);
}

string Toffoli(int qubits, int ctrl1, int ctrl2, int target, int cv){
	vector <string> tf (qubits, "ID");
	tf[ctrl1] = "Control1(0)";
	if (cv>>1) tf[ctrl1] = "Control1(1)";
	tf[ctrl2] = "Control1(0)";
	if (cv&1) tf[ctrl2] = "Control1(1)";
	tf[target] = "Target1(X)";

	return concatena(tf, qubits);
}

string Controlled1(int qubits, int ctrl, int target, string op, int cv){
	vector <string> c1 (qubits, "ID");
	c1[ctrl] = "Control1(0)";
	if (cv) c1[ctrl] = "Control1(1)";
	c1[target] = "Target1(" + op + ")";

	return concatena(c1, qubits);
}


string Controlled2(int qubits, int ctrl1, int ctrl2, int target, string op, int cv){
	vector <string> c2 (qubits, "ID");
	c2[ctrl1] = "Control1(0)";
	if (cv&2) c2[ctrl1] = "Control1(1)";
	c2[ctrl2] = "Control1(0)";
	if (cv&1) c2[ctrl2] = "Control1(1)";
	c2[target] = "Target1(" + op + ")";

	return concatena(c2, qubits);
}

string ControlledN(int qubits, vector <int> ctrls, int target, string op, int cv){
	if (cv == -1) cv = pow(2,ctrls.size()) - 1;
	vector <string> c (qubits, "ID");

	for (int i = ctrls.size() - 1; i >=0; i--){
		c[ctrls[i]] = "Control1(0)";
		if (cv & 1) c[ctrls[i]] = "Control1(1)";
		cv = cv >> 1;
	}

	c[target] = "Target1(" + op + ")";

	return concatena(c, qubits);
}


string Pauli_X(int qubits, int reg, int width){
	vector <string> px (qubits, "ID");
	for (int i = 0; i < width; i++) px[i+reg] = "X";

	return concatena(px, qubits);
}

string Pauli_Z(int qubits, int reg, int width){
	vector <string> px (qubits, "ID");
	for (int i = 0; i < width; i++) px[i+reg] = "Z";

	return concatena(px, qubits);
}


string Hadamard(int qubits, int reg, int width){
	vector <string> hn (qubits, "ID");
	for (int i = 0; i < width; i++) hn[i+reg] = "H";

	return concatena(hn, qubits);
}

string NumberAsX(int qubits, int num){
	vector <string> c (qubits, "ID");
	for (int i = 0; i < qubits; i++)
		if (num >> ( qubits - i - 1) & 1) c[i] = "X";

	return concatena(c, qubits);
}

string concatena(vector <string> vec, int size, bool rev){
	string s;
	if (!rev){
		s = vec[0];
		for (int i = 1; i < size; i++)
			s += "," + vec[i];
	}
	else{
		s = vec[size-1];
		for (int i = size - 2; i >= 0; i--)
			s += "," + vec[i];
	}

	return s;
}
