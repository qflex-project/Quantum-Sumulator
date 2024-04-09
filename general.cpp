#include "lib_general.h"
#include "lib_grover.h"
#include "lib_shor.h"
#include "dgm.h"
#include <vector>
#include <iostream>
#include <sys/time.h>

using namespace std;

int doSub(int bits, int num1, int num2){
	int qubits = bits;

	auto numberCirc = NumberAsX(qubits, num1);
	auto qft = QFT(qubits, 1, 0, qubits-1);
	auto subf = SubF(qubits, 1, 0, num2, qubits - 1);
	auto rqft = RQFT(qubits, 1, 0, qubits-1);

	vector <string> circ;

	circ.push_back(numberCirc);
	circ.insert(circ.end(), qft.begin(), qft.end());
	circ.push_back(subf);
	circ.insert(circ.end(), rqft.begin(), rqft.end());

	DGM dgm;
	dgm.qubits = qubits;
	dgm.exec_type = t_CPU;

	dgm.n_threads = 1;
	dgm.cpu_region = 14;
	dgm.cpu_coales = 11;
	
	dgm.multi_gpu = 1;
	dgm.gpu_region = 8;
	dgm.gpu_coales = 4;
	dgm.tam_block = 64;
	dgm.rept = 2;

	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	dgm.setFunction(circ, 1);
	dgm.execute(1);

	int res = 0;
	for (int i = 0; i < qubits; i++){
		res = (res << 1) | dgm.measure(i);
	}

	dgm.freeMemory();

	return res;
}

int doAdd(int bits, int num1, int num2){
	int qubits = bits;

	auto numberCirc = NumberAsX(qubits, num1);
	auto qft = QFT(qubits, 1, 0, qubits-1);
	auto add = AddF(qubits, 1, 0, num2, qubits - 1);
	auto rqft = RQFT(qubits, 1, 0, qubits-1);

	vector <string> circ;

	circ.push_back(numberCirc);
	circ.insert(circ.end(), qft.begin(), qft.end());
	circ.push_back(add);
	circ.insert(circ.end(), rqft.begin(), rqft.end());

	DGM dgm;
	dgm.qubits = qubits;
	dgm.exec_type = t_CPU;

	dgm.n_threads = 1;
	dgm.cpu_region = 14;
	dgm.cpu_coales = 11;
	
	dgm.multi_gpu = 1;
	dgm.gpu_region = 8;
	dgm.gpu_coales = 4;
	dgm.tam_block = 64;
	dgm.rept = 2;

	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	dgm.setFunction(circ, 1);
	dgm.execute(1);

	int res = 0;
	for (int i = 0; i < qubits; i++){
		res = (res << 1) | dgm.measure(i);
	}

	dgm.freeMemory();

	return res;
}

int doAddMod(int bits, int num1, int num2){
	int qubits = bits + 2; // overflow qubits and ctrl overflow qubit
	int N = 1 << bits;

	auto numberCirc = "ID," + NumberAsX(bits, num1) + ",ID";

	auto qft = QFT(qubits, 1, 0, qubits-2);
	auto rqft = RQFT(qubits, 1, 0, qubits-2);
	auto addMod = AddMod(qubits, 1, 0, qubits-1, bits, num2, N);

	vector <string> circ;

	circ.push_back(numberCirc);
	circ.insert(circ.end(), qft.begin(), qft.end());
	circ.insert(circ.end(), addMod.begin(), addMod.end());
	circ.insert(circ.end(), rqft.begin(), rqft.end());

	DGM dgm;
	dgm.qubits = qubits;
	dgm.exec_type = t_CPU;

	dgm.n_threads = 1;
	dgm.cpu_region = 14;
	dgm.cpu_coales = 11;
	
	dgm.multi_gpu = 1;
	dgm.gpu_region = 8;
	dgm.gpu_coales = 4;
	dgm.tam_block = 64;
	dgm.rept = 2;

	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	dgm.setFunction(circ, 1);
	dgm.execute(1);

	int res = 0;
	for (int i = 0; i < qubits - 1; i++){
		res = (res << 1) | dgm.measure(i);
	}

	dgm.freeMemory();

	return res;
}

int doC2MultMod(int bits, int b, int a, int x){
	int width = bits;
	int qubits = 2*bits + 3;

	int ctrl = 0;
	int regX = 1;
	int regB = width + 2;
	int over = width + 1;
	int over_bool = 2*width + 2;
	int N = 1 << width;

	int numAsX = (x << (width + 1)) + b;

	auto entranceData = "X," + NumberAsX(width*2 + 1, numAsX) + ",ID";

	auto cMultMod = CMultMod2(qubits, ctrl, regX, regB, over, over_bool, width, a, N);

	vector <string> circ;

	circ.push_back(entranceData);
	circ.insert(circ.end(), cMultMod.begin(), cMultMod.end());

	DGM dgm;
	dgm.qubits = qubits;
	dgm.exec_type = t_CPU;
 
	dgm.n_threads = 1;
	dgm.cpu_region = 14;
	dgm.cpu_coales = 11;
	
	dgm.multi_gpu = 1;
	dgm.gpu_region = 8;
	dgm.gpu_coales = 4;
	dgm.tam_block = 64;
	dgm.rept = 2;

	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	dgm.setFunction(circ, 1);
	dgm.execute(1);

	int res = 0;
	for (int i = regB - 1; i < (regB + width); i++){
		res = (res << 1) | dgm.measure(i);
	}

	dgm.freeMemory();

	return res;
}

int doCU(int bits, int N, int a, int x){
	int width = bits;
	int qubits = 2*bits + 3;

	int ctrl = 0;
	int regX = 1;
	int over = regX + width;
	int regB = over + 1;
	int over_bool = regB + width;
	int inv_a = mul_inv(a, N);

	vector <string> circ;

	int numAsX = (x << (width + 1));

	cout << "bits: " << bits << endl;
	cout << "N: " << N << endl;
	cout << "a: " << a << endl;
	cout << "inv_a: " << inv_a << endl;
	cout << "qubits: " << qubits << endl;
	cout << "regX: " << regX << endl;
	cout << "over: " << over << endl;
	cout << "regB: " << regB << endl;
	cout << "over_bool: " << over_bool << endl;
	
	auto entranceData = "H," + NumberAsX(width*2 + 1, numAsX) + ",ID";

	auto cMultMod = CMultMod2(qubits, ctrl, regX, regB, over, over_bool, width, a, N);
	auto swapR = CSwapR(qubits, ctrl, regX, regB, width);
	auto rcMultMod = CRMultMod2(qubits, ctrl, regX, regB, over, over_bool, width, inv_a, N);

	cout << "entranceData" << endl;
	cout << entranceData << endl;

	cout << "\ncMultMod" << endl;
	for (int i = 0; i < cMultMod.size(); i++){
		cout << cMultMod[i] << endl;
	}

	cout << "\nswapR" << endl;
	for (int i = 0; i < swapR.size(); i++){
		cout << swapR[i] << endl;
	}

	cout << "\nrcMultMod" << endl;
	for (int i = 0; i < rcMultMod.size(); i++){
		cout << rcMultMod[i] << endl;
	}

	circ.push_back(entranceData);
	circ.insert(circ.end(), cMultMod.begin(), cMultMod.end());
	circ.insert(circ.end(), swapR.begin(), swapR.end());
	circ.insert(circ.end(), rcMultMod.begin(), rcMultMod.end());

	DGM dgm;
	dgm.qubits = qubits;
	dgm.exec_type = t_CPU;
 
	dgm.n_threads = 1;
	dgm.cpu_region = 14;
	dgm.cpu_coales = 11;
	
	dgm.multi_gpu = 1;
	dgm.gpu_region = 8;
	dgm.gpu_coales = 4;
	dgm.tam_block = 64;
	dgm.rept = 2;

	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	dgm.executeFunction(entranceData);
	cout << "Entrance" << endl;
	for (int i = 0; i < qubits; i++){
		dgm.printProbability(i);
	}
	cout << "###############" << endl;

	dgm.executeFunction(cMultMod);
	cout << "cMultMod" << endl;
	for (int i = 0; i < qubits; i++){
		dgm.printProbability(i);
	}
	cout << "###############" << endl;

	dgm.executeFunction(swapR);
	cout << "swapR" << endl;
	for (int i = 0; i < qubits; i++){
		dgm.printProbability(i);
	}

	dgm.executeFunction(rcMultMod);
	cout << "rcMultMod" << endl;
	for (int i = 0; i < qubits; i++){
		dgm.printProbability(i);
	}

	printMem(dgm.state, dgm.qubits);

	dgm.colapse(ctrl, 1);
	cout << "After Colapse" << endl;
	for (int i = 0; i < qubits; i++){
		dgm.printProbability(i);
	}

	printMem(dgm.state, dgm.qubits);

	int res = 0;
	for (int i = 0; i < qubits; i++){
		res = (res << 1) | dgm.measure(i);
	}

	dgm.freeMemory();

	cout << "Result: " << res << endl;

	return res;
}

int doShor(int N, int a){
	DGM dgm;
	dgm.exec_type = t_CPU;
 
	dgm.n_threads = 1;
	dgm.cpu_region = 14;
	dgm.cpu_coales = 11;
	
	dgm.multi_gpu = 1;
	dgm.gpu_region = 8;
	dgm.gpu_coales = 4;
	dgm.tam_block = 64;
	dgm.rept = 2;

	int bits = log2(N) + 1;
	int qubits = 2*bits + 3;

	dgm.qubits = qubits;
	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	int width = bits;
	int inv_a = mul_inv(a, N);

	int ctrl = 0;
	int regX = 1;
	int over = regX + width;
	int regB = over + 1;
	int over_bool = regB + width;

	cout << "N: " << N << endl;
	cout << "a: " << a << endl;
	cout << "bits: " << bits << endl;
	cout << "qubits: " << qubits << endl;
	cout << "width: " << width << endl;
	cout << "regX: " << regX << endl;
	cout << "over: " << over << endl;
	cout << "regB: " << regB << endl;
	cout << "over_bool: " << over_bool << endl;

	auto entranceData = Pauli_X(qubits, regX + width - 1);

	cout << "\nentranceData" << endl;
	cout << entranceData << endl;

	dgm.executeFunction(entranceData);

	cout << "Entrance probability: " << endl;
	for (int i = 0; i < qubits; i++){
		dgm.printProbability(i);
	}

	string H0 = Hadamard(qubits, ctrl);
	cout << "H0:\n" << H0 << endl;

	vector <string> circ;

	long res = 0;

	printMem(dgm.state, dgm.qubits);

	int L = 2*width;
	for (int i = 0; i < 1; i++){
		int mod_a = modular_pow(a, pow(2,i), N);
		int mod_inv_a = modular_pow(inv_a, pow(2,i), N);

		cout << "\nL: " << i << endl;
		cout << "mod_a: " << mod_a << endl;
		cout << "mod_inv_a: " << mod_inv_a << endl;

		circ.clear();

		auto cMultMod = CMultMod2(qubits, ctrl, regX, regB, over, over_bool, width, a, N);
		auto swapR = CSwapR(qubits, ctrl, regX, regB, width);
		auto rcMultMod = CRMultMod2(qubits, ctrl, regX, regB, over, over_bool, width, inv_a, N);

		circ.push_back(H0);
		dgm.executeFunction(H0);
		cout << "Hadamard 0: " << endl;
		for (int i = 0; i < qubits; i++){
			dgm.printProbability(i);
		}

		dgm.executeFunction(cMultMod);
		cout << "cMultMod" << endl;
		for (int i = 0; i < qubits; i++){
			dgm.printProbability(i);
		}

		dgm.executeFunction(swapR);
		cout << "swapR" << endl;
		for (int i = 0; i < qubits; i++){
			dgm.printProbability(i);
		}

		dgm.executeFunction(rcMultMod);
		cout << "rcMultMod" << endl;
		for (int i = 0; i < qubits; i++){
			dgm.printProbability(i);
		}

		// add rot
		if (res) {
			auto rot = genRot2(qubits, ctrl, res);
			dgm.executeFunction(rot);
			cout << "Rot" << endl;
			for (int i = 0; i < qubits; i++){
				dgm.printProbability(i);
			}
		}

		printMem(dgm.state, dgm.qubits);

		dgm.executeFunction(H0);
		cout << "Hadamard 1: " << endl;
		for (int i = 0; i < qubits; i++){
			dgm.printProbability(i);
		}

		//cout << "\nL Circ" << endl;
		//for (int i = 0; i < circ.size(); i++){
		//	cout << circ[i] << endl;
		//}

		//dgm.executeFunction(circ);

		int m = dgm.measure(ctrl);

		cout << "After Measurement: " << endl;
		for (int i = 0; i < qubits; i++){
			dgm.printProbability(i);
		}

		printMem(dgm.state, dgm.qubits);

		cout << "###################" << endl;

		res = res << 1;

		if (m)
			res |= 1;
	}

	dgm.freeMemory();

	return res;
}

int doFullShor(int N, int a){
	DGM dgm;
	dgm.exec_type = t_CPU;
 
	dgm.n_threads = 1;
	dgm.cpu_region = 14;
	dgm.cpu_coales = 11;
	
	dgm.multi_gpu = 1;
	dgm.gpu_region = 8;
	dgm.gpu_coales = 4;
	dgm.tam_block = 64;
	dgm.rept = 2;

	int bits = log2(N) + 1;
	int qubits = 4*bits + 2;

	dgm.qubits = qubits;
	dgm.allocateMemory();
	dgm.setMemoryValue(0);

	int width = bits;
	int inv_a = mul_inv(a, N);

	int L = 2*width;
	int regX = L;
	int over = regX + width;
	int regB = over + 1;
	int over_bool = regB + width;

	//cout << "N: " << N << endl;
	//cout << "a: " << a << endl;
	//cout << "bits: " << bits << endl;
	//cout << "qubits: " << qubits << endl;
	//cout << "width: " << width << endl;
	//cout << "regX: " << regX << endl;
	//cout << "over: " << over << endl;
	//cout << "regB: " << regB << endl;
	//cout << "over_bool: " << over_bool << endl;

	string HL = Hadamard(qubits, 0, L);
	//cout << "\nHL" << endl;
	//cout << HL << endl;

	dgm.executeFunction(HL);

	auto entranceX = Pauli_X(qubits, regX + width - 1);

	//cout << "\nentranceX" << endl;
	//cout << entranceX << endl;

	dgm.executeFunction(entranceX);

	//cout << "Entrance probability: " << endl;
	//for (int i = 0; i < qubits; i++){
	//	dgm.printProbability(i);
	//}

	vector <string> circ;

	for (int i = 0; i < L; i++){
		int mod_a = modular_pow(a, pow(2,i), N);
		int mod_inv_a = modular_pow(inv_a, pow(2,i), N);

		//cout << "\nL: " << i << endl;
		//cout << "mod_a: " << mod_a << endl;
		//cout << "mod_inv_a: " << mod_inv_a << endl;

		circ.clear();

		auto cMultMod = CMultMod2(qubits, i, regX, regB, over, over_bool, width, a, N);
		auto swapR = CSwapR(qubits, i, regX, regB, width);
		auto rcMultMod = CRMultMod2(qubits, i, regX, regB, over, over_bool, width, inv_a, N);

		circ.insert(circ.end(), cMultMod.begin(), cMultMod.end());
		circ.insert(circ.end(), swapR.begin(), swapR.end());
		circ.insert(circ.end(), rcMultMod.begin(), rcMultMod.end());

		dgm.executeFunction(circ);
	}

	auto rqft = RQFT(qubits, 1, 0, L-1);
	//cout << "RQFT" << endl;
	//for (int i = 0; i < rqft.size(); i++){
	//	cout << rqft[i] << endl;
	//}

	dgm.executeFunction(rqft);

	long res = 0;
	for (int i = 0; i < L; i++){
		int m = dgm.measure(i);
		res = (res << 1) | m;
	}

	dgm.freeMemory();

	return res;
}

int main(int argc, char **argv){
	int n_threads = 1, cpu_region = 14, cpu_coalesc = 11, multi_gpu = 1, gpu_region = 8, gpu_coalesc = 4, tam_block = 64, rept = 2;

	/*
	int bits = atoi(argv[1]);
	int N = atoi(argv[2]);
	int a = atoi(argv[3]);
	int x = atoi(argv[4]);

	int res = doCU(bits, N, a, x);
	*/

	int N = atoi(argv[1]);
	int a = atoi(argv[2]);

	vector<long> results;
	for (int i = 0; i < 30; i++){
		int res = doFullShor(N, a);

		results.push_back(res);
	}

	for (int i = 0; i < results.size(); i++){
		cout << "Result: " << results[i] << " - " << revert_bits(results[i], 8) << endl;
	}

	return 0;
}