#include "lib_shor.h"
#include "dgm.h"
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <map>
#include <omp.h>

using namespace std;

int main(int argc, char** argv){
	map <int, int> qubitsMap;
	qubitsMap[15] = 57;	
	qubitsMap[17] = 119;
	qubitsMap[19] = 253;
	qubitsMap[21] = 485;
	qubitsMap[23] = 1017;
	qubitsMap[25] = 2045;
	qubitsMap[27] = 2863;

	struct timeval timev, tvBegin, tvEnd;
	float t;
	vector <float> amostras;

	int execType = t_CPU, n_threads = 1, cpu_region = 14, cpu_coalesc = 11, multi_gpu = 1, gpu_region = 8, gpu_coalesc = 4, tam_block = 64, rept = 2, seed = 0;

	if (argc < 2){
		cout << "You need to define the execution parameters" << endl;
		return 0;
	}

	int qubits = atoi(argv[1]);
	if (qubitsMap.count(qubits) == 0){
		cout << "The amount of qubits does not map to a valid number to be factored: " << qubits << endl; 
		return 0;
	}
 
	if (argc > 2) {
		execType = atoi(argv[2]);
	}

	if (argc > 3) {
		seed = atoi(argv[3]);
	}
	srand(seed);

	if (execType < t_CPU || execType > t_HYBRID){
		cout << "Invalid execution type: " << execType << endl; 
		return 0;
	}

	if (execType == t_PAR_CPU || execType == t_HYBRID) {
		n_threads = omp_get_max_threads();
	} else if (execType == t_GPU) {
		if (argc > 4) {
			multi_gpu = atoi(argv[4]);
		}
	}

	vector<int> factors;

	cout << "Executing Shor: " << qubits << " qubits" << endl;

	gettimeofday(&tvBegin, NULL);
	factors = Shor(qubitsMap[qubits], execType, n_threads, cpu_region, cpu_coalesc, multi_gpu, gpu_region, gpu_coalesc, tam_block, rept);
	gettimeofday(&tvEnd, NULL);
	timeval_subtract(&timev, &tvEnd, &tvBegin);
	t = timev.tv_sec + (timev.tv_usec / 1000000.0);
	
	cout << "Time: " << t << endl;

	if (factors.size() == 2) {
		cout << "Found factors: " << factors[0] << " -- " << factors[1] << endl;
	}
	else {
		cout << "Failed to find factors" << endl;
	}

	return 0;
}
