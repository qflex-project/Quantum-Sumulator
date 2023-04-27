#include "lib_grover.h"
#include "dgm.h"
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <omp.h>

using namespace std;

int main(int argc, char **argv){

	struct timeval timev, tvBegin, tvEnd;
	float t;
	int execType = t_CPU, n_threads = 1, cpu_region = 14, cpu_coalesc = 11, multi_gpu = 1, gpu_region = 8, gpu_coalesc = 4, tam_block = 64, rept = 2, seed = 0;

	if (argc < 2) {
		cout << "You need to define the execution parameters" << endl;
		return 0;
	}

	int qubits = atoi(argv[1]);


	if (argc > 2) {
		execType = atoi(argv[2]);
	}

  	if (argc > 3) {
		seed = atoi(argv[3]);
	}
	srand(seed);

	if (argc > 4) {
		cpu_region = atoi(argv[4]);
	}

	if (argc > 5) {
		cpu_coalesc = atoi(argv[5]);
	}

	int value = 10;
	if (argc > 6) {
		value = atoi(argv[6]);
	}

	if (execType < t_CPU || execType > t_HYBRID){
		cout << "Invalid execution type: " << execType << endl; 
		return 0;
	}

	if (execType == t_PAR_CPU || execType == t_HYBRID) {
		n_threads = omp_get_max_threads();
	} else if (execType == t_GPU) {
		if (argc > 7) {
			multi_gpu = atoi(argv[7]);
		}
	}
 
	cout << "Executing Grover: " << qubits << " qubits" << endl;

	gettimeofday(&tvBegin, NULL);
	float result = Grover(qubits, value, execType, n_threads, cpu_region, cpu_coalesc, multi_gpu, gpu_region, gpu_coalesc, tam_block, rept);
 	gettimeofday(&tvEnd, NULL);
	timeval_subtract(&timev, &tvEnd, &tvBegin);
	t = timev.tv_sec + (timev.tv_usec / 1000000.0);
	
	long reg_size_per_thread = (1 << (qubits - cpu_region));
	std::cout << "reg_size_per_thread: " << reg_size_per_thread << std::endl;

	cout << "Time: " << t << endl;

	if (result != -1) {
		cout << "Result found successfully: " << result << endl;
	}
	else {
		cout << "Failed to find the value" << endl;
	}
	
}
