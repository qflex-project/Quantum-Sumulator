#ifndef _DGM_H_
#define _DGM_H_

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include <complex>
#include <sys/time.h>
#include <iostream>
#include "common.h"
#include "gates.h"

#include <mutex>
#include <condition_variable>

#define ccomplex _Complex

using namespace std;

std::complex <float>* GenericExecute(std::complex <float> *state, string function, int qubits, int type, int threads, int factor);
std::complex <float>*  GenericExecute(std::complex <float> *state, vector<string> function, int qubits, int type, int threads, int factor);

bool setDevice(int num = 0);

std::complex <float>* GpuExecutionWrapper(std::complex <float>* r_memory, PT **pts, int qubits, int multi_gpu, int coalesc, int qbs_region, int tam_block, int rept, int num_it);
std::complex <float>* GpuExecution(std::complex <float>* r_memory, std::complex <float>* w_memory, PT **pts, int qubits, float *total_time, long MAX_PT, long MAX_QB, int it);
std::complex <float>* GpuExecution2(std::complex <float>* r_memory, PT **pts, int pts_size, int qubits, long MAX_PT, int it);
std::complex <float>* GpuExecution3(std::complex <float>* r_memory, std::complex <float>* w_memory, int sub_size, int shift_write, PT *pt, int qubits, long MAX_PT, long MAX_QB, int it);
bool ProjectState(std::complex <float>* state, int qubits, int region_size, long reg_id, long reg_mask, int multi_gpu);
bool GetState(std::complex <float>* state, int qubits, int region_size, long reg_id, long reg_mask, int multi_gpu);

void PCpuExecution1(std::complex <float> *state, PT **pts, int qubits, long n_threads, int coales, int region, int it);
void PCpuExecution1_0(std::complex <float> *state, PT **pts, int qubits, int start, int end, int pos_count, int reg_id, int reg_mask);

std::string getBinaryString(long num, int n, bool includeNum = false);

inline long LINE (long pos, long shift){
	return ((pos >> shift) & 1) * 2;
}
inline long BASE (long pos, long shift){
	return pos & (~(1 << shift));
}

enum {
	t_CPU,
	t_PAR_CPU,
	t_GPU,
	t_HYBRID,
	t_SPEC
};

class Projection {
	public:
		long qubits = 0;
		long cur_id = 0;
		long next_id = 0;
		long region_mask = 0;
		long region_size = 0;
		long count = 0;
		long total = 0;
		long operator_start = 0;
		long operators_count = 0;

		long inc_mask = 0;
		long parent_proj_id = 0;

		std::mutex mutex_;

		Projection(){};
		void setData(Projection *parent_proj, PT **pts, long start, long end, long qubits, long coales, long region_size, bool include_main_diag = false);
		long getNextProjectionId();
		void printInfo();
};

class Group{
public:
	vector <string> ops;
	vector <long> pos_ops;
	vector <bool> ctrl;
	vector <long> pos_ctrl;

	Group(){};
	bool isAfected(int pos, int afect);
};

class DGM{
public:
	long total_op, dense, main_diag, sec_diag, c_dense, c_main_diag, c_sec_diag; //counters

	vector <string> diag;
	long MAX_QB, MAX_PT, qb_afected;

	long factor;

	int exec_type;
 
 	long n_threads;
	long cpu_coales;
	long cpu_region;

	int multi_gpu;
	long gpu_coales;
	long gpu_region;
	 
	int tam_block;
	int rept;
	
	vector <PT*> vec_pts;
	PT** pts;
	long qubits;

	float measure_value;

	float elapsed_time;
	struct timeval timev;

	std::complex <float> *state;

	DGM();
	~DGM();

	bool en_print;

	void printPTs();
	void erase();
	void setExecType(int type);

	void setCpuStructure(long cpu_region, long cpu_coales);
	void setGpuStructure(long gpu_coales, long gpu_region, int rept = 1);

	void allocateMemory();
	void setMemory(std::complex <float> *mem);
	void freeMemory();
	void setMemoryValue(int pos);
	void setSuperposition();
	
	int measure(int q_pos);
	map <long, float> measure(vector<int> q_pos);
	void colapse(int q_pos, int value);

	void setFunction(string function, int it = 1, bool er = true);
	void setFunction(vector<string> steps, int it = 1, bool er = true);
	map <long, Group> genGroups(string step);
	void genPTs(map<long, Group> &gps, vector <PT*> &step_pts);
	void genMatrix(std::complex <float>* matrix, vector<std::complex <float>*> &matrices, long tam, long current, long line, long column, std::complex <float> cmplx);

	void CountOps(int it = 1);

	void executeFunction(string function, int it = 1);
	void executeFunction(vector<string> steps, int it = 1);
	std::complex <float>* execute(int it);

	void HybridExecution(PT **pts);
	void HybridExecution2(PT **pts);

	void CpuExecution1(int it);
	void CpuExecution1_1(PT *pt, long mem_size);
	void CpuExecution1_2(PT *pt, long mem_size);
	void CpuExecution1_3(PT *pt, long mem_size);

	//void CpuExecution1(int it);
	void CpuExecution2_1(PT *pt, long mem_size);
	void CpuExecution2_2(PT *pt, long mem_size);
	void CpuExecution2_3(PT *pt, long mem_size);

	void CpuExecution3_1(PT *pt, long mem_size);
	void CpuExecution3_2(PT *pt, long mem_size);
	void CpuExecution3_3(PT *pt, long mem_size);

};

class Barrier {
public:
	explicit Barrier(int count) : count_(count), initial_count_(count) {}

	void arrive_and_wait(int point = 0) {
		std::unique_lock<std::mutex> lock(mutex_);
		int generation = generation_;

		if (point_ == 0) {
			point_ = point;
		} else if (point_ != point) {
			std::cout << std::string("ERROR: barrier being called from differ points: " + std::to_string(point) + " and " + std::to_string(point_) + "\n");
		}

		if (--count_ == 0) {
			generation_++;
			count_ = initial_count_;
			point_ = 0;
			cv_.notify_all();
		} else {
			cv_.wait(lock, [this, generation] { return generation != generation_; });
		}
	}

private:
	std::mutex mutex_;
	std::condition_variable cv_;
	int count_;
	int initial_count_;
	int generation_ = 0;
	int point_ = 0;
};

#endif
