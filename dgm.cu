#include <iostream>
#include "dgm.h"
#include <omp.h>
#include <unistd.h>
#include <cstdio>
#include <iterator>
#include <bitset>
#include <random>

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = ",")
{
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////

std::complex <float>* GenericExecute(std::complex <float> *state, string function, int qubits, int type, int threads, int factor = 0){
	DGM dgm;
	dgm.exec_type = type;
	dgm.n_threads = threads;
	dgm.qubits = qubits;
	dgm.factor = factor;

	dgm.setMemory(state);

	dgm.executeFunction(function);

	state = dgm.state;

	dgm.state = NULL;

	return state;
}

std::complex <float>* GenericExecute(std::complex <float> *state, vector<string> function, int qubits, int type, int threads, int factor = 0){
	DGM dgm;
	dgm.exec_type = type;
	dgm.n_threads = threads;
	dgm.qubits = qubits;
	dgm.factor = factor;
	dgm.setMemory(state);

	dgm.executeFunction(function);

	dgm.state = NULL;

	return state;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DGM::DGM(){
	MAX_QB = QB_LIMIT;
	MAX_PT = PT_TAM;

	pts = NULL;
	state = NULL;
	en_print = false;
	exec_type = t_CPU;
	factor = 1;
	multi_gpu = 1;
}

DGM::~DGM(){erase();}

void DGM::setExecType(int type){
	exec_type = type;
}

void DGM::printPTs(){
	for (int i = 0; i < vec_pts.size() -1; i++){
		vec_pts[i]->print();
	}
}

void DGM::erase(){
	if (!pts) return;

	long i = 0;
	while (pts[i] != NULL){
		pts[i]->destructor();
		free(pts[i]);
		i++;
	}

	vec_pts.clear();
	pts = NULL;
}

void DGM::allocateMemory(){
	state = (std::complex <float>*) calloc(pow(2, qubits), sizeof(std::complex <float>));
}

void DGM::setMemory(std::complex <float>* mem){
	freeMemory();
	state = mem;
}

void DGM::freeMemory(){
	if (state) free(state);
	state = NULL;
}

void DGM::setMemoryValue(int pos){
	state[pos] = 1;
}

void DGM::setSuperposition(){
	auto value = sqrt(1.0/(pow(2, qubits)));
	for (long i = 0; i < pow(2, qubits); i++) state[i] = value;
}

int DGM::measure(int q_pos){
	long size = pow(2.0, qubits);

	long shift = (qubits - 1 - q_pos);
	long mask = 1 << shift;

	int count_one, count_zero, num_pb;
	float zero, one, norm_factor, r;
	one = zero = 0;

	//#pragma omp for;
	for (long i = 0; i < size; i++){
		if (i & mask)
			one += norm(state[i]);
		else
			zero += norm(state[i]);
	}

	long m;

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(0.05, 0.95);

	r = dis(gen);

	if (r > zero){
		norm_factor = sqrt(one);
		m = 1;
	}
	else{
		norm_factor = sqrt(zero);
		m = 0;
	}

	//#pragma omp for
	for (long i = 0; i < size/2; i++){
		long pos0 = (i * 2) - (i & (mask-1));
		long pos1 = pos0 | mask;

		//std::cout << "Pos0 " << pos0 << ": " << real(state[pos0]) << " + " << imag(state[pos0]) << "i" << endl;
		//std::cout << "Pos1 " << pos1 << ": " << real(state[pos1]) << " + " << imag(state[pos1]) << "i" << endl;

		if (m) {
			state[pos0]	= state[pos1]/norm_factor;
		}
		else {
			state[pos0]	= state[pos0]/norm_factor;
		}

		state[pos1] = 0.0;

		//std::cout << "Pos0 " << pos0 << ": " << real(state[pos0]) << " + " << imag(state[pos0]) << "i" << endl;
		//std::cout << "Pos1 " << pos1 << ": " << real(state[pos1]) << " + " << imag(state[pos1]) << "i" << endl;

		//cout << "###" << endl;
	}

	std::cout << q_pos << ": " << m << " (" << zero << " , " << one << ") " << r << " " << mask << " " << shift << " (" << sqrt(zero) << " , " << sqrt(one) << ") " << std::endl;

	return m;
}

void DGM::colapse(int q_pos, int value){
	long size = pow(2.0, qubits);

	long shift = (qubits - 1 - q_pos);
	long mask = 1 << shift;

	float norm_factor = 0;
	for (long i = 0; i < size; i++){
		if (((i >> shift) & 1) == value) {
			norm_factor += norm(state[i]);
		}
	}

	norm_factor = sqrt(norm_factor);
	int m = value;

	for (long i = 0; i < size/2; i++){
		long pos0 = (i * 2) - (i & (mask-1));
		long pos1 = pos0 | mask;

		if (m) {
			state[pos0]	= state[pos1]/norm_factor;
		}
		else {
			state[pos0]	= state[pos0]/norm_factor;
		}

		state[pos1] = 0.0;
	}

	std::cout << "colapse " << q_pos << " - " << value << " - " << norm_factor << std::endl;
}

void DGM::printProbability(int q_pos){
	long size = pow(2.0, qubits);

	long shift = (qubits - 1 - q_pos);
	long mask = 1 << shift;

	int count_one, count_zero, num_pb;
	float zero, one, norm_factor, r;
	one = zero = 0;

	//#pragma omp for;
	for (long i = 0; i < size; i++){
		if (i & mask)
			one += norm(state[i]);
		else
			zero += norm(state[i]);
	}

	if (std::abs(zero) < 1e-5) zero = 0;
	if (std::abs(one) < 1e-5) one = 0;

	std::cout << q_pos << " -- (" << zero << " , " << one << ")" << std::endl;
}


map <long, float> DGM::measure(vector<int> q_pos){
	long mask = 0;

	for (int i =0; i < q_pos.size(); i++) mask = mask | (1<<(qubits - 1 - q_pos[i]));

	map <long, float> m;

	long size = pow(2.0, qubits);

	for (long i =0; i < size; i++) m[i&mask] += pow(state[i].real(), 2.0) + pow(state[i].imag(), 2.0);

	return m;
}

void DGM::setFunction(string function, int it, bool er){
	vector <string> steps;

	Tokenize(function, steps, ";");

	setFunction(steps, it, er);
}

void DGM::setFunction(vector <string> steps, int it, bool er){
	if (er) erase();
	else vec_pts.pop_back();


	vector <PT*> step_pts, vec_tmp;
	map<long, Group> gps;

	for (long j = 0; j< it; j++)
	for (long i = 0; i < steps.size(); i++){
		gps = genGroups(steps[i]);
		genPTs(gps, step_pts);

		if (i%2)
			sort(step_pts.begin(), step_pts.end(), increasing);
		else
			sort(step_pts.begin(), step_pts.end(), decreasing);

		vec_pts.insert(vec_pts.end(), step_pts.begin(), step_pts.end());
	}

	vec_pts.push_back(NULL);

	pts = &vec_pts[0];
}

map <long, Group> DGM::genGroups(string step){
	vector <string> ops;
	Tokenize(step, ops); //separa os operadores usando "," como delimitador
	qubits = ops.size();

	size_t found_c, found_t, p;
	string str;
	long pos, ctrl_value, ctrl_num;
	
	map<long, Group> gps;

	char * pEnd;
	pos = 0;
	vector<string>::iterator it;
	for (it = ops.begin() ; it != ops.end(); ++it){ //percorre os operadores
		str = *it;
		//cout << str << endl;
		found_c = str.find("Control"); //tamanho 7
		found_t = str.find("Target");  //tamanho 6
		p = str.find("(") + 1;

		if (found_c != string::npos){ //Controle
			ctrl_num = strtol(str.c_str()+7, &pEnd, 10);
			ctrl_value = strtol(str.c_str()+p, &pEnd, 10);

			gps[ctrl_num].ctrl.push_back(ctrl_value); //adicona o valor do controle
			gps[ctrl_num].pos_ctrl.push_back(pos);  //e a sua posição ao map relacionado ao controle
		}
		else if(found_t != string::npos){ //Target
			ctrl_num = strtol(str.c_str()+6, &pEnd, 10);
			str = str.substr(p, str.size()-p-1);

			gps[ctrl_num].ops.push_back(str);     //adicona o operador
			gps[ctrl_num].pos_ops.push_back(pos); //e a sua posição ao map relacionado ao target
		}
		else{ //operador normal
			if (str != "ID"){ //se for ID ignora
				gps[0].ops.push_back(str);     //adiciona o operador
				gps[0].pos_ops.push_back(pos); //e a sua posição ao map '0'
			}
		}
		pos++;
	}
	
	return gps;
}

void DGM::genPTs(map<long, Group> &gps, vector <PT*> &step_pts){
	step_pts.clear();
	Gates gates;

	map<long,Group>::iterator it;	
	Group gp;
	PT* pt;
	long ctrl_mask, ctrl_value, ctrl_count;
	long size;
	
	for (it = gps.begin(); it != gps.end(); ++it){ //percorre os grupos
		gp = it->second;
		size = gp.ops.size();
		
		ctrl_count = gp.ctrl.size();
		ctrl_value = ctrl_mask = 0;

		for (long i = 0; i < ctrl_count; i++){ //gera a mascara e o valor do controle (em binario)
			gp.pos_ctrl[i] =  qubits - gp.pos_ctrl[i] - 1;
			ctrl_mask += (1 << gp.pos_ctrl[i]);
			if (gp.ctrl[i]) ctrl_value += (1 << gp.pos_ctrl[i]);
		}

		for (int p = 0; p < size; p++){
			
			pt = (PT*) malloc(sizeof(PT));
			pt->affected = false;

			pt->qubits = 1;
			pt->start = qubits - gp.pos_ops[p];
			pt->end = pt->start - 1;
			pt->mat_size = 2;
			
			pt->matrix = gates.getMatrix(gp.ops[p]);

			pt->ctrl_value = ctrl_value;
			pt->ctrl_mask = ctrl_mask;
			pt->ctrl_count = ctrl_count;

			if (ctrl_count){
				pt->ctrl_pos = (long*)malloc(sizeof(long) * ctrl_count);
				copy(gp.pos_ctrl.begin(), gp.pos_ctrl.end(), pt->ctrl_pos);
			}

			step_pts.push_back(pt);
		}
	}
}

void DGM::genMatrix(std::complex <float>* matrix, vector<std::complex <float>*> &matrices, long tam, long current, long line, long column, std::complex <float> cmplx){
	if (cmplx == COMPLEX_ZERO) return;

	if (current == tam){ //percorreu até a ultima matriz
		matrix[line*(1<<tam) + column] = cmplx;
		return;
	}

	for (long l = 0; l < 2; l++)
		for (long c = 0; c < 2; c++)
			genMatrix(matrix, matrices, tam, current+1, (line<<1)|l, (column<<1)|c, cmplx * matrices[current][l*2+c]);
}


void DGM::executeFunction(vector <string> function, int it){
	setFunction(function);
	execute(it);
}

void DGM::executeFunction(string function, int it){
	if (function == "") return;

	setFunction(function);
	execute(it);
}


std::complex <float>* DGM::execute(int it){
	std::complex <float>* result = state;

	switch (exec_type){
		case t_CPU:
			CpuExecution1(it);
			break;
		case t_PAR_CPU:
			PCpuExecution1(state, pts, qubits, n_threads, cpu_coales, cpu_region, it);
			break;
		#ifndef ONLY_CPU
		case t_GPU:
			result = GpuExecutionWrapper(state, pts, qubits, gpu_coales, gpu_region, multi_gpu, tam_block, rept, it);
			break;
		case t_HYBRID:
			HybridExecution2(pts);
			break;
		#endif
		default:
			cout << "Erro exec type" << endl;
			exit(1);
	}

	return result;
}


void DGM::CountOps(int it){
	dense = main_diag = sec_diag = c_dense = c_main_diag = c_sec_diag = 0;

	for (int i =0; pts[i]!=NULL; i++){
		long mt = pts[i]->matrixType();
		switch (mt){
			case DENSE:
				(pts[i]->ctrl_mask) ? c_dense++ : dense++;
				break;
			case DIAG_PRI:
				(pts[i]->ctrl_mask) ? c_main_diag++ : main_diag++;
				break;
			case DIAG_SEC:
				(pts[i]->ctrl_mask) ? c_sec_diag++ : sec_diag++;
				break;
			default:
				cout << "Error on operator type" << endl;
				exit(1);
		}
	}

	dense *= it;
	c_dense *= it;
	main_diag *= it;
	c_main_diag *= it;
	sec_diag *= it;
	c_sec_diag *= it;

	total_op = dense + c_dense + main_diag + c_main_diag + sec_diag + c_sec_diag;
}

void DGM::CpuExecution1(int it){
	long mem_size = pow(2.0, qubits);

	for (int x = 0; x < it; x++){
		long i = 0;
		while (pts[i] != NULL){
			long mt = pts[i]->matrixType();

			switch (mt){
				case DENSE:
					CpuExecution1_1(pts[i], mem_size);
					break;
				case DIAG_PRI:
					CpuExecution1_2(pts[i], mem_size);
					break;
				case DIAG_SEC:
					CpuExecution1_3(pts[i], mem_size);
					break;
				default:
					exit(1);
			}
			i++;
		}
	}
}

void DGM::CpuExecution1_1(PT *pt, long mem_size){ //Denso
	long pos0, pos1, shift;
	
	shift = 1 << pt->end;
	
	std::complex <float> tmp;
		
	if (!pt->ctrl_count){ 			//operador não controlado
		mem_size /= 2;
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
			state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];
			state[pos0] = tmp;
		}
	}
	else{					//operador controlado
		long mask = ~(pt->ctrl_mask | shift);
		long inc = (~mask) + 1;

		for (long pos = 0; pos < mem_size; pos = (pos+inc) & mask){
			pos0 = pos | pt->ctrl_value;
			pos1 = pos0 | shift;

			tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
			state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];			
			state[pos0] = tmp;
		}
	}
}

void DGM::CpuExecution1_2(PT *pt, long mem_size){ //Diagonal Principal
	long pos0, shift = pt->end;
		
	if (!pt->ctrl_count)	//operador não controlado
		for (long pos = 0; pos < mem_size; pos++)
			state[pos] = pt->matrix[((pos >> shift) & 1) * 3] * state[pos];
	else{					//operador controlado
		long mask = ~(pt->ctrl_mask);
		long inc = (~mask) + 1;

		for (long pos = 0; pos < mem_size; pos = (pos+inc) & mask){
			pos0 = pos | pt->ctrl_value;

			state[pos0] = pt->matrix[((pos0 >> shift) & 1) * 3] * state[pos0];
		}
	}
}

void DGM::CpuExecution1_3(PT *pt, long mem_size){ //Diagonal Secundária
	long pos0, pos1, shift;
	
	shift = 1 << pt->end;

	std::complex <float> tmp;
		
	if (!pt->ctrl_count){ 	//operador não controlado
		mem_size /= 2;
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			tmp = pt->matrix[1] * state[pos1];
			state[pos1] = pt->matrix[2] * state[pos0];
			state[pos0] = tmp;
		}
	}
	else{					//operador controlado
		long mask = ~(pt->ctrl_mask | shift);
		long inc = (~mask) + 1;
		
		for (long pos = 0; pos < mem_size; pos = (pos+inc) & mask){
			pos0 = pos | pt->ctrl_value;
			pos1 = pos0 | shift;

			tmp = pt->matrix[1] * state[pos1];
			state[pos1] = pt->matrix[2] * state[pos0];
			state[pos0] = tmp;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DGM::CpuExecution2_1(PT *pt, long mem_size){ //Denso
	long pos0, pos1, shift;
	
	shift = 1 << pt->end;
	mem_size /= 2;

	std::complex <float> tmp;
		
	if (!pt->ctrl_count) 			//operador não controlado
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
			state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];
			state[pos0] = tmp;
		}
	else{					//operador controlado
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;
			if ((pos0 & pt->ctrl_mask) == pt->ctrl_value){
				tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
				state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];			
				state[pos0] = tmp;
			}
		}
	}
}

void DGM::CpuExecution2_2(PT *pt, long mem_size){ //Diagonal Principal
	long shift = pt->end;
		
	if (!pt->ctrl_count)	//operador não controlado
		for (long pos = 0; pos < mem_size; pos++)
			state[pos] = pt->matrix[((pos >> shift) & 1) * 3] * state[pos];
	else					//operador controlado
		for (long pos = 0; pos < mem_size; pos++)
			if ((pos & pt->ctrl_mask) == pt->ctrl_value)
				state[pos] = pt->matrix[((pos >> shift) & 1) * 3] * state[pos];

}



void DGM::CpuExecution2_3(PT *pt, long mem_size){ //Diagonal Secundária
	long pos0, pos1, shift;
	
	shift = 1 << pt->end;
	mem_size /= 2;

	std::complex <float> tmp;
		
	if (!pt->ctrl_count) 	//operador não controlado
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			tmp = pt->matrix[1] * state[pos1];
			state[pos1] = pt->matrix[2] * state[pos0];
			state[pos0] = tmp;
		}
	else					//operador controlado
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;
			if ((pos0 & pt->ctrl_mask) == pt->ctrl_value){
				tmp = pt->matrix[1] * state[pos1];
				state[pos1] = pt->matrix[2] * state[pos0];
				state[pos0] = tmp;
			}
		}
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DGM::CpuExecution3_1(PT *pt, long mem_size){ //Denso
	long pos0, pos1, shift;
	
	shift = 1 << pt->end;

	std::complex <float> tmp;
		
	if (!pt->ctrl_count){ 			//operador não controlado
		mem_size /= 2;
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
			state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];
			state[pos0] = tmp;
		}
	}
	else{					//operador controlado
		vector <long> gap, max;
		long i, c, mask;

		mask = pt->ctrl_mask | shift;

		c = 0;
		for (i = 0; i < qubits; i++){
			if (((mask >> i) & 1) == 0) c++;
			else if (c){
				gap.push_back(1<<(i-c));
				max.push_back(1<<i);
				c = 0;
			}
		}
		if (c){
			gap.push_back(1<<(i-c));
			max.push_back(1<<(qubits+1));
		}
		else{	
			gap.push_back(1<<(qubits+1));
			max.push_back(1<<(qubits+2));
		}

		long pos = 0;

		while (pos < mem_size){
				pos0 = pos | pt->ctrl_value;
				pos1 = pos0 | shift;

				//cout << pos0 <<  " " << pos1 << endl; 

				tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
				state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];			
				state[pos0] = tmp;

				pos += gap[0];
				i = 0;
				while (pos & max[i]){
					pos ^= max[i++];
					pos += gap[i];
				}

		}
		//cout << endl;
	}	
}

void DGM::CpuExecution3_2(PT *pt, long mem_size){ //Diagonal Principal
	long pos0, shift = pt->end;
		
	if (!pt->ctrl_count)	//operador não controlado
		for (long pos = 0; pos < mem_size; pos++)
			state[pos] = pt->matrix[((pos >> shift) & 1) * 3] * state[pos];
	else{					//operador controlado
		vector <long> gap, max;
		long i, c, mask;

		mask = pt->ctrl_mask;

		c = 0;
		for (i = 0; i < qubits; i++){
			if (((mask >> i) & 1) == 0) c++;
			else if (c){
				gap.push_back(1<<(i-c));
				max.push_back(1<<i);
				c = 0;
			}
		}
		if (c){
			gap.push_back(1<<(i-c));
			max.push_back(1<<(qubits+1));
		}
		else{	
			gap.push_back(1<<(qubits+1));
			max.push_back(1<<(qubits+2));
		}

		long pos = 0;

		while (pos < mem_size){
				pos0 = pos | pt->ctrl_value;

				//cout << pos0 << endl; 
				state[pos0] = pt->matrix[((pos0 >> shift) & 1) * 3] * state[pos0];

				pos += gap[0];
				i = 0;
				while (pos & max[i]){
					pos ^= max[i++];
					pos += gap[i];
				}

		}
	}
}

void DGM::CpuExecution3_3(PT *pt, long mem_size){ //Diagonal Secundária
	long pos0, pos1, shift;
	
	shift = 1 << pt->end;


	std::complex <float> tmp;
		
	if (!pt->ctrl_count){ 	//operador não controlado
		mem_size /= 2;
		for (long pos = 0; pos < mem_size; pos++){
			pos0 = (pos * 2) - (pos & (shift-1));
			pos1 = pos0 | shift;

			tmp = pt->matrix[1] * state[pos1];
			state[pos1] = pt->matrix[2] * state[pos0];
			state[pos0] = tmp;
		}
	}
	else{					//operador controlado
		vector <long> gap, max;
		long i, c, mask;

		mask = pt->ctrl_mask | shift;

		c = 0;
		for (i = 0; i < qubits; i++){
			if (((mask >> i) & 1) == 0) c++;
			else if (c){
				gap.push_back(1<<(i-c));
				max.push_back(1<<i);
				c = 0;
			}
		}
		if (c){
			gap.push_back(1<<(i-c));
			max.push_back(1<<(qubits+1));
		}
		else{	
			gap.push_back(1<<(qubits+1));
			max.push_back(1<<(qubits+2));
		}

		long pos = 0;

		while (pos < mem_size){
			pos0 = pos | pt->ctrl_value;
			pos1 = pos0 | shift;

			tmp = pt->matrix[1] * state[pos1];
			state[pos1] = pt->matrix[2] * state[pos0];
			state[pos0] = tmp;

			pos += gap[0];
			i = 0;
			while (pos & max[i]){
				pos ^= max[i++];
				pos += gap[i];
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PCpuExecution1(std::complex <float> *state, PT **pts, int qubits, long n_threads, int coales, int region, int it){
	long i, start, end;
	i = start = 0;
	while (pts[i] != NULL){
		long count = coales;
		long reg_mask = (coales)? (1 << coales) - 1 : 0;

		//Pega os operadores que estão dentro da região coalescida (reg_mask inicial),
		//e acrescenta operadores em qubits fora dela até chegar ao limite da região (region definida)
		start = i;
		while (count < region && pts[i] != NULL){					//Repete enquanto o número de qubits da região não atingir o limite (region) e houver operadores
			if (//pts[i]->matrixType() != DIAG_PRI &&					//O qubit de operadores de diagonal principal não importa para região (sempre podem ser acrescentados)
				!((reg_mask >> pts[i]->end) & 1)){				//Se o qubit do operador estiver fora da região (reg_mask), incrementa o contador de qubits da região
				count++;
			}

			if (count <= region)// && pts[i]->matrixType() != DIAG_PRI)
				reg_mask = reg_mask | (1 << pts[i]->end);			//Acrescenta o qubit do operador na região se ainda não tiver atingido o limite (region)
				
			i++;
		}
		//Segue acerscentado até encontrar um operador que não esteja dentro da região
		while (pts[i] != NULL){
			if (((reg_mask >> pts[i]->end) & 1))// || pts[i]->matrixType() == DIAG_PRI)
				i++;
			else
				break;
		}
		end = i;													//Executa até o operador na posiçao 'i' (exclusive) nesta iteração


		//Se o número de qubits na região (count) não tiver atingido o limite (region),
		//acrescenta os ultimos qubits (final da mascara) à região até completar
		//for (long a = 1<<(qubits-1); count < region; a = a >> 1){
		for (long a = 1; count < region; a = a << 1){
			if (a & ~reg_mask){
				reg_mask = reg_mask | a;
				count++;
			}
		}

		if (count < region)
			region = count;

		long reg_count = (1 << (qubits - region)) + 1; 				//Número de regiões 			-	 +1 para a condição de parada incluir todos
		long pos_count = 1 << (region - 1); 						//Número de posições na região 	-	 -1 porque são duas posições por iteração

		omp_set_num_threads(n_threads);

		long ext_reg_id = 0;	//contador 'global' do número de regiões já computadas

		#pragma omp parallel
		{

			long reg_id;		//indentificador local da região

			//Define a primeira região (reg_id) da thread
			#pragma omp critical (teste)
			{
				reg_id = ext_reg_id;
				ext_reg_id = (ext_reg_id + reg_mask + 1) & ~reg_mask;
				reg_count--;
				if (reg_count <= 0)
					reg_id = -1;
			}

			int print = (omp_get_thread_num()==0);
			
			
			while (reg_id != -1){		
				//Computa os operadores
				PCpuExecution1_0(state, pts, qubits, start, end, pos_count, reg_id, reg_mask);
		
				//Define a próxima região (reg_id) da thread
				#pragma omp critical (teste)
				{
					reg_id = ext_reg_id;
					ext_reg_id = (ext_reg_id + reg_mask + 1) & ~reg_mask;
					reg_count--;
					if (reg_count <= 0)
						reg_id = -1;
				}
			}
		}
	}
}

void PCpuExecution1_0(std::complex <float> *state, PT **pts, int qubits, int start, int end, int pos_count, int reg_id, int reg_mask){
	PT *QG;
	long pos0, pos1;
	std::complex <float> tmp;

	//cout << "\nExecution" << endl;
	//cout << "reg_id " << getBinaryString(reg_id, qubits, true) << endl;
	//cout << "reg_mask " << getBinaryString(reg_mask, qubits, true) << endl;
	//cout << "pos_count " << pos_count << endl;
	//cout << "start " << start << endl;
	//cout << "end " << end << endl;

	for (int op = start; op < end; op++){
		QG = pts[op];
		long shift = (1 << QG->end);						//mascara com a posição do qubit do operador
		long mt = QG->matrixType();
		//if (mt == DIAG_PRI) shift = coalesc;	//se for um operador de diagonal principal, a posição do qubit não é relevante
		long pos_mask = reg_mask & ~shift;			//mascara da posição --- retira o 'shift' da reg_mask, para o 'inc pular sobre ' esse bit também
		long inc = ~pos_mask + 1;						  	//usado para calcular a proxima posição de uma região
		long pos = 0;

		//cout << "OP " << op << endl;
		//cout << "pos_mask " << getBinaryString(pos_mask, qubits, true) << endl;
		//cout << "shift " << getBinaryString(shift, qubits, true) << endl;
		//cout << "inc " << getBinaryString(inc, qubits, true) << endl;
					
		if (!QG->ctrl_count){
			switch (mt){
				case DENSE:
					for (long p = 0; p < pos_count; p++){
						pos0 = pos | reg_id;
						pos1 = pos0 | shift;
						pos = (pos+inc) & pos_mask;

						//cout << getBinaryString(pos0, qubits, true) << " - " << getBinaryString(pos1, qubits, true) << endl;

						tmp 		= QG->matrix[2] * state[pos0] + QG->matrix[3] * state[pos1];
						state[pos0] = QG->matrix[0] * state[pos0] + QG->matrix[1] * state[pos1];
						state[pos1] = tmp;
					}
					break;
				case DIAG_PRI:
					for (long p = 0; p < pos_count; p++){
							pos0 = pos | reg_id;
							pos1 = pos0 | shift;
							pos = (pos+inc) & pos_mask;

							tmp			= QG->matrix[3] * state[pos1];
							state[pos0] *= QG->matrix[0];// * state[pos0];
							state[pos1] = tmp;// * state[pos1];tmp;
					}
					break;
				
				case DIAG_SEC:
					for (long p = 0; p < pos_count; p++){
							pos0 = pos | reg_id;
							pos1 = pos0 | shift;
							pos = (pos+inc) & pos_mask;

							tmp 		= QG->matrix[2] * state[pos0];
							state[pos0] = QG->matrix[1] * state[pos1];
							state[pos1] = tmp;
					}
					break;
				default:
					printf("Erro de Tipo\n");
			}
		}
		//Importante: reg_id é o identificador da região e corresponde ao valor dos qubits externos à região de operação (reg_mask)
		else {			
			if ((QG->ctrl_mask & reg_id & ~reg_mask) == (QG->ctrl_value & ~reg_mask)){		//Verifica se a parte 'global' do controle satisfaz a região (reg_id)

				// É preciso arrumar o reg_mask retirando os qubits de controle que estão dentro da região e arrumar o reg_id para incluir o valor dos controles
				long ctrl_reg_id = reg_id | QG->ctrl_value;				//Esta operação inclui o valor dos controles locais no reg_id (funciona pois os valores globais já deram match)
				long ctrl_reg_mask = reg_mask;							//Valor inicial da mascara da região com controle
				long ctrl_pos_count = pos_count;						//Número inicial de posições a serem calculadas

				for (int i = 0, m = 1; i < qubits; i++, m = m << 1){ 	//percorre os qubits
					if (m & reg_mask & QG->ctrl_mask){					//se o qubit pertencer a região e for um controle:
						ctrl_reg_mask ^= m;								//	remove ele da região(reg_mask) (para não iterar sobre ele)
						ctrl_pos_count /= 2;							//	diminui a quantidade de posições que é preciso calcular.
					}
				}

				pos_mask = ctrl_reg_mask & ~shift;						//mascara da posição --- retira o 'shift' da reg_mask, para o 'inc pular sobre' esse bit também
				inc = ~pos_mask + 1;

				switch (mt){
					case DENSE:
						for (long p = 0; p < ctrl_pos_count; p++){
							pos0 = pos | ctrl_reg_id;
							pos1 = pos0 | shift;
							pos = (pos+inc) & pos_mask;

							tmp 		= QG->matrix[2] * state[pos0] + QG->matrix[3] * state[pos1];
							state[pos0] = QG->matrix[0] * state[pos0] + QG->matrix[1] * state[pos1];
							state[pos1] = tmp;
						}
						break;
					case DIAG_PRI:
						for (long p = 0; p < ctrl_pos_count; p++){
							pos0 = pos | ctrl_reg_id;
							pos1 = pos0 | shift;
							pos = (pos+inc) & pos_mask;

							tmp			= QG->matrix[3] * state[pos1];
							state[pos0] *= QG->matrix[0];
							state[pos1] = tmp;
						}
						break;
					
					case DIAG_SEC:
						for (long p = 0; p < ctrl_pos_count; p++){
							pos0 = pos | ctrl_reg_id;
							pos1 = pos0 | shift;
							pos = (pos+inc) & pos_mask;

							tmp 		= QG->matrix[2] * state[pos0];
							state[pos0] = QG->matrix[1] * state[pos1];
							state[pos1] = tmp;
						}
						break;

					default:
						printf("Erro de Tipo");
				}
			}
		}
	}
}

#ifndef ONLY_CPU

void DGM::HybridExecution(PT **pts){
	long mem_size = pow(2.0, qubits);
	long qubits_limit = 20;
	long global_coales = 15; //(cpu_coales > gpu_coales) ? cpu_coales : gpu_coales;

	long global_region = qubits_limit;
	long global_start, global_end;

	long global_count, global_reg_mask, global_reg_count, global_pos_count, ext_proj_id; 

	omp_set_num_threads(n_threads);

	int i = 0;
	while (pts[i] != NULL){
		global_count = global_coales;
		global_reg_mask = (global_coales)? (1 << global_coales) - 1 : 0;

		//Realiza a projeção dos operadores de acordo com o limite de qubits que podem ser executados
		global_start = i;
		while (global_count < global_region && pts[i] != NULL){			//Repete enquanto o número de qubits da região não atingir o limite (region) e houver operadores
			if (//pts[i]->matrixType() != DIAG_PRI &&					//O qubit de operadores de diagonal principal não importa para região (sempre podem ser acrescentados)
			!((global_reg_mask >> pts[i]->end) & 1)){				
				global_count++;
			}

			if (global_count <= global_region)// && pts[i]->matrixType() != DIAG_PRI)
				global_reg_mask = global_reg_mask | (1 << pts[i]->end);			//Acrescenta o qubit do operador na região se ainda não tiver atingido o limite (region)	

			i++;
		}

		while (pts[i] != NULL){
			if (((global_reg_mask >> pts[i]->end) & 1))// || pts[i]->matrixType() == DIAG_PRI)
				i++;
			else
				break;
		}
		global_end = i;

		//Se o número de qubits na região (count) não tiver atingido o limite (region),
		//acrescenta os ultimos qubits (final da mascara) à região até completar
		//for (long a = 1<<(qubits-1); count < region; a = a >> 1){
		for (long a = 1; global_count < global_region; a = a << 1){
			if (a & ~global_reg_mask){
				global_reg_mask = global_reg_mask | a;
				global_count++;
			}
		}

		if (global_count < global_region)
			global_region = global_count;
	
		global_reg_count = (1 << (qubits - global_region)) + 1; 				//Número de regiões	- +1 para a condição de parada incluir todos
		global_pos_count = 1 << (global_region - 1);

		/////////////////////////////////////////////////////////////////////////////////////////////////////

		ext_proj_id = 0;	//contador 'global' do número de regiões já computadas

		//Define a primeira região (reg_id) da thread

		#pragma omp parallel num_threads(n_threads)
		{
			if (omp_get_thread_num() < n_threads){  //CPU EXECUTION
				long cpu_proj_id;
				
				#pragma omp critical (global_hybrid)
				{
					cpu_proj_id = ext_proj_id;
					ext_proj_id = (ext_proj_id + global_reg_mask + 1) & ~global_reg_mask;
					global_reg_count--;
					if (global_reg_count <= 0)
						cpu_proj_id = -1;
				}
	
				while (cpu_proj_id != -1){
					long cpu_i, cpu_start, cpu_end;

					cpu_start = global_start;
			
					cpu_i = cpu_start;
			
					while (cpu_start < global_end){
						long cpu_count = cpu_coales;
						long cpu_reg_mask = (cpu_coales)? (1 << cpu_coales) - 1 : 0;
			
						while ((cpu_count < cpu_region) && (cpu_i < global_end)){	//Tem que pertencer a região 'global'
							if (!((cpu_reg_mask >> pts[cpu_i]->end) & 1)){			//Se o qubit do operador estiver fora da região (reg_mask), incrementa o contador de qubits da região
								cpu_count++;
							}
		
							if (cpu_count <= cpu_region)// && pts[i]->matrixType() != DIAG_PRI)
								cpu_reg_mask = cpu_reg_mask | (1 << pts[cpu_i]->end);	//Acrescenta o qubit do operador na região se ainda não tiver atingido o limite (region)
						
							cpu_i++;
						}
			
						while (cpu_i < global_end){
							if (((cpu_reg_mask >> pts[cpu_i]->end) & 1))// || pts[i]->matrixType() == DIAG_PRI)
								cpu_i++;
							else
								break;
						}
						cpu_end = cpu_i;
			
						for (long a = 1; cpu_count < cpu_region; a = a << 1){
							if ((a & global_reg_mask) && (a & ~cpu_reg_mask)){ //tem que não estar na região da cpu e estar na global
								cpu_reg_mask = cpu_reg_mask | a;
								cpu_count++;
							}
						}
	
						long cpu_reg_count = (1 << (global_region - cpu_region)) + 1; 		//Número de regiões				-	+1 para a condição de parada incluir todos
						long cpu_pos_count = 1 << (cpu_region - 1); 						//Número de posições na região	-	-1 porque são duas posições por iteração

				
						long cpu_ext_proj_id = 0;
						long inc_ext_proj_id = ~(cpu_reg_mask ^ global_reg_mask) & ((1 << qubits) - 1);
			
						long proj_id;		//indentificador local da região
						proj_id = cpu_ext_proj_id | cpu_proj_id;
						cpu_ext_proj_id = (cpu_ext_proj_id + inc_ext_proj_id + 1) & ~inc_ext_proj_id;
						cpu_reg_count--;
						
						while (proj_id != -1){
							//Computa os operadores
							PCpuExecution1_0(state, pts, qubits, cpu_start, cpu_end, cpu_pos_count, proj_id, cpu_reg_mask);
				
							proj_id = cpu_ext_proj_id | cpu_proj_id;
							cpu_ext_proj_id = (cpu_ext_proj_id + inc_ext_proj_id + 1) & ~inc_ext_proj_id;
							cpu_reg_count--;
							if (cpu_reg_count <= 0)
								proj_id = -1;
						}
			
						cpu_start = cpu_end;
					}
		
					#pragma omp critical (global_hybrid)
					{
						cpu_proj_id = ext_proj_id;
						ext_proj_id = (ext_proj_id + global_reg_mask + 1) & ~global_reg_mask;
						global_reg_count--;
							if (global_reg_count <= 0)
						cpu_proj_id = -1;
					}
				}
				
			}
			//#pragma omp section          //GPU EXECUTION
			else{
				#ifndef ONLY_CPU
					long gpu_proj_id;
					
					#pragma omp critical (global_hybrid)
					{
						gpu_proj_id = ext_proj_id;
						ext_proj_id = (ext_proj_id + global_reg_mask + 1) & ~global_reg_mask;
						global_reg_count--;
						if (global_reg_count <= 0)
							gpu_proj_id = -1;
					}

					while (gpu_proj_id != -1){
						//Project Gates
						vector <PT*> gpu_pts;
						
						int gpu_i;

						int map_qb[qubits];
						memset(map_qb, -1, qubits * sizeof(int));
			
						int m = 0;
						for (gpu_i = 0; gpu_i < qubits; gpu_i++){
							if ((1 << gpu_i) & global_reg_mask){
								map_qb[gpu_i] = m++;
							}
						}
						
						PT *aux;
						gpu_pts.clear();
						for (int gpu_i = global_start; gpu_i < global_end; gpu_i++){
							
							//verifica se o controle do operador satisfaz a parte global da região
							if ((pts[gpu_i]->ctrl_mask & gpu_proj_id & ~global_reg_mask) == (pts[gpu_i]->ctrl_value & ~global_reg_mask)){
								aux = new PT();

								aux->qubits = pts[gpu_i]->qubits;

								aux->matrix = pts[gpu_i]->matrix;
								aux->mat_size = pts[gpu_i]->mat_size;
								aux->ctrl_mask = pts[gpu_i]->ctrl_mask & global_reg_mask;
								aux->ctrl_value = pts[gpu_i]->ctrl_value & global_reg_mask;

								aux->end = map_qb[pts[gpu_i]->end];
								aux->start = aux->end - log2(aux->mat_size);

								aux->ctrl_count = 0;
								for (int c = global_coales; c < qubits; c++){
									if (aux->ctrl_mask & (1<<c)){
										aux->ctrl_count++;

										aux->ctrl_mask &= ~(1<<c);			//retira da mascara o controle do qubit atual (c)
										aux->ctrl_mask |= (1 << map_qb[c]);	//e coloca o qubit que ele mapeia (map_qb[c])

										if (aux->ctrl_value & (1<<c)){ 		//se o valor do controle for zero faz a mesma coisa para ctrl_value;
											aux->ctrl_mask &= ~(1<<c);
											aux->ctrl_mask |= (1 << map_qb[c]);
										}
									}
								}	

								gpu_pts.push_back(aux);
							}
						}
						gpu_pts.push_back(NULL);
						////////////////

						ProjectState(state, qubits, global_region, gpu_proj_id, global_reg_mask, multi_gpu);

						GpuExecutionWrapper(NULL, &gpu_pts[0], global_region, gpu_coales, gpu_region, multi_gpu, tam_block, rept, 1);
		
						GetState(state, qubits, global_region, gpu_proj_id, global_reg_mask, multi_gpu);

						for (int c = 0; c < gpu_pts.size() - 1; c++){
							delete gpu_pts[c];
						}
			
						#pragma omp critical (global_hybrid)
						{
							gpu_proj_id = ext_proj_id;
							ext_proj_id = (ext_proj_id + global_reg_mask + 1) & ~global_reg_mask;
							global_reg_count--;
							if (global_reg_count <= 0)
								gpu_proj_id = -1;
						}
					}
				#endif
			}
		//}
		}
	}
}

void DGM::HybridExecution2(PT **pts){
	long mem_size = pow(2.0, qubits);
	long qubits_limit = 20;
	long global_coales = 11; //(cpu_coales > gpu_coales) ? cpu_coales : gpu_coales;

	long global_region = qubits_limit;
	long global_start, global_end;

	long global_count, global_reg_mask, global_reg_count, global_pos_count, ext_proj_id;

	cout << "HybridExecution2" << endl;
	cout << "n_threads: " << n_threads << endl;
	cout << "qubits: " << qubits << endl;
	cout << "global_region: " << global_region << endl;
	cout << "global_coales: " << global_coales << endl;
	cout << "cpu_coales: " << cpu_coales << endl;
	cout << "cpu_region: " << cpu_region << endl;

	omp_set_num_threads(n_threads);

	Barrier cpu_barrier(n_threads);

	while (pts[global_start] != NULL){		
		Projection global_proj;
		global_proj.setData(nullptr, pts, global_start, -1, qubits, global_coales, global_region);

		global_end = global_start + global_proj.operators_count;

		//cout << "\n#######\nGLOBAL PROJ " << global_start << " - " << global_end << endl;
		//global_proj.printInfo();
		//cout << "#######" << endl;

		long global_cpu_proj_id;

		Projection cpu_proj;
		long cpu_start, cpu_end, cpu_pos_count;

		#pragma omp parallel num_threads(n_threads)
		{
			long thread_id = omp_get_thread_num();
			//cout << string("ThreadID " + std::to_string(thread_id) + "\n");
			if (thread_id < n_threads){  //CPU EXECUTION		
				// get first global projection id	
				if (thread_id == 0){ // master thread for cpu
					global_cpu_proj_id = global_proj.getNextProjectionId();
					cpu_end = global_start; // set start operator

					//cout << string("\n#######\nFIRST GLOBAL CPU PROJ ID: " + getBinaryString(global_cpu_proj_id, qubits, true) + " - " + std::to_string(global_proj.count) + "\n");
				}

				// wait for all cpu threads
				cpu_barrier.arrive_and_wait(1);
	
				while (global_cpu_proj_id != -1){
					cpu_barrier.arrive_and_wait(2);
					// get first cpu sub projection from global projection

					while (cpu_end < global_end) {
						cpu_barrier.arrive_and_wait(3);
						if (thread_id == 0){ // master thread for cpu
							cpu_start = cpu_end;
							cpu_proj.setData(&global_proj, pts, cpu_start, global_end, qubits, cpu_coales, cpu_region);
							cpu_end = cpu_start + cpu_proj.operators_count;

							cpu_pos_count = 1 << (cpu_proj.region_size);

							//cout << string("\n#######\nSUB CPU PROJ " + std::to_string(cpu_start) + " - " + std::to_string(cpu_end) + "\n");
							//cpu_proj.printInfo();
						}
						cpu_barrier.arrive_and_wait(4);
						////////////////////////////////////////////////////
						long cpu_proj_id = cpu_proj.getNextProjectionId();
						
						while (cpu_proj_id != -1){
							//cout << string("ThreadId " + std::to_string(thread_id) + ": cpu_proj_id " + getBinaryString(cpu_proj_id, qubits, true) + " - " + std::to_string(cpu_proj.count) + "\n");
							//Computa os operadores
							PCpuExecution1_0(state, pts, qubits, cpu_start, cpu_end, cpu_pos_count / 2, cpu_proj_id, cpu_proj.region_mask);

							cpu_proj_id = cpu_proj.getNextProjectionId();
						}

						//cpu_barrier.arrive_and_wait();

						////////////////////////////////////////////////////
					}

					cpu_barrier.arrive_and_wait(5);
					// get next global projection id
					if (thread_id == 0){ // master thread for cpu
						global_cpu_proj_id = global_proj.getNextProjectionId();
						cpu_end = global_start; // reset start operator

						//cout << string("\n#######\nGLOBAL CPU PROJ ID " + getBinaryString(global_cpu_proj_id, qubits, true) + " - " + std::to_string(global_proj.count) + "\n");
					}
					cpu_barrier.arrive_and_wait(6);

				}
			}
			//#pragma omp section          //GPU EXECUTION
			else{
				cout << "\n#######\nGPU EXECUTION" << endl;
				long gpu_proj_id = global_proj.getNextProjectionId();

				while (gpu_proj_id != -1){
					//Project Gates
					vector <PT*> gpu_pts;
					
					int gpu_i;

					int map_qb[qubits];
					memset(map_qb, -1, qubits * sizeof(int));
		
					int m = 0;
					for (gpu_i = 0; gpu_i < qubits; gpu_i++){
						if ((1 << gpu_i) & global_reg_mask){
							map_qb[gpu_i] = m++;
						}
					}
					
					PT *aux;
					gpu_pts.clear();
					for (int gpu_i = global_start; gpu_i < global_end; gpu_i++){
						
						//verifica se o controle do operador satisfaz a parte global da região
						if ((pts[gpu_i]->ctrl_mask & gpu_proj_id & ~global_reg_mask) == (pts[gpu_i]->ctrl_value & ~global_reg_mask)){
							aux = new PT();

							aux->qubits = pts[gpu_i]->qubits;

							aux->matrix = pts[gpu_i]->matrix;
							aux->mat_size = pts[gpu_i]->mat_size;
							aux->ctrl_mask = pts[gpu_i]->ctrl_mask & global_reg_mask;
							aux->ctrl_value = pts[gpu_i]->ctrl_value & global_reg_mask;

							aux->end = map_qb[pts[gpu_i]->end];
							aux->start = aux->end - log2(aux->mat_size);

							aux->ctrl_count = 0;
							for (int c = global_coales; c < qubits; c++){
								if (aux->ctrl_mask & (1<<c)){
									aux->ctrl_count++;

									aux->ctrl_mask &= ~(1<<c);			//retira da mascara o controle do qubit atual (c)
									aux->ctrl_mask |= (1 << map_qb[c]);	//e coloca o qubit que ele mapeia (map_qb[c])

									if (aux->ctrl_value & (1<<c)){ 		//se o valor do controle for zero faz a mesma coisa para ctrl_value;
										aux->ctrl_mask &= ~(1<<c);
										aux->ctrl_mask |= (1 << map_qb[c]);
									}
								}
							}	

							gpu_pts.push_back(aux);
						}
					}
					gpu_pts.push_back(NULL);
					////////////////

					ProjectState(state, qubits, global_region, gpu_proj_id, global_reg_mask, multi_gpu);

					GpuExecutionWrapper(NULL, &gpu_pts[0], global_region, gpu_coales, gpu_region, multi_gpu, tam_block, rept, 1);
	
					GetState(state, qubits, global_region, gpu_proj_id, global_reg_mask, multi_gpu);

					for (int c = 0; c < gpu_pts.size() - 1; c++){
						delete gpu_pts[c];
					}
		
					gpu_proj_id = global_proj.getNextProjectionId();
				}
			}
		//}
		}

		global_start = global_end;
	}
}

#endif

void DGM::setCpuStructure(long cpu_region, long cpu_coales){
	this->cpu_region = cpu_region;
	this->cpu_coales = cpu_coales;
}

void DGM::setGpuStructure(long gpu_region, long gpu_coales, int rept){
	this->gpu_region = gpu_region;
	this->gpu_coales = gpu_coales;
	this->rept = rept;
	this->tam_block = 1 << gpu_region / 2 / rept;
}

void Projection::setData(Projection *parent_proj, PT **pts, long start, long end, long qubits, long coales, long region_size, bool include_main_diag) {
	region_mask = (coales > 0) ? (1 << coales) - 1 : 0; // mascara da região, inicializada com os quibts de coalescencia
	
	long count = coales; // conta o número de qubits na região

	//Realiza a projeção dos operadores de acordo com o limite de qubits que podem ser executados

	long i = start;
	while (count < region_size && (end < 0 || i < end) && pts[i] != NULL){	//Repete enquanto o número de qubits da região não atingir o limite (region) e houver operadores
		if (!(include_main_diag && pts[i]->matrixType() == DIAG_PRI)) { // verifica se pode incluir operadores de diagonal principal automaticamente, pois em alguns cenários eles não importam para a projeção e sempre podem ser acrescentados
			if (!((region_mask >> pts[i]->end) & 1)) { // se o qubit do operador não estiver sido acrescentado ainda, então incrementa o contador de qubits da região
				count++;

				if (count <= region_size) // se não tiver atingido o limite
					region_mask = region_mask | (1 << pts[i]->end);			// Acrescenta o qubit do operador na região se ainda não tiver atingido o limite (region)
			}
		}

		i++;
	}

	// segue acrescentando operadores enquanto eles estiverem dentro da região
	while ((end < 0 || i < end) && pts[i] != NULL){
		if ((include_main_diag && pts[i]->matrixType() == DIAG_PRI) || ((region_mask >> pts[i]->end) & 1))
			i++;
		else
			break;
	}

	//Se o número de qubits na região (count) não tiver atingido o limite (region_size),
	//acrescenta os ultimos qubits (final da mascara) à região até completar
	//for (long a = 1<<(qubits-1); count < region; a = a >> 1){
	for (long b = 1; count < region_size; b = b << 1){
		if ((b & ~region_mask) && (!parent_proj || (b & parent_proj->region_mask))){
			region_mask = region_mask | b;
			count++;
		}
	}

	this->qubits = qubits;
	this->region_size = region_size;
	this->operators_count = i - start;
	this->count = 0;
	this->cur_id = 0;
	this->next_id = 0;

	if (parent_proj) {
		this->total = (1 << (parent_proj->region_size - region_size));

		this->parent_proj_id = parent_proj->cur_id;
		this->inc_mask = (region_mask | (~parent_proj->region_mask)) & ((1 << qubits) - 1);
	}
	else {
		this->total = (1 << (qubits - region_size));

		this->parent_proj_id = 0;
		this->inc_mask = region_mask;
	}
}

long Projection::getNextProjectionId() {
	std::lock_guard<std::mutex> lock(mutex_);

	count++;
	if (count > total)
		return -1;

	cur_id = next_id;
	next_id = ((cur_id + inc_mask + 1) & ~inc_mask); // & ((1 << qubits) - 1);

	return cur_id | parent_proj_id;
}

void Projection::printInfo() {
	cout << "Projection Info" << endl;
	cout << "region_size: " << region_size << endl;
	cout << "operators_count: " << operators_count << endl;
	cout << "count: " << count << endl;
	cout << "total: " << total << endl;
	cout << "cur_id: " << getBinaryString(cur_id, qubits, true) << endl;
	cout << "inc_mask: " << getBinaryString(inc_mask, qubits, true) << endl;
	cout << "region_mask: " << getBinaryString(region_mask, qubits, true) << endl;
	cout << "parent_proj_id: " << getBinaryString(parent_proj_id, qubits, true) << endl;
}

std::string getBinaryString(long num, int n, bool includeNum) {
    if (n < 1) {
        return "Invalid number of bits";
    }

    std::bitset<sizeof(long) * 8> bits(num); // Create a bitset with the integer value

    // Convert the bitset to a string
    std::string binaryString = bits.to_string();

    // Ensure the string has 'n' bits
    if (binaryString.length() < n) {
        // Pad the string with leading zeros
        binaryString = std::string(n - binaryString.length(), '0') + binaryString;
    } else if (binaryString.length() > n) {
        // Trim excess bits if the string has more than 'n' bits
        binaryString = binaryString.substr(binaryString.length() - n);
    }

	if (includeNum)
		binaryString = std::to_string(num) + " (" + binaryString + ")";

    return binaryString;
}

void report_num_threads(int level){
	#pragma omp single
	{
		printf("Level %d: number of threads in the team - %d\n", level, omp_get_num_threads());
	}
}