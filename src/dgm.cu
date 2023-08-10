#include <omp.h>
#include <unistd.h>

#include <cstdio>
#include <iostream>
#include <iterator>

#include "dgm.h"
#include "gpu.h"
#include "pcpu.h"

void Tokenize(const std::string &str, std::vector<std::string> &tokens,
              const std::string &delimiters = ",") {
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the std::vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////

float complex *GenericExecute(float complex *state, std::string function,
                              int qubits, int type, int threads,
                              int factor = 0) {
  DGM dgm;
  dgm.exec_type = type;
  dgm.cpu_params.n_threads = threads;
  dgm.qubits = qubits;
  dgm.factor = factor;

  dgm.setMemory(state);

  dgm.executeFunction(function);

  state = dgm.state;

  dgm.state = NULL;

  return state;
}

float complex *GenericExecute(float complex *state,
                              std::vector<std::string> function, int qubits,
                              int type, int threads, int factor = 0) {
  DGM dgm;
  dgm.exec_type = type;
  dgm.cpu_params.n_threads = threads;
  dgm.qubits = qubits;
  dgm.factor = factor;
  dgm.setMemory(state);

  dgm.executeFunction(function);

  dgm.state = NULL;

  return state;
}

///////////////////////////////////////////////////////////////////////////////////////////////

DGM::DGM() {
  MAX_QB = QB_LIMIT;
  MAX_PT = PT_TAM;

  pts = NULL;
  state = NULL;
  en_print = false;
  exec_type = t_CPU;
  factor = 1;
  gpu_params.multi_gpu = 1;
}

DGM::~DGM() { erase(); }

void DGM::setExecType(int type) { exec_type = type; }

void DGM::printPTs() {
  for (int i = 0; i < vec_pts.size() - 1; i++) {
    vec_pts[i]->print();
  }
}

void DGM::erase() {
  if (!pts) return;

  long i = 0;
  while (pts[i] != NULL) {
    pts[i]->destructor();
    free(pts[i]);
    i++;
  }

  vec_pts.clear();
  pts = NULL;
}

void DGM::allocateMemory() {
  state = (float complex *)calloc(pow(2, qubits), sizeof(float complex));
}

void DGM::setMemory(float complex *mem) {
  freeMemory();
  state = mem;
}

void DGM::freeMemory() {
  if (state) free(state);
  state = NULL;
}

void DGM::setMemoryValue(int pos) { state[pos] = 1; }

int DGM::measure(int q_pos) {
  long size = pow(2.0, qubits);

  long shift = (qubits - 1 - q_pos);

  int count_one, count_zero, num_pb;
  float zero, one, norm, r;
  one = zero = 0;

  // #pragma omp for;
  for (long i = 0; i < size; i++) {
    if ((i >> shift) & 1)
      one += pow(crealf(state[i]), 2.0) + pow(cimagf(state[i]), 2.0);
    else
      zero += pow(crealf(state[i]), 2.0) + pow(cimagf(state[i]), 2.0);
  }

  long m;
  srand(time(NULL));
  count_one = 0;
  count_zero = 0;
  num_pb = 1;

  for (int i = 0; i < num_pb; i++) {
    r = (double)rand() / RAND_MAX;
    if (zero > r)
      count_zero++;
    else
      count_one++;
  }

  if (count_one > count_zero) {
    measure_value = one;
    norm = sqrt(one);
    m = 1;
  } else {
    measure_value = zero;
    norm = sqrt(zero);
    m = 0;
  }

  long mask;
  mask = pow(2, shift) - 1;
#pragma omp for
  for (long i = 0; i < size / 2; i++) {
    long pos = (i << 1) - (i & mask);
    state[pos] = state[pos | (m << shift)] / norm;
    state[pos | (1 << shift)] = 0.0;
  }

  return m;
}

void DGM::colapse(int q_pos, int value) {
  long size = pow(2.0, qubits);
  long mask = (qubits - 1 - q_pos);

  float m;
  m = 0;

  for (long i = 0; i < size; i++)
    if (((i >> mask) & 1) == value)
      m += pow(crealf(state[i]), 2.0) + pow(cimagf(state[i]), 2.0);

  std::cout << m << std::endl;

  m = sqrt(m);
  for (long i = 0; i < size; i++) {
    if (((i >> mask) & 1) == value)
      state[i] = state[i] / m;
    else
      state[i] = 0.0;
  }
}

std::map<long, float> DGM::measure(std::vector<int> q_pos) {
  long mask = 0;

  for (int i = 0; i < q_pos.size(); i++)
    mask = mask | (1 << (qubits - 1 - q_pos[i]));

  std::map<long, float> m;

  long size = pow(2.0, qubits);

  for (long i = 0; i < size; i++)
    m[i & mask] += pow(crealf(state[i]), 2.0) + pow(cimagf(state[i]), 2.0);

  return m;
}

void DGM::setFunction(std::string function, int it, bool er) {
  std::vector<std::string> steps;

  Tokenize(function, steps, ";");

  setFunction(steps, it, er);
}

void DGM::setFunction(std::vector<std::string> steps, int it, bool er) {
  if (er)
    erase();
  else
    vec_pts.pop_back();

  std::vector<PT *> step_pts, vec_tmp;
  std::map<long, Group> gps;

  for (long j = 0; j < it; j++)
    for (long i = 0; i < steps.size(); i++) {
      gps = genGroups(steps[i]);
      genPTs(gps, step_pts);

      if (i % 2)
        sort(step_pts.begin(), step_pts.end(), increasing);
      else
        sort(step_pts.begin(), step_pts.end(), decreasing);

      vec_pts.insert(vec_pts.end(), step_pts.begin(), step_pts.end());
    }

  vec_pts.push_back(NULL);

  pts = &vec_pts[0];
}

std::map<long, Group> DGM::genGroups(std::string step) {
  std::vector<std::string> ops;
  Tokenize(step, ops);  // separa os operadores usando "," como delimitador
  qubits = ops.size();

  size_t found_c, found_t, p;
  std::string str;
  long pos, ctrl_value, ctrl_num;

  std::map<long, Group> gps;

  char *pEnd;
  pos = 0;
  std::vector<std::string>::iterator it;
  for (it = ops.begin(); it != ops.end(); ++it) {  // percorre os operadores
    str = *it;
    // std::cout << str << std::endl;
    found_c = str.find("Control");  // tamanho 7
    found_t = str.find("Target");   // tamanho 6
    p = str.find("(") + 1;

    if (found_c != std::string::npos) {  // Controle
      ctrl_num = strtol(str.c_str() + 7, &pEnd, 10);
      ctrl_value = strtol(str.c_str() + p, &pEnd, 10);

      gps[ctrl_num].ctrl.push_back(ctrl_value);  // adicona o valor do controle
      gps[ctrl_num].pos_ctrl.push_back(
          pos);  // e a sua posição ao map relacionado ao controle
    } else if (found_t != std::string::npos) {  // Target
      ctrl_num = strtol(str.c_str() + 6, &pEnd, 10);
      str = str.substr(p, str.size() - p - 1);

      gps[ctrl_num].ops.push_back(str);  // adicona o operador
      gps[ctrl_num].pos_ops.push_back(
          pos);           // e a sua posição ao map relacionado ao target
    } else {              // operador normal
      if (str != "ID") {  // se for ID ignora
        gps[0].ops.push_back(str);      // adiciona o operador
        gps[0].pos_ops.push_back(pos);  // e a sua posição ao map '0'
      }
    }
    pos++;
  }

  return gps;
}

void DGM::genPTs(std::map<long, Group> &gps, std::vector<PT *> &step_pts) {
  step_pts.clear();
  Gates gates;

  std::map<long, Group>::iterator it;
  Group gp;
  PT *pt;
  long ctrl_mask, ctrl_value, ctrl_count;
  long size;

  for (it = gps.begin(); it != gps.end(); ++it) {  // percorre os grupos
    gp = it->second;
    size = gp.ops.size();

    ctrl_count = gp.ctrl.size();
    ctrl_value = ctrl_mask = 0;

    for (long i = 0; i < ctrl_count;
         i++) {  // gera a mascara e o valor do controle (em binario)
      gp.pos_ctrl[i] = qubits - gp.pos_ctrl[i] - 1;
      ctrl_mask += (1 << gp.pos_ctrl[i]);
      if (gp.ctrl[i]) ctrl_value += (1 << gp.pos_ctrl[i]);
    }

    for (int p = 0; p < size; p++) {
      pt = (PT *)malloc(sizeof(PT));
      pt->affected = false;

      pt->qubits = 1;
      pt->start = qubits - gp.pos_ops[p];
      pt->end = pt->start - 1;
      pt->mat_size = 2;

      pt->matrix = gates.getMatrix(gp.ops[p]);

      pt->ctrl_value = ctrl_value;
      pt->ctrl_mask = ctrl_mask;
      pt->ctrl_count = ctrl_count;

      if (ctrl_count) {
        pt->ctrl_pos = (long *)malloc(sizeof(long) * ctrl_count);
        copy(gp.pos_ctrl.begin(), gp.pos_ctrl.end(), pt->ctrl_pos);
      }

      step_pts.push_back(pt);
    }
  }
}

void DGM::genMatrix(float complex *matrix,
                    std::vector<float complex *> &matrices, long tam,
                    long current, long line, long column, float complex cmplx) {
  if (cmplx == 0.0) return;

  if (current == tam) {  // percorreu até a ultima matriz
    matrix[line * (1 << tam) + column] = cmplx;
    return;
  }

  for (long l = 0; l < 2; l++)
    for (long c = 0; c < 2; c++)
      genMatrix(matrix, matrices, tam, current + 1, (line << 1) | l,
                (column << 1) | c, cmplx * matrices[current][l * 2 + c]);
}

void DGM::executeFunction(std::vector<std::string> function, int it) {
  setFunction(function);
  execute(it);
}

void DGM::executeFunction(std::string function, int it) {
  if (function == "") return;

  setFunction(function);
  execute(it);
}

float complex *DGM::execute(int it) {
  float complex *result = state;

  switch (exec_type) {
    case t_CPU:
      CpuExecution1(it);
      break;
    case t_PAR_CPU:
      PCpuExecution1(state, pts, qubits, cpu_params.n_threads,
                     cpu_params.cpu_coales, cpu_params.cpu_region, it);
      break;
    case t_GPU:
      result = GpuExecutionWrapper(state, pts, qubits, gpu_params.gpu_coales,
                                   gpu_params.gpu_region, gpu_params.multi_gpu,
                                   gpu_params.tam_block, gpu_params.rept, it);
      break;
    case t_HYBRID:
      HybridExecution(pts);
      break;
    // case t_DIST:
    //   DistributedExecution(state, pts, qubits, n_threads, cpu_coales,
    //                        cpu_region, it);
    //   break;
    default:
      std::cout << "Erro exec type" << std::endl;
      exit(1);
  }

  return result;
}

OPSCounter DGM::CountOps(int it) {
  OPSCounter counter;

  for (int i = 0; pts[i] != NULL; i++) {
    long mt = pts[i]->matrixType();
    switch (mt) {
      case DENSE:
        (pts[i]->ctrl_mask) ? counter.c_dense++ : counter.dense++;
        break;
      case DIAG_PRI:
        (pts[i]->ctrl_mask) ? counter.c_main_diag++ : counter.main_diag++;
        break;
      case DIAG_SEC:
        (pts[i]->ctrl_mask) ? counter.c_sec_diag++ : counter.sec_diag++;
        break;
      default:
        std::cout << "Error on operator type" << std::endl;
        exit(1);
    }
  }

  counter.dense *= it;
  counter.c_dense *= it;
  counter.main_diag *= it;
  counter.c_main_diag *= it;
  counter.sec_diag *= it;
  counter.c_sec_diag *= it;

  counter.total_op = counter.dense + counter.c_dense + counter.main_diag +
                     counter.c_main_diag + counter.sec_diag +
                     counter.c_sec_diag;
  return counter;
}

void DGM::CpuExecution1(int it) {
  long mem_size = pow(2.0, qubits);

  for (int x = 0; x < it; x++) {
    long i = 0;
    while (pts[i] != NULL) {
      long mt = pts[i]->matrixType();

      switch (mt) {
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

void DGM::CpuExecution1_1(PT *pt, long mem_size) {  // Denso
  long pos0, pos1, shift;

  shift = 1 << pt->end;

  float complex tmp;

  if (!pt->ctrl_count) {  // operador não controlado
    mem_size /= 2;
    for (long pos = 0; pos < mem_size; pos++) {
      pos0 = (pos * 2) - (pos & (shift - 1));
      pos1 = pos0 | shift;

      tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
      state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];
      state[pos0] = tmp;
    }
  } else {  // operador controlado
    long mask = ~(pt->ctrl_mask | shift);
    long inc = (~mask) + 1;

    for (long pos = 0; pos < mem_size; pos = (pos + inc) & mask) {
      pos0 = pos | pt->ctrl_value;
      pos1 = pos0 | shift;

      tmp = pt->matrix[0] * state[pos0] + pt->matrix[1] * state[pos1];
      state[pos1] = pt->matrix[2] * state[pos0] + pt->matrix[3] * state[pos1];
      state[pos0] = tmp;
    }
  }
}

void DGM::CpuExecution1_2(PT *pt, long mem_size) {  // Diagonal Principal
  long pos0, shift = pt->end;

  if (!pt->ctrl_count)  // operador não controlado
    for (long pos = 0; pos < mem_size; pos++)
      state[pos] = pt->matrix[((pos >> shift) & 1) * 3] * state[pos];
  else {  // operador controlado
    long mask = ~(pt->ctrl_mask);
    long inc = (~mask) + 1;

    for (long pos = 0; pos < mem_size; pos = (pos + inc) & mask) {
      pos0 = pos | pt->ctrl_value;

      state[pos0] = pt->matrix[((pos0 >> shift) & 1) * 3] * state[pos0];
    }
  }
}

void DGM::CpuExecution1_3(PT *pt, long mem_size) {  // Diagonal Secundária
  long pos0, pos1, shift;

  shift = 1 << pt->end;

  float complex tmp;

  if (!pt->ctrl_count) {  // operador não controlado
    mem_size /= 2;
    for (long pos = 0; pos < mem_size; pos++) {
      pos0 = (pos * 2) - (pos & (shift - 1));
      pos1 = pos0 | shift;

      tmp = pt->matrix[1] * state[pos1];
      state[pos1] = pt->matrix[2] * state[pos0];
      state[pos0] = tmp;
    }
  } else {  // operador controlado
    long mask = ~(pt->ctrl_mask | shift);
    long inc = (~mask) + 1;

    for (long pos = 0; pos < mem_size; pos = (pos + inc) & mask) {
      pos0 = pos | pt->ctrl_value;
      pos1 = pos0 | shift;

      tmp = pt->matrix[1] * state[pos1];
      state[pos1] = pt->matrix[2] * state[pos0];
      state[pos0] = tmp;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void report_num_threads(int level) {
#pragma omp single
  {
    printf("Level %d: number of threads in the team - %d\n", level,
           omp_get_num_threads());
  }
}

void DGM::HybridExecution(PT **pts) {
  long mem_size = pow(2.0, qubits);
  long qubits_limit = 20;
  long global_coales =
      15;  //(cpu_coales > gpu_coales) ? cpu_coales : gpu_coales;

  long global_region = qubits_limit;
  long global_start, global_end;

  long global_count, global_reg_mask, global_reg_count,
      ext_proj_id;  //, global_pos_count; //atualmente não utilizada

  omp_set_num_threads(cpu_params.n_threads);

  int i = 0;
  while (pts[i] != NULL) {
    global_count = global_coales;
    global_reg_mask = (global_coales) ? (1 << global_coales) - 1 : 0;

    // Realiza a projeção dos operadores de acordo com o limite de qubits que
    // podem ser executados
    global_start = i;
    while (global_count < global_region &&
           pts[i] !=
               NULL) {  // Repete enquanto o número de qubits da região não
                        // atingir o limite (region) e houver operadores
      if (              // pts[i]->matrixType() != DIAG_PRI &&
            // //O qubit de operadores de diagonal principal não importa para
            // região (sempre podem ser acrescentados)
          !((global_reg_mask >> pts[i]->end) & 1)) {
        global_count++;
      }

      if (global_count <=
          global_region)  // && pts[i]->matrixType() != DIAG_PRI)
        global_reg_mask =
            global_reg_mask |
            (1 << pts[i]->end);  // Acrescenta o qubit do operador na região se
                                 // ainda não tiver atingido o limite (region)

      i++;
    }

    while (pts[i] != NULL) {
      if (((global_reg_mask >> pts[i]->end) &
           1))  // || pts[i]->matrixType() == DIAG_PRI)
        i++;
      else
        break;
    }
    global_end = i;

    // Se o número de qubits na região (count) nãoo tiver atingido o limite
    // (region), acrescenta os ultimos qubits (final da mascara) à região até
    // completar for (long a = 1<<(qubits-1); count < region; a = a >> 1){
    for (long a = 1; global_count < global_region; a = a << 1) {
      if (a & ~global_reg_mask) {
        global_reg_mask = global_reg_mask | a;
        global_count++;
      }
    }

    if (global_count < global_region) global_region = global_count;

    global_reg_count = (1 << (qubits - global_region)) +
                       1;  // Número de regiões	- +1 para a condição de parada
                           // incluir todos
    // global_pos_count = 1 << (global_region - 1); // Atualmente não utilizada

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    ext_proj_id = 0;  // contador 'global' do número de regiões já computadas

    // Define a primeira região (reg_id) da thread

#pragma omp parallel num_threads(n_threads)
    {
      if (omp_get_thread_num() != 0) {  // CPU EXECUTION
        long cpu_proj_id;

#pragma omp critical(global_teste)
        {
          cpu_proj_id = ext_proj_id;
          ext_proj_id = (ext_proj_id + global_reg_mask + 1) & ~global_reg_mask;
          global_reg_count--;
          if (global_reg_count <= 0) cpu_proj_id = -1;
        }

        while (cpu_proj_id != -1) {
          long cpu_i, cpu_start, cpu_end;

          cpu_start = global_start;

          cpu_i = cpu_start;

          while (cpu_start < global_end) {
            long cpu_count = cpu_params.cpu_coales;
            long cpu_reg_mask =
                (cpu_params.cpu_coales) ? (1 << cpu_params.cpu_coales) - 1 : 0;

            while (
                (cpu_count < cpu_params.cpu_region) &&
                (cpu_i < global_end)) {  // Tem que pertencer a região 'global'
              if (!((cpu_reg_mask >> pts[cpu_i]->end) &
                    1)) {  // Se o qubit do operador estiver fora da região
                           // (reg_mask), incrementa o contador de qubits da
                           // região
                cpu_count++;
              }

              if (cpu_count <= cpu_params.cpu_region)
              // && pts[i]->matrixType() != DIAG_PRI)
              {
                cpu_reg_mask = cpu_reg_mask | (1 << pts[cpu_i]->end);
              }
              // Acrescenta o qubit do operador
              // na região se ainda não tiver
              // atingido o limite (region)

              cpu_i++;
            }

            while (cpu_i < global_end) {
              if (((cpu_reg_mask >> pts[cpu_i]->end) &
                   1))  // || pts[i]->matrixType() == DIAG_PRI)
                cpu_i++;
              else
                break;
            }
            cpu_end = cpu_i;

            for (long a = 1; cpu_count < cpu_params.cpu_region; a = a << 1) {
              if ((a & global_reg_mask) &&
                  (a & ~cpu_reg_mask)) {  // tem que não estar na região da cpu
                                          // e estar na global
                cpu_reg_mask = cpu_reg_mask | a;
                cpu_count++;
              }
            }

            long cpu_reg_count =
                (1 << (global_region - cpu_params.cpu_region)) +
                1;  // Número de regiões 			      -	 +1 para
                    // a condição de parada incluir todos
            long cpu_pos_count = 1 << (cpu_params.cpu_region - 1);
            // Número de posições na região 	-	 -1
            // porque são duas posições por iteração

            long cpu_ext_proj_id = 0;
            long inc_ext_proj_id =
                ~(cpu_reg_mask ^ global_reg_mask) & ((1 << qubits) - 1);

            long proj_id;  // indentificador local da região
            proj_id = cpu_ext_proj_id | cpu_proj_id;
            cpu_ext_proj_id =
                (cpu_ext_proj_id + inc_ext_proj_id + 1) & ~inc_ext_proj_id;
            cpu_reg_count--;

            while (proj_id != -1) {
              // Computa os operadores
              PCpuExecution1_0(state, pts, qubits, cpu_start, cpu_end,
                               cpu_pos_count, proj_id, cpu_reg_mask);

              proj_id = cpu_ext_proj_id | cpu_proj_id;
              cpu_ext_proj_id =
                  (cpu_ext_proj_id + inc_ext_proj_id + 1) & ~inc_ext_proj_id;
              cpu_reg_count--;
              if (cpu_reg_count <= 0) proj_id = -1;
            }

            cpu_start = cpu_end;
          }

#pragma omp critical(global_teste)
          {
            cpu_proj_id = ext_proj_id;
            ext_proj_id =
                (ext_proj_id + global_reg_mask + 1) & ~global_reg_mask;
            global_reg_count--;
            if (global_reg_count <= 0) cpu_proj_id = -1;
          }
        }

      }
      // #pragma omp section          //GPU EXECUTION
      else {
        long gpu_proj_id;

#pragma omp critical(global_teste)
        {
          gpu_proj_id = ext_proj_id;
          ext_proj_id = (ext_proj_id + global_reg_mask + 1) & ~global_reg_mask;
          global_reg_count--;
          if (global_reg_count <= 0) gpu_proj_id = -1;
        }

        while (gpu_proj_id != -1) {
          // Project Gates
          std::vector<PT *> gpu_pts;

          int gpu_i;

          int map_qb[qubits];
          memset(map_qb, -1, qubits * sizeof(int));

          int m = 0;
          for (gpu_i = 0; gpu_i < qubits; gpu_i++) {
            if ((1 << gpu_i) & global_reg_mask) {
              map_qb[gpu_i] = m++;
            }
          }

          PT *aux;
          gpu_pts.clear();
          for (int gpu_i = global_start; gpu_i < global_end; gpu_i++) {
            // verifica se o controle do operador satisfaz a parte global da
            // região
            if ((pts[gpu_i]->ctrl_mask & gpu_proj_id & ~global_reg_mask) ==
                (pts[gpu_i]->ctrl_value & ~global_reg_mask)) {
              aux = new PT();

              aux->qubits = pts[gpu_i]->qubits;

              aux->matrix = pts[gpu_i]->matrix;
              aux->mat_size = pts[gpu_i]->mat_size;
              aux->ctrl_mask = pts[gpu_i]->ctrl_mask & global_reg_mask;
              aux->ctrl_value = pts[gpu_i]->ctrl_value & global_reg_mask;

              aux->end = map_qb[pts[gpu_i]->end];
              aux->start = aux->end - log2((float)aux->mat_size);

              aux->ctrl_count = 0;
              for (int c = global_coales; c < qubits; c++) {
                if (aux->ctrl_mask & (1 << c)) {
                  aux->ctrl_count++;

                  aux->ctrl_mask &= ~(
                      1
                      << c);  // retira da mascara o controle do qubit atual (c)
                  aux->ctrl_mask |= (1 << map_qb[c]);  // e coloca o qubit que
                                                       // ele mapeia (map_qb[c])

                  if (aux->ctrl_value &
                      (1 << c)) {  // se o valor do controle for zero faz a
                                   // mesma coisa para ctrl_value;
                    aux->ctrl_mask &= ~(1 << c);
                    aux->ctrl_mask |= (1 << map_qb[c]);
                  }
                }
              }

              gpu_pts.push_back(aux);
            }
          }
          gpu_pts.push_back(NULL);
          ////////////////

          ProjectState(state, qubits, global_region, gpu_proj_id,
                       global_reg_mask, gpu_params.multi_gpu);

          GpuExecutionWrapper(NULL, &gpu_pts[0], global_region,
                              gpu_params.gpu_coales, gpu_params.gpu_region,
                              gpu_params.multi_gpu, gpu_params.tam_block,
                              gpu_params.rept, 1);

          GetState(state, qubits, global_region, gpu_proj_id, global_reg_mask,
                   gpu_params.multi_gpu);

          for (int c = 0; c < gpu_pts.size() - 1; c++) {
            delete gpu_pts[c];
          }

#pragma omp critical(global_teste)
          {
            gpu_proj_id = ext_proj_id;
            ext_proj_id =
                (ext_proj_id + global_reg_mask + 1) & ~global_reg_mask;
            global_reg_count--;
            if (global_reg_count <= 0) gpu_proj_id = -1;
          }
        }
      }
      //}
    }
  }
}

void DGM::setCpuStructure(long cpu_region, long cpu_coales) {
  this->cpu_params.cpu_region = cpu_region;
  this->cpu_params.cpu_coales = cpu_coales;
}

void DGM::setGpuStructure(long gpu_region, long gpu_coales, int rept) {
  this->gpu_params.gpu_region = gpu_region;
  this->gpu_params.gpu_coales = gpu_coales;
  this->gpu_params.rept = rept;
  this->gpu_params.tam_block = 1 << gpu_region / 2 / rept;
}

// Coalescimento
void MPI_coalesc(float complex *state, int qubits, int proj_qubits, long reg_id,
                 long reg_mask, int world_size) {
  int qbs_coales = 0;
  for (int i = 0; i < qubits; i++) {
    if ((reg_mask >> i) & 1) {
      qbs_coales++;
    } else {
      break;
    }
  }

  int mem_portions = pow(2.0, proj_qubits - qbs_coales);
  int portion_size = 1 << qbs_coales;

  float malloc_size = (1 << proj_qubits) * sizeof(float complex);

  float complex *new_state =
      (float complex *)(malloc(sizeof(float complex) * pow(2, qubits)));
  int *chunk_sizes = (int *)(malloc(sizeof(int) * world_size));
  int *displ = (int *)(malloc(sizeof(int) * world_size));

  long inc = ~(reg_mask >> qbs_coales);

  long dest_pos, src_pos, base = 0;
  for (int d = 0; d < world_size; d++) {
    displ[d] = dest_pos;
    for (int b = mem_portions / world_size * d;
         b < mem_portions / world_size * (d + 1); b++) {
      src_pos = (base << qbs_coales) | reg_id;

      memcpy(new_state + dest_pos, state + src_pos,
             portion_size * sizeof(float complex));

      base = (base + inc + 1) & ~inc;
      dest_pos += portion_size;
    }
    chunk_sizes[d] = dest_pos - displ[d];
  }
}

// void DistributedExecution(float complex *state, PT **pts, int qubits,
//                           long n_threads, int coales, int region, int it) {}