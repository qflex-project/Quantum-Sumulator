#include "cpu.h"

#include <math.h>
#include <omp.h>
#include <vector>

#include <string.h>

#include "pcpu.h"
#include "gpu.h"


void HybridExecution(float complex *state, PT **pts, int qubits, const CPUParams& cpu_params, const GPUParams& gpu_params) {
  int mem_size = (int)pow(2.0, qubits);
  int qubits_limit = 20;
  int global_coales = 15;
  // (cpu_coales > gpu_coales) ? cpu_coales : gpu_coales;

  int global_region = qubits_limit;
  int global_start;
  int global_end;

  int global_count;
  int global_reg_mask;
  int global_reg_count;
  int ext_proj_id;
  //, global_pos_count; //atualmente não utilizada

  omp_set_num_threads(cpu_params.n_threads);

  int i = 0;
  while (pts[i] != NULL) {
    global_count = global_coales;
    global_reg_mask = global_coales ? (1 << global_coales) - 1 : 0;

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
      if ((global_reg_mask >> pts[i]->end) & 1)
           // || pts[i]->matrixType() == DIAG_PRI)
        i++;
      else
        break;
    }
    global_end = i;

    // Se o número de qubits na região (count) nãoo tiver atingido o limite
    // (region), acrescenta os ultimos qubits (final da mascara) à região até
    // completar for (long a = 1<<(qubits-1); count < region; a = a >> 1){
    for (int a = 1; global_count < global_region; a = a << 1) {
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

    #pragma omp parallel num_threads(cpu_params.n_threads)
    {
      if (omp_get_thread_num() != 0) {  // CPU EXECUTION
        int cpu_proj_id;

        #pragma omp critical(global_teste)
        {
          cpu_proj_id = ext_proj_id;
          ext_proj_id = (ext_proj_id + global_reg_mask + 1) & ~global_reg_mask;
          global_reg_count--;
          if (global_reg_count <= 0) cpu_proj_id = -1;
        }

        while (cpu_proj_id != -1) {
          int cpu_i;
          int cpu_start;
          int cpu_end;

          cpu_start = global_start;

          cpu_i = cpu_start;

          while (cpu_start < global_end) {
            int cpu_count = cpu_params.cpu_coales;
            int cpu_reg_mask =
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
              if ((cpu_reg_mask >> pts[cpu_i]->end) & 1)
                // || pts[i]->matrixType() == DIAG_PRI)
                cpu_i++;
              else
                break;
            }
            cpu_end = cpu_i;

            for (int a = 1; cpu_count < cpu_params.cpu_region; a = a << 1) {
              if ((a & global_reg_mask) &&
                  (a & ~cpu_reg_mask)) {  // tem que não estar na região da cpu
                                          // e estar na global
                cpu_reg_mask = cpu_reg_mask | a;
                cpu_count++;
              }
            }

            int cpu_reg_count =
                (1 << (global_region - cpu_params.cpu_region)) +
                1;  // Número de regiões 			      -	 +1 para
                    // a condição de parada incluir todos
            int cpu_pos_count = 1 << (cpu_params.cpu_region - 1);
            // Número de posições na região 	-	 -1
            // porque são duas posições por iteração

            int cpu_ext_proj_id = 0;
            int inc_ext_proj_id =
                ~(cpu_reg_mask ^ global_reg_mask) & ((1 << qubits) - 1);

            int proj_id;  // indentificador local da região
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
            if (global_reg_count <= 0) {cpu_proj_id = -1;}
          }
        }

      }
      // #pragma omp section          //GPU EXECUTION
      else {
        int gpu_proj_id;

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