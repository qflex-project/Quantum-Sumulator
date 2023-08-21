#include "pcpu.h"

#include <omp.h>

void PCpuExecution1(float complex *state, PT **pts, int qubits, long n_threads,
                    int coales, int region, int it) {
  long i = 0;
  long start = 0;
  long end = 0;

  long max_reg_count = (1 << (qubits - region));
  long *reg_ids = (long *)malloc(max_reg_count * sizeof(long));

  while (pts[i] != NULL) {
    long count = coales;
    long reg_mask = coales ? (1 << coales) - 1 : 0;

    // Pega os operadores que estão dentro da região coalescida (reg_mask
    // inicial), e acrescenta operadores em qubits fora dela até chegar ao
    // limite da região (region definida)
    start = i;
    while (count < region &&
           pts[i] !=
               NULL) {  // Repete enquanto o número de qubits da região não
                        // atingir o limite (region) e houver operadores
      if (              // pts[i]->matrixType() != DIAG_PRI &&
            // //O qubit de operadores de diagonal principal não importa para
            // região (sempre podem ser acrescentados)
          !((reg_mask >> pts[i]->end) &
            1)) {  // Se o qubit do operador estiver fora da região (reg_mask),
                   // incrementa o contador de qubits da região
        count++;
      }

      if (count <= region)  // && pts[i]->matrixType() != DIAG_PRI)
        reg_mask =
            reg_mask |
            (1 << pts[i]->end);  // Acrescenta o qubit do operador na região se
                                 // ainda não tiver atingido o limite (region)

      i++;
    }
    // Segue acerscentado até encontrar um operador que não esteja dentro da
    // região
    while (pts[i] != NULL) {
      if ((reg_mask >> pts[i]->end) & 1) {
        i++;
      } else {
        break;
      }
    }
    end = i;  // Executa até o operador na posiçao 'i' (exclusive) nesta
              // iteração

    // Se o número de qubits na região (count) não tiver atingido o limite
    // (region), acrescenta os ultimos qubits (final da mascara) à região até
    // completar for (long a = 1<<(qubits-1); count < region; a = a >> 1){
    long a = 1;
    while (count < region) {
      if (a & ~reg_mask) {
        reg_mask = reg_mask | a;
        count++;
      }
      a = a << 1;
    }
    // for (long a = 1; count < region; a = a << 1) {
    //   if (a & ~reg_mask) {
    //     reg_mask = reg_mask | a;
    //     count++;
    //   }
    // }

    if (count < region) region = count;

    long reg_count = (1 << (qubits - region));  // Número de regiões
    long pos_count =
        1 << (region - 1);  // Número de posições na região: -1 porque são duas
                            // posições por iteração

    omp_set_num_threads(n_threads);

    long ext_reg_id = 0;  // contador 'global' do número de regiões já
                          // computadas

    // em reg_ids temos is ids das regiões a serem executadas em paralelo
    for (size_t j = 0; j < reg_count; j++) {
      reg_ids[j] = ext_reg_id;
      ext_reg_id = (ext_reg_id + reg_mask + 1) & ~reg_mask;
    }

#pragma omp parallel for schedule(runtime)
    for (size_t j = 0; j < reg_count; j++) {
      PCpuExecution1_0(state, pts, qubits, start, end, pos_count, reg_ids[j],
                       reg_mask);
    }
  }
  free(reg_ids);
}

void PCpuExecution1_0(float complex *state, PT **pts, int qubits, int start,
                      int end, int pos_count, int reg_id, int reg_mask) {
  PT *QG;
  long pos0, pos1;
  float complex tmp;

  for (int op = start; op < end; op++) {
    QG = pts[op];
    long shift = (1 << QG->end);  // mascara com a posição do qubit do operador
    long mt = QG->matrixType();
    // if (mt == DIAG_PRI) shift = coalesc;	//se for um operador de diagonal
    // principal, a posição do qubit não é relevante
    long pos_mask =
        reg_mask &
        ~shift;  // mascara da posição --- retira o 'shift' da reg_mask, para o
                 // 'inc pular sobre ' esse bit também
    long inc =
        ~pos_mask + 1;  // usado para calcular a proxima posição de uma região
    long pos = 0;

    if (!QG->ctrl_count) {
      switch (mt) {
        case DENSE:
          for (long p = 0; p < pos_count; p++) {
            pos0 = pos | reg_id;
            pos1 = pos0 | shift;

            pos = (pos + inc) & pos_mask;

            tmp = QG->matrix[2] * state[pos0] + QG->matrix[3] * state[pos1];
            state[pos0] =
                QG->matrix[0] * state[pos0] + QG->matrix[1] * state[pos1];
            state[pos1] = tmp;
            // considerando a portatibilidade para uso de MPI:
            // pos0 e pos1 são alterados pos_count vezes
            // cada processo terá um state diferente após a execução
            // nas posições que pos0 e pos1 indicaram durante a execução
            // como juntar os states de cada processo em um state único?
            // pois atualmente o state é compartilhado
            // o ideal seria enviar para cada processo exatamente as posições
            // que serão operadas e com isso usar o scatter e o gather para
            // unificar (o gatherAll talvez) o state
          }
          break;
        case DIAG_PRI:
          for (long p = 0; p < pos_count; p++) {
            pos0 = pos | reg_id;
            pos1 = pos0 | shift;

            pos = (pos + inc) & pos_mask;

            tmp = QG->matrix[3] * state[pos1];
            state[pos0] *= QG->matrix[0];  // * state[pos0];
            state[pos1] = tmp;             // * state[pos1];tmp;
          }
          break;

        case DIAG_SEC:
          for (long p = 0; p < pos_count; p++) {
            pos0 = pos | reg_id;
            pos1 = pos0 | shift;

            pos = (pos + inc) & pos_mask;

            tmp = QG->matrix[2] * state[pos0];
            state[pos0] = QG->matrix[1] * state[pos1];
            state[pos1] = tmp;
          }
          break;
        default:
          printf("Erro de Tipo\n");
      }
    }
    // Importante: reg_id é o identificador da região e corresponde ao valor dos
    // qubits externos à região de operação (reg_mask)
    else {
      if ((QG->ctrl_mask & reg_id & ~reg_mask) ==
          (QG->ctrl_value &
           ~reg_mask)) {  // Verifica se a parte 'global' do controle satisfaz a
                          // região (reg_id)

        // É preciso arrumar o reg_mask retirando os qubits de controle que
        // estão dentro da região e arrumar o reg_id para incluir o valor dos
        // controles
        long ctrl_reg_id =
            reg_id | QG->ctrl_value;  // Esta operação inclui o valor dos
                                      // controles locais no reg_id (funciona
                                      // pois os valores globais já deram match)
        long ctrl_reg_mask =
            reg_mask;  // Valor inicial da mascara da região com controle
        long ctrl_pos_count =
            pos_count;  // Número inicial de posições a serem calculadas

        for (int i = 0, m = 1; i < qubits;
             i++, m = m << 1) {                // percorre os qubits
          if (m & reg_mask & QG->ctrl_mask) {  // se o qubit pertencer a região
                                               // e for um controle:
            ctrl_reg_mask ^= m;  //	remove ele da região(reg_mask) (para não
                                 // iterar sobre ele)
            ctrl_pos_count /=
                2;  //	diminui a quantidade de posições que é preciso calcular.
          }
        }

        pos_mask =
            ctrl_reg_mask &
            ~shift;  // mascara da posição --- retira o 'shift' da reg_mask,
                     // para o 'inc pular sobre' esse bit também
        inc = ~pos_mask + 1;

        switch (mt) {
          case DENSE:
            for (long p = 0; p < ctrl_pos_count; p++) {
              pos0 = pos | ctrl_reg_id;
              pos1 = pos0 | shift;

              pos = (pos + inc) & pos_mask;

              tmp = QG->matrix[2] * state[pos0] + QG->matrix[3] * state[pos1];
              state[pos0] =
                  QG->matrix[0] * state[pos0] + QG->matrix[1] * state[pos1];
              state[pos1] = tmp;
            }
            break;
          case DIAG_PRI:
            for (long p = 0; p < ctrl_pos_count; p++) {
              pos0 = pos | ctrl_reg_id;
              pos1 = pos0 | shift;

              pos = (pos + inc) & pos_mask;

              tmp = QG->matrix[3] * state[pos1];
              state[pos0] *= QG->matrix[0];
              state[pos1] = tmp;
            }
            break;

          case DIAG_SEC:
            for (long p = 0; p < ctrl_pos_count; p++) {
              pos0 = pos | ctrl_reg_id;
              pos1 = pos0 | shift;

              pos = (pos + inc) & pos_mask;

              tmp = QG->matrix[2] * state[pos0];
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