#include "pcpu.h"

#include <vector>
#include <omp.h>

typedef struct MaskAndRegion {
  e_size reg_mask;
  e_size region;
} MaskNewRegion;

MaskNewRegion getMaskAndRegion(PT **pts, e_size coales, e_size region, e_size &i) {
  e_size count = coales;
  e_size reg_mask = coales ? (1LL << coales) - 1LL : 0LL;
  // Repete enquanto o número de qubits da região não
  // atingir o limite (region) e houver operadores
  while (count < region && pts[i] !=NULL) {
    // Se o qubit do operador estiver fora da região (reg_mask),
    // incrementa o contador de qubits da região
    if (!((reg_mask >> pts[i]->end) & 1LL)) {              
      count++;
    }
    // Acrescenta o qubit do operador na região se
    // ainda não tiver atingido o limite (region)
    if (count <= region){
      reg_mask = reg_mask | (1LL << pts[i]->end);
    }
    i++;
  }
  // Segue acerscentado até encontrar um operador que não esteja dentro da região
  while (pts[i] != NULL) {
    if ((reg_mask >> pts[i]->end) & 1LL) {
      i++;
    } else {
      break;
    }
  }
  // Se o número de qubits na região (count) não tiver atingido o limite (region),
  // acrescenta os ultimos qubits (final da mascara) à região até completar
  e_size a = 1LL;
  while (count < region) {
    if (a & ~reg_mask) {
      reg_mask = reg_mask | a;
      count++;
    }
    a = a << 1LL;
  }

  if (count < region) region = count;

  MaskNewRegion result;
  result.reg_mask = reg_mask;
  result.region = region;
  return result;
}

void PCpuExecution1(float complex *state, PT **pts, e_size qubits, int n_threads,
                    e_size coales, e_size region) {
  e_size i = 0LL;
  e_size start = 0LL;
  e_size end = 0LL;

  e_size max_reg_count = (1LL << (qubits - region));
  std::vector<e_size> reg_ids(max_reg_count);
  // Pega os operadores que estão dentro da região coalescida (reg_mask
  // inicial), e acrescenta operadores em qubits fora dela até chegar ao
  // limite da região (region definida)
  while (pts[i] != NULL) {
    start = i;
    MaskNewRegion maskNewRegion = getMaskAndRegion(pts, coales, region, i);
    // Executa até o operador na posiçao 'i' (exclusive) nesta iteração
    end = i;  
    // Número de regiões
    e_size reg_count = (1LL << (qubits - maskNewRegion.region));  
    // Número de posições na região: -1 porque são duas posições por iteração
    e_size pos_count = 1LL << (maskNewRegion.region - 1LL);
    // contador 'global' do número de regiões já computadas
    e_size ext_reg_id = 0LL;  

    // em reg_ids temos is ids das regiões a serem executadas em paralelo
    for (e_size j = 0LL; j < reg_count; j++) {
      reg_ids[j] = ext_reg_id;
      ext_reg_id = (ext_reg_id + maskNewRegion.reg_mask + 1LL) & ~maskNewRegion.reg_mask;
    }

    omp_set_num_threads(n_threads);

#pragma omp parallel for schedule(runtime)
    for (size_t j = 0; j < reg_count; j++) {
      PCpuExecution1_0(state, pts, qubits, start, end, pos_count, reg_ids[j],
                       maskNewRegion.reg_mask);
    }
  }
}

void PCpuExecution1_0(float complex *state, PT **pts, e_size qubits, e_size start,
                      e_size end, e_size pos_count, e_size reg_id, e_size reg_mask) {
  PT *QG;
  e_size pos0;
  e_size pos1;
  float complex tmp;

  for (e_size op = start; op < end; op++) {
    QG = pts[op];
    // mascara com a posição do qubit do operador
    e_size shift = (1LL << QG->end);  
    e_size mt = QG->matrixType();
    // mascara da posição --- retira o 'shift' da reg_mask, para o
    // 'inc pular sobre ' esse bit também
    e_size pos_mask = reg_mask & ~shift;
    // usado para calcular a proxima posição de uma região
    e_size inc = ~pos_mask + 1LL;  
    e_size pos = 0LL;

    if (!QG->ctrl_count) {
      switch (mt) {
        case DENSE:
          for (e_size p = 0LL; p < pos_count; p++) {
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
          for (e_size p = 0LL; p < pos_count; p++) {
            pos0 = pos | reg_id;
            pos1 = pos0 | shift;

            pos = (pos + inc) & pos_mask;

            tmp = QG->matrix[3] * state[pos1];
            state[pos0] *= QG->matrix[0];
            state[pos1] = tmp;
          }
          break;

        case DIAG_SEC:
          for (e_size p = 0LL; p < pos_count; p++) {
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
      // Verifica se a parte 'global' do controle satisfaz a região (reg_id)
      // É preciso arrumar o reg_mask retirando os qubits de controle que
      // estão dentro da região e arrumar o reg_id para incluir o valor dos
      // controles
      if ((QG->ctrl_mask & reg_id & ~reg_mask) == (QG->ctrl_value & ~reg_mask)) {
        // Esta operação inclui o valor dos controles locais no reg_id (funciona
        // pois os valores globais já deram match)
        e_size ctrl_reg_id = reg_id | QG->ctrl_value;
        // Valor inicial da mascara da região com controle
        e_size ctrl_reg_mask = reg_mask;  
        // Número inicial de posições a serem calculadas
        e_size ctrl_pos_count = pos_count;  

        // percorre os qubits
        for (e_size i = 0LL, m = 1LL; i < qubits; i++, m = m << 1LL) {
          // se o qubit pertencer a região e for um controle:
          if (m & reg_mask & QG->ctrl_mask) {  
            //	remove ele da região(reg_mask) (para não iterar sobre ele)
            ctrl_reg_mask ^= m;  
            //	diminui a quantidade de posições que é preciso calcular.
            ctrl_pos_count /= 2LL;  
          }
        }
        // mascara da posição --- retira o 'shift' da reg_mask, 
        // para o 'inc pular sobre' esse bit também
        pos_mask = ctrl_reg_mask & ~shift;  
        inc = ~pos_mask + 1LL;

        switch (mt) {
          case DENSE:
            for (e_size p = 0LL; p < ctrl_pos_count; p++) {
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
            for (e_size p = 0LL; p < ctrl_pos_count; p++) {
              pos0 = pos | ctrl_reg_id;
              pos1 = pos0 | shift;

              pos = (pos + inc) & pos_mask;

              tmp = QG->matrix[3] * state[pos1];
              state[pos0] *= QG->matrix[0];
              state[pos1] = tmp;
            }
            break;

          case DIAG_SEC:
            for (e_size p = 0LL; p < ctrl_pos_count; p++) {
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