#include "pcpu.h"

#include <omp.h>

#include <vector>

void PCpuExecution1(float complex *state, PT **pts, e_size qubits,
                    int n_threads, e_size coales, e_size region) {
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
    MaskNewRegion data = getMaskAndRegion(pts, coales, region, i);
    // Executa até o operador na posiçao 'i' (exclusive) nesta iteração
    end = i;
    // Número de regiões
    e_size reg_count = (1LL << (qubits - data.region));
    // Número de posições na região: -1 porque são duas posições por iteração
    e_size pos_count = 1LL << (data.region - 1LL);
    // contador 'global' do número de regiões já computadas
    e_size ext_reg_id = 0LL;

    // em reg_ids temos is ids das regiões a serem executadas em paralelo
    for (e_size j = 0LL; j < reg_count; j++) {
      reg_ids[j] = ext_reg_id;
      ext_reg_id = (ext_reg_id + data.reg_mask + 1LL) & ~data.reg_mask;
    }

    omp_set_num_threads(n_threads);

#pragma omp parallel for schedule(runtime)
    for (size_t j = 0; j < reg_count; j++) {
      PCpuExecution1_0(state, pts, qubits, start, end, pos_count, reg_ids[j],
                       data.reg_mask);
    }
  }
}

void handle_dense(e_size &pos, e_size &pos0, e_size &pos1,
                  const e_size pos_count, const e_size reg_id,
                  const e_size shift, const e_size inc, const e_size pos_mask,
                  float complex *state, const PT *QG) {
  for (e_size p = 0LL; p < pos_count; p++) {
    pos0 = pos | reg_id;
    pos1 = pos0 | shift;
    pos = (pos + inc) & pos_mask;
    float complex tmp =
        QG->matrix[2] * state[pos0] + QG->matrix[3] * state[pos1];
    state[pos0] = QG->matrix[0] * state[pos0] + QG->matrix[1] * state[pos1];
    state[pos1] = tmp;
  }
}

void handle_primary(e_size &pos, e_size &pos0, e_size &pos1,
                    const e_size pos_count, const e_size reg_id,
                    const e_size shift, const e_size inc, const e_size pos_mask,
                    float complex *state, const PT *QG) {
  for (e_size p = 0LL; p < pos_count; p++) {
    pos0 = pos | reg_id;
    pos1 = pos0 | shift;
    pos = (pos + inc) & pos_mask;
    float complex tmp = QG->matrix[3] * state[pos1];
    state[pos0] *= QG->matrix[0];
    state[pos1] = tmp;
  }
}

void handle_secondary(e_size &pos, e_size &pos0, e_size &pos1,
                      const e_size pos_count, const e_size reg_id,
                      const e_size shift, const e_size inc,
                      const e_size pos_mask, float complex *state,
                      const PT *QG) {
  for (e_size p = 0LL; p < pos_count; p++) {
    pos0 = pos | reg_id;
    pos1 = pos0 | shift;
    pos = (pos + inc) & pos_mask;
    float complex tmp = QG->matrix[2] * state[pos0];
    state[pos0] = QG->matrix[1] * state[pos1];
    state[pos1] = tmp;
  }
}

typedef struct CtrlInfo {
  e_size reg_id;
  e_size pos_count;
  e_size reg_mask;
} CtrlInfo;

CtrlInfo prepare_controlled(const e_size reg_id, const e_size reg_mask,
                            const e_size pos_count, const e_size qubits,
                            const PT *QG) {
  // Verifica se a parte 'global' do controle satisfaz a região (reg_id)
  // É preciso arrumar o reg_mask retirando os qubits de controle que
  // estão dentro da região e arrumar o reg_id para incluir o valor dos
  // controles
  // Esta operação inclui o valor dos controles locais no reg_id (funciona
  // pois os valores globais já deram match)
  CtrlInfo result;
  result.reg_id = reg_id | QG->ctrl_value;
  // Valor inicial da mascara da região com controle
  result.reg_mask = reg_mask;
  // Número inicial de posições a serem calculadas
  result.pos_count = pos_count;

  // percorre os qubits
  for (e_size i = 0LL, m = 1LL; i < qubits; i++, m = m << 1LL) {
    // se o qubit pertencer a região e for um controle:
    if (m & reg_mask & QG->ctrl_mask) {
      //	remove ele da região(reg_mask) (para não iterar sobre ele)
      result.reg_mask ^= m;
      //	diminui a quantidade de posições que é preciso calcular.
      result.pos_count /= 2LL;
    }
  }

  return result;
}

void PCpuExecution1_0(float complex *state, PT **pts, e_size qubits,
                      e_size start, e_size end, e_size pos_count, e_size reg_id,
                      e_size reg_mask) {
  const PT *QG;
  e_size pos0;
  e_size pos1;

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
          handle_dense(pos, pos0, pos1, pos_count, reg_id, shift, inc, pos_mask,
                       state, QG);
          break;
        case DIAG_PRI:
          handle_primary(pos, pos0, pos1, pos_count, reg_id, shift, inc,
                         pos_mask, state, QG);
          break;
        case DIAG_SEC:
          handle_secondary(pos, pos0, pos1, pos_count, reg_id, shift, inc,
                           pos_mask, state, QG);
          break;
        default:
          printf("Erro de Tipo\n");
      }
    }
    // Importante: reg_id é o identificador da região e corresponde ao valor dos
    // qubits externos à região de operação (reg_mask)
    else if ((QG->ctrl_mask & reg_id & ~reg_mask) ==
             (QG->ctrl_value & ~reg_mask)) {
      CtrlInfo ctrl =
          prepare_controlled(reg_id, reg_mask, pos_count, qubits, QG);

      // mascara da posição --- retira o 'shift' da reg_mask,
      // para o 'inc pular sobre' esse bit também
      pos_mask = ctrl.reg_mask & ~shift;
      inc = ~pos_mask + 1LL;
      switch (mt) {
        case DENSE:
          handle_dense(pos, pos0, pos1, ctrl.pos_count, ctrl.reg_id, shift, inc,
                       pos_mask, state, QG);
          break;
        case DIAG_PRI:
          handle_primary(pos, pos0, pos1, ctrl.pos_count, ctrl.reg_id, shift,
                         inc, pos_mask, state, QG);
          break;
        case DIAG_SEC:
          handle_secondary(pos, pos0, pos1, ctrl.pos_count, ctrl.reg_id, shift,
                           inc, pos_mask, state, QG);
          break;
        default:
          printf("Erro de Tipo\n");
      }
    }
  }
}