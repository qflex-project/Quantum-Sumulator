#include "pt.h"

#include <math.h>

void PT::destructor() {
  if ((mat_size != 1) && !matrix) {
    delete matrix;
  }
  if (!ctrl_pos) {
    delete ctrl_pos;
  }
  if (!ctrl_rest) {
    delete ctrl_rest;
  }
}

long PT::ctrlAffect(long qubit) const {
  for (long i = 0; i < ctrl_count; i++)
    if (ctrl_pos[i] < qubit) {
      return i;
    }

  return ctrl_count;
}

long PT::matrixType() const {
  if ((matrix[1] == 0.0) && (matrix[2] == 0.0)) {
    return DIAG_PRI;
  } else if ((matrix[0] == 0.0) && (matrix[3] == 0.0)) {
    return DIAG_SEC;
  }

  return DENSE;
}

void PT::setArgs(long *arg, long /*affect*/) const {
  arg[SHIFT] = end;
  arg[CTRL_MASK] = ctrl_mask;
  arg[CTRL_VALUE] = ctrl_value;
}

void PT::setArgs_soft(long *arg, long /*affect*/) const {
  arg[SHIFT] = end;
  arg[CTRL_MASK] = ctrl_mask;
  arg[CTRL_VALUE] = ctrl_value;
}

void PT::setArgsGPU(long *arg, int region_start, int region_size,
                    int coalesc) const {
  long mask_coalesc = (1 << coalesc) - 1;
  long mask_region = ((1 << (region_size - coalesc)) - 1) << region_start;
  long mask_rest = ~(mask_coalesc | mask_region);

  arg[SHIFT] = end;
  if (end >= coalesc) arg[SHIFT] = end - region_start + coalesc;

  if (ctrl_mask) {
    arg[CTRL_MASK] = ctrl_mask & mask_rest;
    arg[CTRL_VALUE] = ctrl_value & mask_rest;

    arg[CTRL_REG_MASK] =
        (ctrl_mask & mask_coalesc) |
        ((ctrl_mask & mask_region) >> (region_start - coalesc));

    arg[CTRL_REG_VALUE] =
        (ctrl_value & mask_coalesc) |
        ((ctrl_value & mask_region) >> (region_start - coalesc));
  } else {
    arg[CTRL_MASK] = arg[CTRL_VALUE] = arg[CTRL_REG_MASK] =
        arg[CTRL_REG_VALUE] = 0;
  }
}

void PT::print() const {
  printf(
      "qubits: %d\nmat_size: %d\nstart: %d\nend: %d\nAffect: %d\nCtrl-(value: "
      "%ld, mask: %ld, count: %ld)\n",
      qubits, mat_size, start, end, affected, ctrl_value, ctrl_mask,
      ctrl_count);

  for (int i = 0; i < ctrl_count; i++) {
    printf("%d: %ld\n", i, ctrl_pos[i]);
  }
}

void PT::printMatrix() const {
  for (int i = 0; i < mat_size; i++) {
    for (int j = 0; j < mat_size; j++) {
      printf("%d: %.4f, %.4f  \t", i * mat_size + j,
             crealf(matrix[i * mat_size + j]),
             cimagf(matrix[i * mat_size + j]));
    }
    printf("\n");
  }
}

bool increasing(const PT *pt1, const PT *pt2) {
  if (pt1->affected == pt2->affected)
    return pt1->start > pt2->start;
  else
    return pt2->affected;
}

bool decreasing(const PT *pt1, const PT *pt2) {
  if (pt1->affected == pt2->affected)
    return pt1->start < pt2->start;
  else
    return pt1->affected;
}

MaskNewRegion getMaskAndRegion(PT **pts, e_size coales, e_size region,
                               e_size &i) {
  e_size count = coales;
  e_size reg_mask = coales ? (1LL << coales) - 1LL : 0LL;
  // Repete enquanto o número de qubits da região não
  // atingir o limite (region) e houver operadores
  while (count < region && pts[i] != NULL) {
    // Se o qubit do operador estiver fora da região (reg_mask),
    // incrementa o contador de qubits da região
    if (!((reg_mask >> pts[i]->end) & 1LL)) {
      count++;
    }
    // Acrescenta o qubit do operador na região se
    // ainda não tiver atingido o limite (region)
    if (count <= region) {
      reg_mask = reg_mask | (1LL << pts[i]->end);
    }
    i++;
  }
  // Segue acerscentado até encontrar um operador que não esteja dentro da
  // região
  while (pts[i] != NULL) {
    if ((reg_mask >> pts[i]->end) & 1LL) {
      i++;
    } else {
      break;
    }
  }
  // Se o número de qubits na região (count) não tiver atingido o limite
  // (region), acrescenta os ultimos qubits (final da mascara) à região até
  // completar
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