#include "lib_shor.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "common.h"
#include "dgm.h"
#include "gates.h"

int revert_bits(int res, int n) {
  int c = 0;

  for (int i = 0; i < n; i++) {
    c = (c << 1) | (res & 1);
    res = res >> 1;
  }
  return c;
}

int quantum_ipow(int a, int b) {
  int r = 1;

  for (int i = 0; i < b; i++) {
    r *= a;
  }

  return r;
}

/* Calculate the greatest common divisor with Euclid's algorithm */

int quantum_gcd(int u, int v) {
  int r;

  while (v) {
    r = v;
    v = u % v;
    u = r;
  }
  return u;
}

void quantum_frac_approx(int *a, int *b, int width) {
  float f = (float)*a / (float)*b;
  float g = f;
  int i;
  int num2 = 0;
  int den2 = 1;
  int num1 = 1;
  int den1 = 0;
  int num = 0;
  int den = 0;

  do {
    i = (int)(g + 0.000005);

    g -= (float)i - 0.000005f;
    g = 1.0f / g;

    if (i * den1 + den2 > 1 << width) break;

    num = i * num1 + num2;
    den = i * den1 + den2;

    num2 = num1;
    den2 = den1;
    num1 = num;
    den1 = den;

  } while (fabs(((double)num / den) - f) > 1.0 / (2 * (1 << width)));

  *a = num;
  *b = den;

  return;
}

void ApplyQFT(int qubits, int type, int multi_gpu, int qbs_region, int coalesc,
              int tam_block, int rept) {
  DGM dgm;
  dgm.exec_type = type;
  dgm.gpu_params.multi_gpu = multi_gpu;

  dgm.cpu_params.cpu_region = qbs_region;
  dgm.cpu_params.cpu_coales = coalesc;
  dgm.gpu_params.tam_block = tam_block;
  dgm.gpu_params.rept = rept;

  dgm.qubits = qubits;
  dgm.allocateMemory();
  dgm.setMemoryValue(0);

  std::vector<std::string> qft = QFT2(qubits, 0, qubits);

  dgm.executeFunction(qft);
}

std::string genRot(int qubits, int reg, int value) {
  std::vector<std::string> func(qubits, "ID");
  std::string name;

  int k = 2;
  float complex rot;
  float complex eps;
  eps = M_E;

  rot = 1;
  while (value) {
    if (value & 1) rot *= cpowf(eps, -2.0f * M_PI * I / (float)pow(2.0, k));
    value = value >> 1;
    k++;
  }

  if (rot != 1) {
    Gates g;
    name = "Rot_" + int2str(value);
    g.addGate(name, 1.0, 0.0, 0.0, rot);
    func[reg] = name;

    return concatena(func, qubits);
  }

  return "";
}

std::vector<std::string> CMultMod(int qubits, int ctrl, int reg1, int reg2,
                                  int over, int over_bool, int width, long a,
                                  long N) {
  int ctrl2;
  std::vector<std::string> qft = QFT(qubits, reg2, over, width);

  std::string HN = Hadamard(qubits, reg2, width);

  std::vector<std::string> rqft = RQFT(qubits, reg2, over, width);

  //////////////////////////////////////////////////////////////

  std::vector<std::string> mult_mod;
  std::vector<std::string> am;
  mult_mod.push_back(Hadamard(qubits, over, 1));
  mult_mod.push_back(HN);

  ctrl2 = reg1 + width - 1;
  for (int i = 0; i < width; i++) {
    am = C2AddMod(qubits, ctrl, ctrl2 - i, reg2, over, over_bool, width, a, N);
    mult_mod.insert(mult_mod.end(), am.begin(), am.end());
    a = (a * 2) % N;
  }

  mult_mod.insert(mult_mod.end(), rqft.begin(), rqft.end());

  return mult_mod;
}

std::vector<std::string> CRMultMod(int qubits, int ctrl, int reg1, int reg2,
                                   int over, int over_bool, int width, long a,
                                   long N) {
  int ctrl2;
  std::vector<std::string> qft = QFT(qubits, reg2, over, width);
  std::vector<std::string> rqft = RQFT(qubits, reg2, over, width);

  //////////////////////////////////////////////////////////////

  std::vector<std::string> mult_mod;
  std::vector<std::string> am;

  ctrl2 = reg1 + width - 1;
  for (int i = 0; i < width; i++) {
    am = C2SubMod(qubits, ctrl, ctrl2 - i, reg2, over, over_bool, width, a, N);
    mult_mod.insert(mult_mod.begin(), am.begin(), am.end());

    a = (a * 2) % N;
  }

  mult_mod.insert(mult_mod.begin(), qft.begin(), qft.end());
  mult_mod.insert(mult_mod.end(), rqft.begin(), rqft.end());

  return mult_mod;
}

std::vector<std::string> C2AddMod(int qubits, int ctrl1, int ctrl2, int reg,
                                  int over, int over_bool, int width, long a,
                                  long N) {
  std::vector<std::string> qft = QFT(qubits, reg, over, width);
  std::vector<std::string> rqft = RQFT(qubits, reg, over, width);

  std::string c2_add_a = C2AddF(qubits, ctrl1, ctrl2, reg, over, a, width);
  std::string c2_sub_a = C2SubF(qubits, ctrl1, ctrl2, reg, over, a, width);

  std::string sub_N = SubF(qubits, reg, over, N, width);
  std::string c_add_N = CAddF(qubits, over_bool, reg, over, N, width);

  std::string n_over = Pauli_X(qubits, over, 1);
  std::string c_over = CNot(qubits, over, over_bool);

  std::vector<std::string> func;

  func.push_back(c2_add_a);
  func.push_back(sub_N);
  func.insert(func.end(), rqft.begin(), rqft.end());
  func.push_back(c_over);
  func.insert(func.end(), qft.begin(), qft.end());
  func.push_back(c_add_N);
  func.push_back(c2_sub_a);
  func.insert(func.end(), rqft.begin(), rqft.end());
  func.push_back(n_over);
  func.push_back(c_over);
  func.push_back(n_over);
  func.insert(func.end(), qft.begin(), qft.end());
  func.push_back(c2_add_a);

  return func;
}

std::vector<std::string> C2SubMod(int qubits, int ctrl1, int ctrl2, int reg,
                                  int over, int over_bool, int width, long a,
                                  long N) {
  std::vector<std::string> qft = QFT(qubits, reg, over, width);
  std::vector<std::string> rqft = RQFT(qubits, reg, over, width);

  std::string c2_add_a = C2AddF(qubits, ctrl1, ctrl2, reg, over, a, width);
  std::string c2_sub_a = C2SubF(qubits, ctrl1, ctrl2, reg, over, a, width);

  std::string add_N = AddF(qubits, reg, over, N, width);
  std::string c_add_N = CAddF(qubits, over_bool, reg, over, N, width);
  std::string c_sub_N = CSubF(qubits, over_bool, reg, over, N, width);

  std::string n_over = Pauli_X(qubits, over, 1);
  std::string c_over = CNot(qubits, over, over_bool);

  std::vector<std::string> func;

  func.push_back(c2_sub_a);
  func.insert(func.end(), rqft.begin(), rqft.end());
  func.push_back(n_over);
  func.push_back(c_over);
  func.push_back(n_over);
  func.insert(func.end(), qft.begin(), qft.end());
  func.push_back(c2_add_a);
  func.push_back(c_sub_N);
  func.insert(func.end(), rqft.begin(), rqft.end());
  func.push_back(c_over);
  func.insert(func.end(), qft.begin(), qft.end());
  func.push_back(add_N);
  func.push_back(c2_sub_a);

  return func;
}

//////////////////////////////////////////////////////////////////////////

std::string CAddF(int qubits, int ctrl1, int reg, int over, long num,
                  int width) {
  std::vector<std::string> caf = AddF(qubits, reg, over, num, width, true);

  caf[ctrl1] = "Control1(1)";

  return concatena(caf, qubits);
}

std::string C2AddF(int qubits, int ctrl1, int ctrl2, int reg, int over,
                   long num, int width) {
  std::vector<std::string> caf = AddF(qubits, reg, over, num, width, true);

  caf[ctrl1] = "Control1(1)";
  caf[ctrl2] = "Control1(1)";

  return concatena(caf, qubits);
}

std::string AddF(int qubits, int reg, int over, long num, int width) {
  return concatena(AddF(qubits, reg, over, num, width, false), qubits);
}

std::vector<std::string> AddF(int qubits, int reg, int over, long num,
                              int width, bool controlled) {
  int size = width + 1;
  std::vector<float complex> rot(size, 1);
  float complex c;

  Gates g;

  long aux = num;

  float complex eps = M_E;

  for (int i = 0; i < size; i++) {
    if (aux & 1)
      for (int j = i; j < size; j++)
        rot[j] *= cpowf(eps, 2 * M_PI * I / pow(2.0, j - i + 1));
    aux = aux >> 1;
  }

  std::vector<std::string> add(qubits, "ID");
  std::string name;

  aux = reg + width - 1;
  c = 1;
  for (int i = 0; i < size; i++) {
    if (rot[i] != c) {
      name = "ADD_" + int2str(num) + "_" + int2str(i);
      g.addGate(name, 1.0, 0.0, 0.0, rot[i]);
      if (controlled) name = "Target1(" + name + ")";
      add[aux - i] = name;
    }
  }

  name = add[reg - 1];
  add[reg - 1] = "ID";
  add[over] = name;

  return add;
}

std::string CSubF(int qubits, int ctrl1, int reg, int over, long num,
                  int width) {
  std::vector<std::string> csf = SubF(qubits, reg, over, num, width, true);

  csf[ctrl1] = "Control1(1)";

  return concatena(csf, qubits);
}

std::string C2SubF(int qubits, int ctrl1, int ctrl2, int reg, int over,
                   long num, int width) {
  std::vector<std::string> csf = SubF(qubits, reg, over, num, width, true);

  csf[ctrl1] = "Control1(1)";
  csf[ctrl2] = "Control1(1)";

  return concatena(csf, qubits);
}

std::string SubF(int qubits, int reg, int over, long num, int width) {
  return concatena(SubF(qubits, reg, over, num, width, false), qubits);
}

std::vector<std::string> SubF(int qubits, int reg, int over, long num,
                              int width, bool controlled) {
  long size = width + 1;
  std::vector<float complex> rot(size, 1);
  float complex c;

  Gates g;

  long aux = num;

  float complex eps = M_E;
  for (int i = 0; i < size; i++) {
    if (aux & 1)
      for (int j = i; j < size; j++)
        rot[j] *= cpowf(eps, -2 * M_PI * I / pow(2.0, j - i + 1));
    aux = aux >> 1;
  }

  std::vector<std::string> sub(qubits, "ID");
  std::string name;

  aux = reg + width - 1;
  c = 1;
  for (int i = 0; i < size; i++) {
    if (rot[i] != c) {
      name = "SUB_" + int2str(num) + "_" + int2str(i);
      g.addGate(name, 1.0, 0.0, 0.0, rot[i]);
      if (controlled) name = "Target1(" + name + ")";
      sub[aux - i] = name;
    }
  }

  name = sub[reg - 1];
  sub[reg - 1] = "ID";
  sub[over] = name;

  return sub;
}

std::vector<std::string> QFT(int qubits, int reg, int over, int width) {
  std::string s, name;
  std::vector<std::string> qft;

  Gates g;
  float complex c;
  for (int i = 1; i <= width + 1; i++) {
    name = "R" + int2str(i);

    c = M_E;
    c = cpowf(c, 2 * M_PI * I / pow(2.0, i));

    g.addGate(name, 1.0, 0.0, 0.0, c);
  }

  std::vector<std::string> base(qubits, "ID");

  qft.push_back(Hadamard(qubits, over, 1));
  for (int j = 0; j < width; j++) {
    base[j + reg] = "Control1(1)";
    base[over] = "Target1(R" + int2str(j + 2) + ")";

    s = concatena(base, qubits);
    qft.push_back(s);
    base[j + reg] = "ID";
  }
  base[over] = "ID";

  for (int i = 0; i < width; i++) {
    qft.push_back(Hadamard(qubits, i + reg, 1));

    for (int j = i + 1; j < width; j++) {
      base[j + reg] = "Control1(1)";
      base[i + reg] = "Target1(R" + int2str(j - i + 1) + ")";

      s = concatena(base, qubits);
      qft.push_back(s);

      base[j + reg] = "ID";
    }
    base[i + reg] = "ID";
  }

  return qft;
}

std::vector<std::string> QFT2(int qubits, int reg, int width) {
  std::string s;
  std::vector<std::string> qft;

  Gates g;
  float complex c;
  for (int i = 1; i <= width + 1; i++) {
    c = M_E;
    c = cpowf(c, 2 * M_PI * I / pow(2.0, i));
    g.addGate("R-" + int2str(i), 1.0, 0.0, 0.0, c);
  }

  std::vector<std::string> base(qubits, "ID");

  for (int i = 0; i < width; i++) {
    base[i + reg] = "H";
    s = concatena(base, qubits);
    qft.push_back(s);

    for (int j = i + 1; j < width; j++) {
      base[j + reg] = "Control1(1)";
      base[i + reg] = "Target1(R-" + int2str(j - i + 1) + ")";

      s = concatena(base, qubits);
      qft.push_back(s);

      base[j + reg] = "ID";
    }
    base[i + reg] = "ID";
  }

  return qft;
}

std::vector<std::string> RQFT(int qubits, int reg, int over, int width) {
  std::string s;
  std::vector<std::string> rqft;

  Gates g;
  float complex c;
  for (int i = 1; i <= width + 1; i++) {
    c = M_E;
    c = cpowf(c, -2 * M_PI * I / pow(2.0, i));
    g.addGate("R'" + int2str(i), 1.0, 0.0, 0.0, c);
  }

  std::vector<std::string> base(qubits, "ID");

  rqft.push_back(Hadamard(qubits, over, 1));
  for (int j = 0; j < width; j++) {
    base[j + reg] = "Control1(1)";
    base[over] = "Target1(R'" + int2str(j + 2) + ")";

    s = concatena(base, qubits);
    rqft.push_back(s);
    base[j + reg] = "ID";
  }
  base[over] = "ID";

  for (int i = 0; i < width; i++) {
    base[i + reg] = "H";
    s = concatena(base, qubits);
    rqft.push_back(s);

    for (int j = i + 1; j < width; j++) {
      base[j + reg] = "Control1(1)";
      base[i + reg] = "Target1(R'" + int2str(j - i + 1) + ")";

      s = concatena(base, qubits);
      rqft.push_back(s);

      base[j + reg] = "ID";
    }
    base[i + reg] = "ID";
  }
  reverse(rqft.begin(), rqft.end());

  return rqft;
}

std::vector<std::string> CSwapR(int qubits, int ctrl, int reg1, int reg2,
                                int width) {
  std::vector<std::string> sw;
  std::vector<std::string> base(qubits, "ID");
  std::string s1, s2;

  for (int i = 0; i < width; i++) {
    base[ctrl] = "Control1(1)";
    base[i + reg1] = "Target1(X)";
    base[i + reg2] = "Control1(1)";
    s1 = concatena(base, qubits);

    base[ctrl] = "ID";
    base[i + reg1] = "Control1(1)";
    base[i + reg2] = "Target1(X)";
    s2 = concatena(base, qubits);

    base[i + reg1] = base[i + reg2] = "ID";

    sw.push_back(s2);
    sw.push_back(s1);
    sw.push_back(s2);
  }

  return sw;
}

std::vector<std::string> SwapOver(int qubits, int reg, int width) {
  std::vector<std::string> so;

  for (int i = 0; i < width / 2; i++) {
    so.push_back(CNot(qubits, reg + width - i - 1, reg + i));
    so.push_back(CNot(qubits, reg + i, reg + width - i - 1));
    so.push_back(CNot(qubits, reg + width - i - 1, reg + i));
  }

  return so;
}

long mul_inv(long a, long b) {
  long b0 = b, t, q;
  long x0 = 0, x1 = 1;
  if (b == 1) return 1;
  while (a > 1) {
    q = a / b;
    t = b, b = a % b, a = t;
    t = x0, x0 = x1 - q * x0, x1 = t;
  }
  if (x1 < 0) x1 += b0;

  return x1;
}

std::string int2str(int number) {
  std::stringstream ss;
  ss << number;

  std::string str = ss.str();

  return str;
}

//////////////////////////////////////////////////////

// N - Number to ne factored
// type - Execution Type
// threads - Number of threads to be used in case of a parallel execution on CPU
std::vector<int> Shor(long N, int type, int n_threads, int cpu_region,
                      int cpu_coalesc, int multi_gpu, int gpu_region,
                      int gpu_coalesc, int tam_block, int rept) {
  long a, n, mod_a, mod_inv_a, aux, m, res;

  int qubits, qft_qb, reg1, reg2, over, over_bool;
  int f1, f2, factor;

  //////////////////////////////////
  DGM dgm;
  dgm.exec_type = type;
  dgm.gpu_params.multi_gpu = multi_gpu;
  dgm.cpu_params.n_threads = n_threads;

  dgm.cpu_params.cpu_region = cpu_region;
  dgm.cpu_params.cpu_coales = cpu_coalesc;

  dgm.gpu_params.gpu_region = gpu_region;
  dgm.gpu_params.gpu_coales = gpu_coalesc;
  dgm.gpu_params.tam_block = tam_block;
  dgm.gpu_params.rept = rept;
  //-----------------------------//

  aux = N;
  a = n = 0;
  while (aux) {
    n++;
    aux = aux >> 1;
  }
  qubits = 2 * n + 3;

  //////////////////////////////////////
  dgm.qubits = qubits;
  dgm.allocateMemory();
  dgm.setMemoryValue((1 << (n + 2)));
  //----------------------------------//

  qft_qb = 0;
  reg1 = 1;
  reg2 = n + 2;
  over = n + 1;
  over_bool = qubits - 1;

  while ((quantum_gcd(N, a) > 1) || (a < 2)) {
    a = rand() % N;
  }

  std::string X0 = Pauli_X(qubits, 0, 1);
  std::string H0 = Hadamard(qubits, qft_qb, 1);

  res = 0;
  int L = 2 * n - 1;
  long inv_a = mul_inv(a, N);

  std::vector<std::string> func, f;

  for (int i = L; i >= 0; i--) {
    mod_a = modular_pow(a, pow(2, i), N);
    mod_inv_a = modular_pow(inv_a, pow(2, i), N);

    func.clear();

    func.push_back(H0);

    f = CMultMod(qubits, qft_qb, reg1, reg2, over, over_bool, n, mod_a, N);

    func.insert(func.end(), f.begin(), f.end());
    f = CSwapR(qubits, qft_qb, reg1, reg2, n);
    func.insert(func.end(), f.begin(), f.end());

    f = CRMultMod(qubits, qft_qb, reg1, reg2, over, over_bool, n, mod_a, N);
    func.insert(func.end(), f.begin(), f.end());
    func.push_back(H0);

    if (res) func.push_back(genRot(qubits, qft_qb, res));

    dgm.executeFunction(func);

    m = dgm.measure(qft_qb);

    res = (res << 1) | m;
  }

  /*
          dgm.setFunction(func);
          //cout << "Aqui" << endl;
          dgm.CountOps();


          cout << "Shor " << qubits << " qubits" << endl;
          cout << "Dense: " << dgm.dense << endl;
          cout << "Main Diagonal: " << dgm.main_diag << endl;
          cout << "Secondary Diagonal: " << dgm.sec_diag << endl;
          cout << "C-Dense: " << dgm.c_dense << endl;
          cout << "C-Main Diagonal: " << dgm.c_main_diag << endl;
          cout << "C-Secondary Diagonal: " << dgm.c_sec_diag << endl;
          cout << "Total: " << dgm.total_op << endl << endl;

          return;
  */

  // Result check
  std::vector<int> factors;

  int c = revert_bits(res, 2 * n);

  if (c == 0) {
    return factors;
  }

  int q = 1 << (2 * n);

  quantum_frac_approx(&c, &q, n);

  int r = q;
  int i = 1;
  while ((r * i) < (1 << n)) {
    if (modular_pow(a, r * i, N) == 1) {
      q = r * i;
      break;
    }
    i++;
  }
  if (q >= N) q = r;

  /*
  if((q % 2 == 1) && (2*q<(1<<n)))
  {
          //printf("Odd denominator, trying to expand by 2.\n");
          q *= 2;
  }

  if(q % 2 == 1)
  {
          //printf("Odd period, try again.\n");
          return;
  }
  */

  i = modular_pow(a, q / 2, N);
  f1 = quantum_gcd(N, i + 1);
  f2 = quantum_gcd(N, i - 1);

  if (f1 > f2)
    factor = f1;
  else
    factor = f2;

  if ((factor < N) && (factor > 1)) {
    factors.push_back(factor);
    factors.push_back((int)N / factor);
    return factors;
  }

  if (r != q) {
    i = modular_pow(a, r / 2, N);
    f1 = quantum_gcd(N, i + 1);
    f2 = quantum_gcd(N, i - 1);

    if (f1 > f2)
      factor = f1;
    else
      factor = f2;

    if ((factor < N) && (factor > 1)) {
      factors.push_back(factor);
      factors.push_back((int)N / factor);
      return factors;
    }
  }

  return factors;
}