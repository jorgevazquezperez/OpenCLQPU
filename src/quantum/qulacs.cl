
//! unsigned integer
typedef unsigned int UINT;

// ---------------------------------- DEFINITION OF COMPLEX
// NUMBER------------------------------------------
typedef struct {
  double real;
  double imag;
} Complex;

Complex double_to_complex(double real, double imag) {
  Complex c;

  c.real = real;
  c.imag = imag;

  return c;
}

double real(Complex c) { return c.real; }

double imag(Complex c) { return c.imag; }

Complex conj(Complex c) {

  c.imag = -c.imag;

  return c;
}

Complex add(Complex c1, Complex c2) {
  Complex c3;

  c3.real = c1.real + c2.real;
  c3.imag = c1.imag + c2.imag;

  return c3;
}

Complex sub(Complex c1, Complex c2) {
  Complex c3;

  c3.real = c1.real - c2.real;
  c3.imag = c1.imag - c2.imag;

  return c3;
}
Complex mul(Complex c1, Complex c2) {
  Complex c3;

  c3.real = (c1.real * c2.real) - (c1.imag * c2.imag);
  c3.imag = (c1.real * c2.imag) + (c1.imag * c2.real);

  return c3;
}

Complex div(Complex c1, Complex c2) {
  Complex c3;

  c3.real = ((c1.real * c2.real) + (c1.imag * c2.imag)) /
            (pow(c2.real, 2) + pow(c2.imag, 2));
  c3.imag = -1 * (((c1.real * c2.imag) - (c1.imag * c2.real)) /
                  (pow(c2.real, 2) + pow(c2.imag, 2)));

  if (c3.imag == 0)
    c3.imag = 0;

  return c3;
}

double abs(Complex c) { return sqrt(c.real * c.real + c.imag * c.imag); }

void printComplex(Complex complex) {
  printf("%f + %fj ", complex.real, complex.imag);
}

//--------------------------------------------------------------------------------------------------
//------------------------------------ Matrix4cd y
//Vector4cd----------------------------------------
typedef Complex Matrix4cd[4][4];
typedef Complex Vector4cd[4];

void printMatrix4cd(Matrix4cd matrix) {
  int i, j;
  printf("MATRIZ 4x4: \n");
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      printf("(%f, %f) ", matrix[i][j].real, matrix[i][j].imag);
    }
    printf("\n");
  }
}

void printVector4cd(Vector4cd vector) {
  int i;
  printf("Vector: \n");
  for (i = 0; i < 4; i++) {
    printf("(%f, %f) ", vector[i].real, vector[i].imag);
  }
  printf("\n");
}

void initialise_vector4cd(Vector4cd vector) {
  int i;
  for (i = 0; i < 4; i++) {
    vector[i] = double_to_complex(0, 0);
  }
}

void mul4cd(Matrix4cd matrix, Vector4cd vector, Vector4cd result) {
  int i, j;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      result[i] = add(result[i], mul(matrix[i][j], vector[j]));
    }
  }
}

//----------------------------------------------------------------------------------------
//------------------------------MatrixXcd y
//VectorXcd-------------------------------------
typedef Complex __global *MatrixXcd;
typedef Complex __global *VectorXcd;

void printMatrixXcd(MatrixXcd matrix, size_t rows, size_t cols) {
  int i, j;
  printf("MATRIZ XxX: \n");
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      printf("(%f, %f) ", matrix[i * cols + j].real, matrix[i * cols + j].imag);
    }
    printf("\n");
  }
}

void printVectorXcd(VectorXcd vector, size_t length) {
  int i;
  printf("Vector X: \n");
  for (i = 0; i < length; i++) {
    printf("(%f, %f) ", vector[i].real, vector[i].imag);
  }
  printf("\n");
}

MatrixXcd initialize_matrixXcd(struct clheap __global *heap, size_t rows,
                               size_t cols) {
  MatrixXcd matrix = (MatrixXcd)malloc(heap, sizeof(Complex) * rows * cols);
  int i, j;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      matrix[i * cols + j] = double_to_complex(0, 0);
    }
  }

  return matrix;
}

VectorXcd initialize_vectorXcd(struct clheap __global *heap, size_t length) {
  VectorXcd vector = (VectorXcd)malloc(heap, sizeof(Complex) * length);
  int i;
  for (i = 0; i < length; i++) {
    vector[i] = double_to_complex(0, 0);
  }

  return vector;
}

VectorXcd mulXcd(struct clheap __global *heap, MatrixXcd matrix,
                 VectorXcd vector, size_t rows, size_t cols) {
  VectorXcd result = initialize_vectorXcd(heap, rows);

  int i, j;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      result[i] = add(result[i], mul(matrix[i * cols + j], vector[j]));
    }
  }

  return result;
}

typedef Complex CTYPE;
inline static double _cabs(CTYPE val) { return abs(val); }
inline static double _creal(CTYPE val) { return real(val); }
inline static double _cimag(CTYPE val) { return imag(val); }

//! dimension index
typedef unsigned long long ITYPE;

#define _USE_MATH_DEFINES

//! PI value
#ifndef PI
#ifdef M_PI
#define PI M_PI
#else
#define PI 3.141592653589793
#endif
#endif

//! square root of 2
#define SQRT2 1.414213562373095

//! cosine pi/8
#define COSPI8 0.923879532511287

//! sine pi/8
#define SINPI8 0.382683432365090

//! Maximum number of threads expected
#define MAX_NUM_THREADS 1024

//! Maximum number of threashold
#define PARALLEL_NQUBIT_THRESHOLD 64

/**
 * @file state.h
 * @brief Definition and basic functions for state vector
 */

/**
 * intiialize quantum state to zero state
 *
 * intiialize quantum state to zero state
 * @param[out] psi pointer of quantum state
 * @param[in] dim dimension
 */
void initialize_quantum_state(CTYPE __global *state, ITYPE dim);

void initialize_Haar_random_state(CTYPE __global *state, ITYPE dim);

void initialize_Haar_random_state_with_seed(CTYPE __global *state, ITYPE dim,
                                            UINT seed);

/**
 * Insert 0 to qubit_index-th bit of basis_index. basis_mask must be 1ULL <<
 * qubit_index.
 */
inline static ITYPE insert_zero_to_basis_index(ITYPE basis_index,
                                               ITYPE basis_mask,
                                               UINT qubit_index) {
  ITYPE temp_basis = (basis_index >> qubit_index) << (qubit_index + 1);
  return temp_basis + basis_index % basis_mask;
}

/**
 * Swap amplitude of state[basis_index_0] and state[basis_index_1]
 */
inline static void swap_amplitude(CTYPE __global *state, ITYPE basis_index_0,
                                  ITYPE basis_index_1) {
  CTYPE temp = state[basis_index_0];
  state[basis_index_0] = state[basis_index_1];
  state[basis_index_1] = temp;
}
/**
 * min for ITYPE
 */
inline static ITYPE get_min_ll(ITYPE index_0, ITYPE index_1) {
  return (index_0 < index_1) ? index_0 : index_1;
}
/**
 * max for ITYPE
 */
inline static ITYPE get_max_ll(ITYPE index_0, ITYPE index_1) {
  return (index_0 > index_1) ? index_0 : index_1;
}
/**
 * min for UINT
 */
inline static UINT get_min_ui(UINT index_0, UINT index_1) {
  return (index_0 < index_1) ? index_0 : index_1;
}
/**
 * max for UINT
 */
inline static UINT get_max_ui(UINT index_0, UINT index_1) {
  return (index_0 > index_1) ? index_0 : index_1;
}

/**
 * Count population (number of bit with 1) in 64bit unsigned integer.
 *
 * See http://developer.cybozu.co.jp/takesako/2006/11/binary_hacks.html for why
 * it works
 */
inline static UINT count_population(ITYPE x) {
  x = ((x & 0xaaaaaaaaaaaaaaaaUL) >> 1) + (x & 0x5555555555555555UL);
  x = ((x & 0xccccccccccccccccUL) >> 2) + (x & 0x3333333333333333UL);
  x = ((x & 0xf0f0f0f0f0f0f0f0UL) >> 4) + (x & 0x0f0f0f0f0f0f0f0fUL);
  x = ((x & 0xff00ff00ff00ff00UL) >> 8) + (x & 0x00ff00ff00ff00ffUL);
  x = ((x & 0xffff0000ffff0000UL) >> 16) + (x & 0x0000ffff0000ffffUL);
  x = ((x & 0xffffffff00000000UL) >> 32) + (x & 0x00000000ffffffffUL);
  return (UINT)x;
}

void sort_ui(UINT __global *v, int size, HEAP heap);
UINT __global *create_sorted_ui_list(UINT __global *array, size_t size,
                                     HEAP heap);
UINT __global *create_sorted_ui_list_value(UINT __global *array, size_t size,
                                           UINT value, HEAP heap);
UINT __global *create_sorted_ui_list_list(UINT __global *array1, size_t size1,
                                          UINT __global *array2, size_t size2,
                                          HEAP heap);

/**
 * Generate mask for all the combination of bit-shifting in qubit_index_list
 *
 * For example, when qubit_index_list = [0,3], the result is [0000, 0001, 1000,
 * 1001]. When two-qubit gate on 0- and 3-qubit, state[0010] is a linear
 * combination of amplitude about [0010, 0011, 1010, 1011], which can be
 * generated with bit-shifts with the above results.
 */
ITYPE __global *create_matrix_mask_list(UINT __global *qubit_index_list,
                                        UINT qubit_index_count, HEAP heap);
ITYPE
create_control_mask(UINT __global *qubit_index_list, UINT *value_list,
                    UINT qubit_index_count);

/**
 * Calculate bit-flip mask, phase-flip mask, global_phase, and pivot-qubit.
 *
 * Two masks are generated from Pauli operator.
 * Bit-flip mask is a bit-array of which the i-th bit has 1 if X or Y occurs on
 * the i-th qubit. Phase-flip mask is a bit-array of which the i-th bit has 1 if
 * Y or Z occurs on the i-th qubit.
 *
 * Suppose P is Pauli operator, |x> is computational basis, P|x> = alpha |y>,
 * and P|0> = alpha_0 |y_0>. We see y = x ^ bit_mask. alpha/alpha_0 is (-1)**|x
 * & phase_mask|, where || means the number of 1 in bit array, up to global
 * phase. alpha_0 is i**(the number of Pauli-Y in P).
 *
 * global_phase_90rot_count is the number of Pauli-Y in P.
 * pivot_qubit_index is the last qubit which is flipped by Pauli operator, which
 * is required in the next computation.
 */
void get_Pauli_masks_partial_list(UINT __global *target_qubit_index_list,
                                  UINT *Pauli_operator_type_list,
                                  UINT target_qubit_index_count,
                                  ITYPE *bit_flip_mask, ITYPE *phase_flip_mask,
                                  UINT *global_phase_90rot_count,
                                  UINT *pivot_qubit_index);
void get_Pauli_masks_whole_list(UINT *Pauli_operator_type_list,
                                UINT target_qubit_index_count,
                                ITYPE *bit_flip_mask, ITYPE *phase_flip_mask,
                                UINT *global_phase_90rot_count,
                                UINT *pivot_qubit_index);

/**
 * @file state.h
 * @brief Definition and basic functions for state vector
 */

/**
 * allocate quantum state in memory
 *
 * allocate quantum state in memory
 * @param[in] dim dimension, i.e. size of vector
 * @return pointer to allocated vector
 */
CTYPE __global *allocate_quantum_state(ITYPE dim, HEAP heap);

/**
 * release allocated quantum state
 *
 * release allocated quantum state
 * @param[in] psi quantum state
 */
void release_quantum_state(CTYPE __global *state, HEAP heap);

/**
 * @file state.h
 * @brief Definition and basic functions for state vector
 */

/**
 * allocate quantum state in memory
 *
 * allocate quantum state in memory
 * @param[in] dim dimension, i.e. size of vector
 * @return pointer to allocated vector
 */
CTYPE __global *dm_allocate_quantum_state(ITYPE dim, HEAP heap);

/**
 * intiialize quantum state to zero state
 *
 * intiialize quantum state to zero state
 * @param[out] state pointer of quantum state
 * @param[in] dim dimension
 */
void dm_initialize_quantum_state(CTYPE __global *state, ITYPE dim);

/**
 * release allocated quantum state
 *
 * release allocated quantum state
 * @param[in] state quantum state
 */
void dm_release_quantum_state(CTYPE __global *state, HEAP heap);

/**
 * initialize density matrix from pure state
 *
 * initialize density matrix from pure state
 * @param[out] state pointer of quantum state
 * @param[out] pure_state pointer of quantum state
 * @param[in] dim dimension
 */
void dm_initialize_with_pure_state(CTYPE __global *state,
                                   CTYPE __global *pure_state, ITYPE dim);

double state_norm_squared(CTYPE __global *state, ITYPE dim);
double state_norm_squared_single_thread(CTYPE __global *state, ITYPE dim);
double measurement_distribution_entropy(CTYPE __global *state, ITYPE dim);
CTYPE state_inner_product(CTYPE __global *state_bra, CTYPE __global *state_ket,
                          ITYPE dim);

void state_tensor_product(CTYPE __global *state_left, ITYPE dim_left,
                          CTYPE __global *state_right, ITYPE dim_right,
                          CTYPE *state_dst);
void state_permutate_qubit(UINT *qubit_order, CTYPE __global *state_src,
                           CTYPE *state_dst, UINT qubit_count, ITYPE dim);
void state_drop_qubits(UINT __global *target, UINT *projection,
                       UINT target_count, CTYPE __global *state_src,
                       CTYPE *state_dst, ITYPE dim, HEAP heap);

double M0_prob(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
double M1_prob(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
double marginal_prob( UINT* sorted_target_qubit_index_list,
     UINT __global* measured_value_list, UINT target_qubit_index_count,
     CTYPE __global* state, ITYPE dim) ;

double expectation_value_single_qubit_Pauli_operator(UINT target_qubit_index,
                                                     UINT Pauli_operator_type,
                                                     CTYPE __global *state,
                                                     ITYPE dim);
double expectation_value_multi_qubit_Pauli_operator_whole_list(
    UINT *Pauli_operator_type_list, UINT qubit_count, CTYPE __global *state,
    ITYPE dim);
double expectation_value_multi_qubit_Pauli_operator_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim);

CTYPE transition_amplitude_multi_qubit_Pauli_operator_whole_list(
    UINT *Pauli_operator_type_list, UINT qubit_count, CTYPE __global *state_bra,
    CTYPE __global *state_ket, ITYPE dim);
CTYPE transition_amplitude_multi_qubit_Pauli_operator_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state_bra,
    CTYPE __global *state_ket, ITYPE dim);

double expectation_value_multi_qubit_Pauli_operator_XZ_mask_single_thread(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, CTYPE __global *state, ITYPE dim);
double expectation_value_multi_qubit_Pauli_operator_Z_mask_single_thread(
    ITYPE phase_flip_mask, CTYPE __global *state, ITYPE dim);
double expectation_value_multi_qubit_Pauli_operator_partial_list_single_thread(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim);

double dm_state_norm_squared(CTYPE __global *state, ITYPE dim);
double dm_measurement_distribution_entropy(CTYPE __global *state, ITYPE dim);
void dm_state_add(CTYPE __global *state_added, CTYPE __global *state,
                  ITYPE dim);
void dm_state_add_with_coef(CTYPE coef, CTYPE __global *state_added,
                            CTYPE __global *state, ITYPE dim);
void dm_state_multiply(CTYPE coef, CTYPE __global *state, ITYPE dim);

void dm_state_tensor_product(CTYPE __global *state_left, ITYPE dim_left,
                             CTYPE __global *state_right, ITYPE dim_right,
                             CTYPE *state_dst);
void dm_state_permutate_qubit(UINT *qubit_order, CTYPE __global *state_src,
                              CTYPE *state_dst, UINT qubit_count, ITYPE dim);
void dm_state_partial_trace_from_density_matrix(UINT __global *target,
                                                UINT target_count,
                                                CTYPE __global *state_src,
                                                CTYPE *state_dst, ITYPE dim,
                                                HEAP heap);
void dm_state_partial_trace_from_state_vector(UINT __global *target,
                                              UINT target_count,
                                              CTYPE __global *state_src,
                                              CTYPE *state_dst, ITYPE dim,
                                              HEAP heap);

double dm_M0_prob(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
double dm_M1_prob(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
double dm_marginal_prob(UINT *sorted_target_qubit_index_list,
                        UINT *measured_value_list,
                        UINT target_qubit_index_count, CTYPE __global *state,
                        ITYPE dim);

double dm_expectation_value_multi_qubit_Pauli_operator_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim);

/*
 rules for arguments of update functions:
   The order of arguments must be
      information_about_applying_qubits -> Pauli_operator ->
 information_about_applying_gate_matrix_elements -> a_kind_of_rotation_angle ->
 state_vector -> dimension If there is control-qubit and target-qubit,
 control-qubit is the first argument. If an array of which the size is not known
 comes, the size of that array follows.

  Definition of update function is divided to named_gates,
 single_target_qubit_gates, multiple_target_qubit_gates, QFT_gates
 */

/** X gate **/

/**
 * \~english
 * Apply the Pauli X gate to the quantum state.
 *
 * Apply the Pauli X gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * パウリX演算を作用させて状態を更新。
 *
 * パウリX演算を作用させて状態を更新。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void X_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void X_gate_parallel_unroll(UINT target_qubit_index, CTYPE __global *state,
                            ITYPE dim);
void X_gate_parallel_simd(UINT target_qubit_index, CTYPE __global *state,
                          ITYPE dim);
void X_gate_parallel_sve(UINT target_qubit_index, CTYPE __global *state,
                         ITYPE dim);

/**
 * \~english
 * Apply the Pauli Y gate to the quantum state.
 *
 * Apply the Pauli Y gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * パウリY演算を作用させて状態を更新。
 *
 * パウリY演算を作用させて状態を更新。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void Y_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void Y_gate_parallel_unroll(UINT target_qubit_index, CTYPE __global *state,
                            ITYPE dim);
void Y_gate_parallel_simd(UINT target_qubit_index, CTYPE __global *state,
                          ITYPE dim);
void Y_gate_parallel_sve(UINT target_qubit_index, CTYPE __global *state,
                         ITYPE dim);

/**
 * \~english
 * Apply the Pauli Z gate to the quantum state.
 *
 * Apply the Pauli Z gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * パウリZ演算を作用させて状態を更新。
 *
 * パウリZ演算を作用させて状態を更新。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void Z_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void Z_gate_parallel_unroll(UINT target_qubit_index, CTYPE __global *state,
                            ITYPE dim);
void Z_gate_parallel_simd(UINT target_qubit_index, CTYPE __global *state,
                          ITYPE dim);
void Z_gate_parallel_sve(UINT target_qubit_index, CTYPE __global *state,
                         ITYPE dim);

/**
 * \~english
 * Apply S gate to the quantum state.
 *
 * Apply S gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 位相演算 S を作用させて状態を更新。
 *
 * 位相演算 S = diag(1,i) を作用させて状態を更新。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void S_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply S gate to the quantum state.
 *
 * Apply S gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 位相演算 S^dag を作用させて状態を更新。
 *
 * 位相演算 S^dag = diag(1,-i) を作用させて状態を更新。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void Sdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply T gate to the quantum state.
 *
 * Apply T gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * T 演算を作用させて状態を更新。
 *
 * T 演算（pi/8演算）、 T = diag(1,exp(i pi/4))
 * を作用させて状態を更新。非クリフォード演算。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void T_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply T^dag gate to the quantum state.
 *
 * Apply T^dag gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * T^dag 演算を作用させて状態を更新。
 *
 * T 演算のエルミート共役、 T^dag = diag(1,exp(-i pi/4))
 * を作用させて状態を更新。非クリフォード演算。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void Tdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply the square root of the X gate to the quantum state.
 *
 * Apply the square root of the X gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * パウリ X 演算子の平方根の演算子を作用させて状態を更新。
 *
 * パウリ X 演算子の平方根の演算子を作用させて状態を更新。非クリフォード演算。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void sqrtX_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply hermitian conjugate of the square root of the X gate to the quantum
 * state.
 *
 * Apply hermitian conjugate of the square root of the X gate to the quantum
 * state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * パウリ X 演算子の平方根に対してエルミート共役な演算子を作用させて状態を更新。
 *
 * パウリ X
 * 演算子の平方根に対してエルミート共役な演算子を作用させて状態を更新。非クリフォード演算。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void sqrtXdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply the square root of the Y gate to the quantum state.
 *
 * Apply the square root of the Y gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * パウリ Y 演算子の平方根の演算子を作用させて状態を更新。
 *
 * パウリ Y 演算子の平方根の演算子を作用させて状態を更新。非クリフォード演算。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void sqrtY_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply hermitian conjugate of the square root of the Y gate to the quantum
 * state.
 *
 * Apply hermitian conjugate of the square root of the Y gate to the quantum
 * state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * パウリ Y 演算子の平方根に対してエルミート共役な演算子を作用させて状態を更新。
 *
 * パウリ Y
 * 演算子の平方根に対してエルミート共役な演算子を作用させて状態を更新。非クリフォード演算。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void sqrtYdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply the Hadamard gate to the quantum state.
 *
 * Apply the Hadamard gate to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * アダマール演算子を作用させて状態を更新。
 *
 * アダマール演算子を作用させて状態を更新。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void H_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void H_gate_parallel_unroll(UINT target_qubit_index, CTYPE __global *state,
                            ITYPE dim);
void H_gate_parallel_simd(UINT target_qubit_index, CTYPE __global *state,
                          ITYPE dim);
void H_gate_parallel_sve(UINT target_qubit_index, CTYPE __global *state,
                         ITYPE dim);

/** Hadamard gate multiplied sqrt(2) **/
// void H_gate_unnormalized(UINT target_qubit_index, CTYPE __global* state,
// ITYPE dim);

/**
 * \~english
 * Apply the CNOT gate to the quantum state.
 *
 * Apply the CNOT gate to the quantum state.
 * @param[in] control_qubit_index index of control qubit
 * @param[in] target_qubit_index index of target qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 * \~japanese-en
 * CNOT演算を作用させて状態を更新。
 *
 * 2量子ビット演算、CNOT演算を作用させて状態を更新。
 * @param[in] control_qubit_index 制御量子ビットのインデックス
 * @param[in] target_qubit_index ターゲット量子ビットのインデックス
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 */
void CNOT_gate(UINT control_qubit_index, UINT target_qubit_index,
               CTYPE __global *state, ITYPE dim);
void CNOT_gate_parallel_unroll(UINT control_qubit_index,
                               UINT target_qubit_index, CTYPE __global *state,
                               ITYPE dim);
void CNOT_gate_parallel_simd(UINT control_qubit_index, UINT target_qubit_index,
                             CTYPE __global *state, ITYPE dim);
void CNOT_gate_parallel_sve(UINT control_qubit_index, UINT target_qubit_index,
                            CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply the CZ gate to the quantum state.
 *
 * Apply the CZ gate to the quantum state.
 * @param[in] control_qubit_index index of control qubit
 * @param[in] target_qubit_index index of target qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 * \~japanese-en
 * CZ演算を作用させて状態を更新。
 *
 * 2量子ビット演算、CZ演算を作用させて状態を更新。制御量子ビットとターゲット量子ビットに対して対称に作用（インデックスを入れ替えても同じ作用）。
 * @param[in] control_qubit_index 制御量子ビットのインデックス
 * @param[in] target_qubit_index ターゲット量子ビットのインデックス
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 */
void CZ_gate(UINT control_qubit_index, UINT target_qubit_index,
             CTYPE __global *state, ITYPE dim);
void CZ_gate_parallel_unroll(UINT control_qubit_index, UINT target_qubit_index,
                             CTYPE __global *state, ITYPE dim);
void CZ_gate_parallel_simd(UINT control_qubit_index, UINT target_qubit_index,
                           CTYPE __global *state, ITYPE dim);
void CZ_gate_parallel_sve(UINT control_qubit_index, UINT target_qubit_index,
                          CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply the SWAP to the quantum state.
 *
 * Apply the SWAP to the quantum state.
 * @param[in] control_qubit_index index of control qubit
 * @param[in] target_qubit_index index of target qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 * \~japanese-en
 * SWAP演算を作用させて状態を更新。
 *
 * 2量子ビット演算、SWAP演算を作用させて状態を更新。２つの量子ビットに対して対称に作用する（インデックスを入れ替えても同じ作用）。
 * @param[in] target_qubit_index_0 作用する量子ビットのインデックス
 * @param[in] target_qubit_index_1 作用する量子ビットのインデックス
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 */
void SWAP_gate(UINT target_qubit_index_0, UINT target_qubit_index_1,
               CTYPE __global *state, ITYPE dim);
void SWAP_gate_parallel_unroll(UINT target_qubit_index_0,
                               UINT target_qubit_index_1, CTYPE __global *state,
                               ITYPE dim);
void SWAP_gate_parallel_simd(UINT target_qubit_index_0,
                             UINT target_qubit_index_1, CTYPE __global *state,
                             ITYPE dim);
void SWAP_gate_parallel_sve(UINT target_qubit_index_0,
                            UINT target_qubit_index_1, CTYPE __global *state,
                            ITYPE dim);

/**
 * \~english
 * Project the quantum state to the 0 state.
 *
 * Project the quantum state to the 0 state. The output state is not normalized.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 * \~japanese-en
 * 0状態への射影
 *
 * 0状態への射影演算子 |0><0|
 * を作用させて状態を更新する。ノルムは規格化されない。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void P0_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void P0_gate_single(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void P0_gate_parallel(UINT target_qubit_index, CTYPE __global *state,
                      ITYPE dim);

/**
 * \~english
 * Project the quantum state to the 1 state.
 *
 * Project the quantum state to the 1 state. The output state is not normalized.
 * @param[in] target_qubit_index index of the qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 * \~japanese-en
 * 1状態への射影
 *
 * 1状態への射影演算子 |1><1|
 * を作用させて状態を更新する。ノルムは規格化されない。
 * @param[in] target_qubit_index 作用する量子ビットのインデックス
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void P1_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void P1_gate_single(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void P1_gate_parallel(UINT target_qubit_index, CTYPE __global *state,
                      ITYPE dim);

/**
 * \~english
 * Normalize the quantum state.
 *
 * Normalize the quantum state by multiplying the normalization factor
 * 1/sqrt(norm).
 * @param[in] norm norm of the state
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 * \~japanese-en
 * 状態を規格化
 *
 * 状態に対して 1/sqrt(norm) 倍をする。norm
 * が量子状態のノルムである場合にはノルムが1になるように規格化される。
 * @param[in] norm ノルム
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 */
void normalize(double squared_norm, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Normalize the quantum state.
 *
 * Normalize the quantum state by multiplying the normalization factor
 * 1/sqrt(norm).
 * @param[in] norm norm of the state
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 * \~japanese-en
 * 状態を規格化
 *
 * 状態に対して 1/sqrt(norm) 倍をする。norm
 * が量子状態のノルムである場合にはノルムが1になるように規格化される。
 * @param[in] norm ノルム
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 */
void normalize_single_thread(double squared_norm, CTYPE __global *state,
                             ITYPE dim);

/**
 * \~english
 * Apply a X rotation gate by angle to the quantum state.
 *
 * Apply a X rotation gate by angle to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in] angle angle of the rotation
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * X軸の回転演算を作用させて状態を更新
 *
 * X軸の回転演算
 *
 * cos (angle/2) I + i sin (angle/2) X
 *
 * を作用させて状態を更新。angleは回転角。
 * @param[in] target_qubit_index 量子ビットのインデックス
 * @param[in] angle 回転角
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void RX_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
             ITYPE dim);

/**
 * \~english
 * Apply a Y rotation gate by angle to the quantum state.
 *
 * Apply a Y rotation gate by angle to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in] angle angle of the rotation
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * Y軸の回転演算を作用させて状態を更新
 *
 * Y軸の回転演算
 *
 * cos (angle/2) I + i sin (angle/2) Y
 *
 * を作用させて状態を更新。angleは回転角。
 * @param[in] target_qubit_index 量子ビットのインデックス
 * @param[in] angle 回転角
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void RY_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
             ITYPE dim);

/**
 * \~english
 * Apply a Z rotation gate by angle to the quantum state.
 *
 * Apply a Z rotation gate by angle to the quantum state.
 * @param[in] target_qubit_index index of the qubit
 * @param[in] angle angle of the rotation
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * Z軸の回転演算を作用させて状態を更新
 *
 * Z軸の回転演算
 *
 * cos (angle/2) I + i sin (angle/2) Z
 *
 * を作用させて状態を更新。angleは回転角。
 * @param[in] target_qubit_index 量子ビットのインデックス
 * @param[in] angle 回転角
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void RZ_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
             ITYPE dim);

/**
 * \~english
 * Apply a single-qubit Pauli operator to the quantum state.
 *
 * Apply a single-qubit Pauli operator to the quantum state. Pauli_operator_type
 * must be 0,1,2,3 corresponding to the Pauli I, X, Y, Z operators respectively.
 *
 *
 * @param[in] target_qubit_index index of the qubit
 * @param[in] Pauli_operator_type type of the Pauli operator 0,1,2,3
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * パウリ演算子を作用させて状態を更新
 *
 * パウリ演算子を作用させて状態を更新。Pauli_operator_type はパウリ演算子 I, X,
 * Y, Z に対応して 0,1,2,3 を指定。
 * @param[in] target_qubit_index 量子ビットのインデックス
 * @param[in] Pauli_operator_type
 * 作用するパウリ演算子のタイプ、0,1,2,3。それぞれI, X, Y, Zに対応。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 *
 */
void single_qubit_Pauli_gate(UINT target_qubit_index, UINT Pauli_operator_type,
                             CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply a single-qubit Pauli rotation operator to the quantum state.
 *
 * Apply a single-qubit Pauli rotation operator to the quantum state.
 * Pauli_operator_type must be 0,1,2,3 corresponding to the Pauli I, X, Y, Z
 * operators respectively.
 *
 *
 * @param[in] target_qubit_index index of the qubit
 * @param[in] Pauli_operator_type type of the Pauli operator 0,1,2,3
 * @param[in] angle rotation angle
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * パウリ演算子の回転演算を作用させて状態を更新
 *
 * パウリ演算子 A = I, X, Y, Zに対する回転角 angle の回転演算
 *
 * cos (angle/2) I + i sin (angle/2) A
 *
 * を作用させて状態を更新。Pauli_operator_type はパウリ演算子 I, X, Y, Z
 * に対応して 0,1,2,3 を指定。
 * @param[in] target_qubit_index 量子ビットのインデックス
 * @param[in] Pauli_operator_type
 * 作用するパウリ演算子のタイプ、0,1,2,3。それぞれI, X, Y, Zに対応。
 * @param[in] angle 回転角
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 *
 */
void single_qubit_Pauli_rotation_gate(UINT target_qubit_index,
                                      UINT Pauli_operator_index, double angle,
                                      CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply a single-qubit dense operator to the quantum state.
 *
 * Apply a single-qubit dense operator to the quantum state.
 *
 * @param[in] target_qubit_index index of the qubit
 * @param[in] matrix[4] description of the dense matrix as one-dimensional array
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 任意の１量子ビット演算を作用させて状態を更新
 *
 * 任意の１量子ビット演算を作用させて状態を更新。１量子ビット演算は、4つの要素をもつ１次元配列
 * matrix[4] によって指定される。
 *
 * 例）パウリX演算子：{0, 1, 1,
 * 0}、アダマール演算子：{1/sqrt(2),1/sqrt(2),1/sqrt(2),-1/sqrt(2)}。
 *
 * @param[in] target_qubit_index 量子ビットのインデックス
 * @param[in] matrix[4] １量子ビット演算を指定する4次元配列
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void single_qubit_dense_matrix_gate(UINT target_qubit_index,
                                    __constant CTYPE matrix[4],
                                    CTYPE __global *state, ITYPE dim);
void single_qubit_dense_matrix_gate_non_constant(UINT target_qubit_index,
                                                 CTYPE matrix[4],
                                                 CTYPE __global *state,
                                                 ITYPE dim);
void single_qubit_dense_matrix_gate_parallel_unroll(UINT target_qubit_index,
                                                    CTYPE matrix[4],
                                                    CTYPE __global *state,
                                                    ITYPE dim);
void single_qubit_dense_matrix_gate_parallel(UINT target_qubit_index,
                                             __constant CTYPE matrix[4],
                                             CTYPE __global *state, ITYPE dim);
void single_qubit_dense_matrix_gate_parallel_no_constant(
    UINT target_qubit_index, CTYPE matrix[4], CTYPE __global *state, ITYPE dim);
void single_qubit_dense_matrix_gate_parallel_simd(UINT target_qubit_index,
                                                  CTYPE matrix[4],
                                                  CTYPE __global *state,
                                                  ITYPE dim);
void single_qubit_dense_matrix_gate_parallel_sve(UINT target_qubit_index,
                                                 CTYPE matrix[4],
                                                 CTYPE __global *state,
                                                 ITYPE dim);

/**
 * \~english
 * Apply a single-qubit diagonal operator to the quantum state.
 *
 * Apply a single-qubit diagonal operator to the quantum state.
 *
 * @param[in] target_qubit_index index of the qubit
 * @param[in] diagonal_matrix[2] description of the single-qubit diagonal
 * elements
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * １量子ビットの対角演算を作用させて状態を更新
 *
 * １量子ビットの対角演算を作用させて状態を更新。１量子ビットの対角演算は、その対角成分を定義する2つの要素から成る１次元配列
 * diagonal_matrix[2] によって指定される。
 *
 * @param[in] target_qubit_index 量子ビットのインデックス
 * @param[in] diagonal_matrix[2] ２つの対角成分を定義する１次元配列
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void single_qubit_diagonal_matrix_gate(UINT target_qubit_index,
                                       CTYPE diagonal_matrix[2],
                                       CTYPE __global *state, ITYPE dim);
void single_qubit_diagonal_matrix_gate_parallel_unroll(UINT target_qubit_index,
                                                       CTYPE diagonal_matrix[2],
                                                       CTYPE __global *state,
                                                       ITYPE dim);
void single_qubit_diagonal_matrix_gate_parallel_simd(UINT target_qubit_index,
                                                     CTYPE diagonal_matrix[2],
                                                     CTYPE __global *state,
                                                     ITYPE dim);
void single_qubit_diagonal_matrix_gate_parallel_sve(UINT target_qubit_index,
                                                    CTYPE diagonal_matrix[2],
                                                    CTYPE __global *state,
                                                    ITYPE dim);

/**
 * \~english
 * Apply a single-qubit phase operator to the quantum state.
 *
 * Apply a single-qubit phase operator, diag(1,phsae), to the quantum state.
 *
 * @param[in] target_qubit_index index of the qubit
 * @param[in] phase phase factor
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * １量子ビットの一般位相演算を作用させて状態を更新
 *
 * １量子ビットの一般位相演算 diag(1,phase) を作用させて状態を更新。|1>状態が
 * phase 倍される。
 *
 * @param[in] target_qubit_index 量子ビットのインデックス
 * @param[in] phase 位相因子
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void single_qubit_phase_gate(UINT target_qubit_index, CTYPE phase,
                             CTYPE __global *state, ITYPE dim);
void single_qubit_phase_gate_parallel_unroll(UINT target_qubit_index,
                                             CTYPE phase, CTYPE __global *state,
                                             ITYPE dim);
void single_qubit_phase_gate_parallel_simd(UINT target_qubit_index, CTYPE phase,
                                           CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply a single-qubit controlled single-qubit gate.
 *
 * Apply a single-qubit controlled single-qubit gate.
 *
 * @param[in] control_qubit_index index of the control qubit
 * @param[in] control_value value of the control qubit
 * @param[in] target_qubit_index index of the target qubit
 * @param[in] matrix[4] description of the single-qubit dense matrix
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 単一量子ビットを制御量子ビットとする単一量子ビット演算
 *
 * 単一量子ビットを制御量子ビットとする単一量子ビット演算。制御量子ビット |0>
 * もしくは |1>のどちらの場合に作用するかは control_value = 0,1
 * によって指定。作用する単一量子ビット演算は matrix[4]
 * によって１次元配列として定義。
 *
 * @param[in] control_qubit_index 制御量子ビットのインデックス
 * @param[in] control_value 制御量子ビットの値
 * @param[in] target_qubit_index ターゲット量子ビットのインデックス
 * @param[in] matrix[4] 単一量子ビットの記述を与える１次元配列
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void single_qubit_control_single_qubit_dense_matrix_gate(
    UINT control_qubit_index, UINT control_value, UINT target_qubit_index,
    CTYPE matrix[4], CTYPE __global *state, ITYPE dim);
void single_qubit_control_single_qubit_dense_matrix_gate_unroll(
    UINT control_qubit_index, UINT control_value, UINT target_qubit_index,
    CTYPE matrix[4], CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply a multi-qubit controlled single-qubit gate.
 *
 * Apply a multi-qubit controlled single-qubit gate.
 *
 * @param[in] control_qubit_index_list list of the indexes of the control qubits
 * @param[in] control_value_list list of the vlues of the control qubits
 * @param[in] control_qubit_index_count the number of the control qubits
 * @param[in] target_qubit_index index of the target qubit
 * @param[in] matrix[4] description of the single-qubit dense matrix
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 複数量子ビットを制御量子ビットとする単一量子ビット演算
 *
 * 複数量子ビットを制御量子ビットとする単一量子ビット演算。control_qubit_index_count
 * 個の制御量子ビットについて、どの状態の場合に制御演算が作用するかは
 * control_value_list によって指定。作用する単一量子ビット演算は matrix[4]
 * によって１次元配列として定義。
 *
 * @param[in] control_qubit_index_list 制御量子ビットのインデックスのリスト
 * @param[in] control_value_list 制御演算が作用する制御量子ビットの値のリスト
 * @param[in] control_qubit_index_count 制御量子ビットの数
 * @param[in] target_qubit_index ターゲット量子ビットのインデックス
 * @param[in] matrix[4] 単一量子ビットの記述を与える１次元配列
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_control_single_qubit_dense_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT target_qubit_index, CTYPE matrix[4],
    CTYPE __global *state, ITYPE dim, HEAP heap);
void multi_qubit_control_single_qubit_dense_matrix_gate_unroll(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT target_qubit_index, CTYPE matrix[4],
    CTYPE __global *state, ITYPE dim, HEAP heap);

/**
 * \~english
 * Apply multi-qubit Pauli operator to the quantum state with a whole list.
 *
 * Apply multi-qubit Pauli operator to the quantum state with a whole list of
 * the Pauli operators. Pauli_operator_type_list must be a list of n single
 * Pauli operators.
 *
 * @param[in] Pauli_operator_type_list list of {0,1,2,3} corresponding to I, X,
 * Y, Z for all qubits
 * @param[in] qubit_count the number of the qubits
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * すべての量子ビットに多体のパウリ演算子を作用させて状態を更新
 *
 * 全ての量子ビットに対するパウリ演算子を与えて、多体のパウリ演算子を作用させて状態を更新。Pauli_operator_type_list
 * には全ての量子ビットに対するパウリ演算子を指定。
 *
 * 例）５量子ビット系の３つめの量子ビットへのX演算、４つめの量子ビットへのZ演算、IZXII：{0,0,1,3,0}（量子ビットの順番は右から数えている）
 *
 * @param[in] Pauli_operator_type_list 長さ qubit_count の {0,1,2,3} のリスト。
 * 0,1,2,3 はそれぞれパウリ演算子 I, X, Y, Z に対応。
 * @param[in] qubit_count 量子ビットの数
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_Pauli_gate_whole_list(UINT *Pauli_operator_type_list,
                                       UINT qubit_count, CTYPE __global *state,
                                       ITYPE dim_);

/**
 * \~english
 * Apply multi-qubit Pauli operator to the quantum state with a whole list.
 *
 * Apply multi-qubit Pauli operator to the quantum state with a whole list of
 * the Pauli operators. Pauli_operator_type_list must be a list of n single
 * Pauli operators.
 *
 * @param[in] Pauli_operator_type_list list of {0,1,2,3} corresponding to I, X,
 * Y, Z for all qubits
 * @param[in] qubit_count the number of the qubits
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * すべての量子ビットに多体のパウリ演算子を作用させて状態を更新
 *
 * 全ての量子ビットに対するパウリ演算子を与えて、多体のパウリ演算子を作用させて状態を更新。Pauli_operator_type_list
 * には全ての量子ビットに対するパウリ演算子を指定。
 *
 * 例）５量子ビット系の３つめの量子ビットへのX演算、４つめの量子ビットへのZ演算、IZXII：{0,0,1,3,0}（量子ビットの順番は右から数えている）
 *
 * @param[in] Pauli_operator_type_list 長さ qubit_count の {0,1,2,3} のリスト。
 * 0,1,2,3 はそれぞれパウリ演算子 I, X, Y, Z に対応。
 * @param[in] qubit_count 量子ビットの数
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_Pauli_gate_whole_list_single_thread(
    UINT *Pauli_operator_type_list, UINT qubit_count, CTYPE __global *state,
    ITYPE dim_);

/**
 * \~english
 * Apply multi-qubit Pauli operator to the quantum state with a partial list.
 *
 * Apply multi-qubit Pauli operator to the quantum state with a partial list of
 * the Pauli operators. Pauli_operator_type_list must be a list of n single
 * Pauli operators.
 *
 * @param[in] target_qubit_index_list list of the target qubits
 * @param[in] Pauli_operator_type_list list of {0,1,2,3} corresponding to I, X,
 * Y, Z for the target qubits
 * @param[in] target_qubit_index_count the number of the target qubits
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * すべての量子ビットに多体のパウリ演算子を作用させて状態を更新
 *
 * 全ての量子ビットに対するパウリ演算子を与えて、多体のパウリ演算子を作用させて状態を更新。Pauli_operator_type_list
 * には全ての量子ビットに対するパウリ演算子を指定。
 *
 * 例）５量子ビット系の３つめの量子ビットへのX演算、４つめの量子ビットへのZ演算、IZXII：Pauli_operator_type_list
 * ={2,3}, Pauli_operator_type_list
 * ={1,3}（量子ビットの順番は右から数えている）.
 *
 * @param[in] target_qubit_index_list 作用する量子ビットのインデックスのリスト。
 * @param[in] Pauli_operator_type_list
 * 作用する量子ビットのみに対するパウリ演算子を指定する、長さ
 * target_qubit_index_count の{0,1,2,3}のリスト。0,1,2,3 はそれぞれパウリ演算子
 * I, X, Y, Z に対応。
 * @param[in] target_qubit_index_count 作用する量子ビットの数
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_Pauli_gate_partial_list(UINT __global *target_qubit_index_list,
                                         UINT *Pauli_operator_type_list,
                                         UINT target_qubit_index_count,
                                         CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply multi-qubit Pauli operator to the quantum state with a partial list.
 *
 * Apply multi-qubit Pauli operator to the quantum state with a partial list of
 * the Pauli operators. Pauli_operator_type_list must be a list of n single
 * Pauli operators.
 *
 * @param[in] target_qubit_index_list list of the target qubits
 * @param[in] Pauli_operator_type_list list of {0,1,2,3} corresponding to I, X,
 * Y, Z for the target qubits
 * @param[in] target_qubit_index_count the number of the target qubits
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * すべての量子ビットに多体のパウリ演算子を作用させて状態を更新
 *
 * 全ての量子ビットに対するパウリ演算子を与えて、多体のパウリ演算子を作用させて状態を更新。Pauli_operator_type_list
 * には全ての量子ビットに対するパウリ演算子を指定。
 *
 * 例）５量子ビット系の３つめの量子ビットへのX演算、４つめの量子ビットへのZ演算、IZXII：Pauli_operator_type_list
 * ={2,3}, Pauli_operator_type_list
 * ={1,3}（量子ビットの順番は右から数えている）.
 *
 * @param[in] target_qubit_index_list 作用する量子ビットのインデックスのリスト。
 * @param[in] Pauli_operator_type_list
 * 作用する量子ビットのみに対するパウリ演算子を指定する、長さ
 * target_qubit_index_count の{0,1,2,3}のリスト。0,1,2,3 はそれぞれパウリ演算子
 * I, X, Y, Z に対応。
 * @param[in] target_qubit_index_count 作用する量子ビットの数
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_Pauli_gate_partial_list_single_thread(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply multi-qubit Pauli rotation operator to the quantum state with a whole
 * list.
 *
 * Apply multi-qubit Pauli rotation operator to state with a whole list of the
 * Pauli operators. Pauli_operator_type_list must be a list of n single Pauli
 * operators. Update a quantum state with a mutliple qubits Pauli rotation, cos
 * (angle/2 ) I + i sin ( angle/2 ) A, where A is the Pauli operator specified
 * by Pauli_operator.
 *
 * @param[in] Pauli_operator_type_list array of {0,1,2,3} of the length
 * n_qubits. 0,1,2,3 corresponds to i,x,y,z
 * @param[in] qubit_count number of the qubits
 * @param[in] angle rotation angle
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 複数量子ビット（全て指定）のパウリ演算子による回転演算を用いて状態を更新。
 *
 * 複数量子ビット（全て指定）のパウリ演算子 A による回転演算
 *
 * cos ( angle/2 ) I + i sin ( angle/2 ) A
 *
 * を用いて状態を更新。Pauli_operator_type_list
 * は全ての量子ビットに対するパウリ演算子のリスト。パウリ演算子はI, X, Y,
 * Zがそれぞれ 0,1,2,3 に対応。
 *
 * 例) ５量子ビットに対する YIZXI
 * の場合は、{0,1,3,0,2}（量子ビットの順番は右から数えている）。
 *
 * @param[in] Pauli_operator_type_list 長さ qubit_count の {0,1,2,3} のリスト。
 * 0,1,2,3 はそれぞれパウリ演算子 I, X, Y, Z に対応。
 * @param[in] qubit_count 量子ビットの数
 * @param[in] angle 回転角度
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_Pauli_rotation_gate_whole_list(UINT *Pauli_operator_type_list,
                                                UINT qubit_count, double angle,
                                                CTYPE __global *state,
                                                ITYPE dim_);

/**
 * \~english
 * Apply multi-qubit Pauli rotation operator to the quantum state with a whole
 * list.
 *
 * Apply multi-qubit Pauli rotation operator to state with a whole list of the
 * Pauli operators. Pauli_operator_type_list must be a list of n single Pauli
 * operators. Update a quantum state with a mutliple qubits Pauli rotation, cos
 * (angle/2 ) I + i sin ( angle/2 ) A, where A is the Pauli operator specified
 * by Pauli_operator.
 *
 * @param[in] Pauli_operator_type_list array of {0,1,2,3} of the length
 * n_qubits. 0,1,2,3 corresponds to i,x,y,z
 * @param[in] qubit_count number of the qubits
 * @param[in] angle rotation angle
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 複数量子ビット（全て指定）のパウリ演算子による回転演算を用いて状態を更新。
 *
 * 複数量子ビット（全て指定）のパウリ演算子 A による回転演算
 *
 * cos ( angle/2 ) I + i sin ( angle/2 ) A
 *
 * を用いて状態を更新。Pauli_operator_type_list
 * は全ての量子ビットに対するパウリ演算子のリスト。パウリ演算子はI, X, Y,
 * Zがそれぞれ 0,1,2,3 に対応。
 *
 * 例) ５量子ビットに対する YIZXI
 * の場合は、{0,1,3,0,2}（量子ビットの順番は右から数えている）。
 *
 * @param[in] Pauli_operator_type_list 長さ qubit_count の {0,1,2,3} のリスト。
 * 0,1,2,3 はそれぞれパウリ演算子 I, X, Y, Z に対応。
 * @param[in] qubit_count 量子ビットの数
 * @param[in] angle 回転角度
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_Pauli_rotation_gate_whole_list_single_thread(
    UINT *Pauli_operator_type_list, UINT qubit_count, double angle,
    CTYPE __global *state, ITYPE dim_);

/**
 * \~english
 * Apply multi-qubit Pauli rotation operator to the quantum state with a partial
 * list.
 *
 * Apply multi-qubit Pauli rotation operator to state with a partial list of the
 * Pauli operators.
 *
 * @param[in] target_qubit_index_list list of the target qubits
 * @param[in] Pauli_operator_type_list list of {0,1,2,3} of length
 * target_qubit_index_count. {0,1,2,3} corresponds to I, X, Y, Z for the target
 * qubits.
 * @param[in] target_qubit_index_count the number of the target qubits
 * @param[in] angle rotation angle
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 複数量子ビット（部分的に指定）のパウリ演算子による回転演算を用いて状態を更新。
 *
 * 複数量子ビット（部分的に指定）のパウリ演算子 A による回転演算
 *
 * cos ( angle/2 ) I + i sin ( angle/2 ) A
 *
 * を用いて状態を更新。パウリ演算子Aは、target_qubit_index_list
 * で指定される一部の量子ビットに対するパウリ演算子のリスト
 * Pauli_operator_type_list
 * によって定義。回転角には(1/2)の因子はついていない。パウリ演算子はI, X, Y,
 * Zがそれぞれ 0,1,2,3 に対応。
 *
 * 例) ５量子ビットに対する YIZXI の場合は、target_qubit_index_list = {1,2,4},
 * Pauli_operator_type_list = {1,3,2}（量子ビットの順番は右から数えている）。
 *
 * @param[in] target_qubit_index_list 作用する量子ビットのインデックス
 * @param[in] Pauli_operator_type_list 長さ target_qubit_index_countの {0,1,2,3}
 * のリスト。0,1,2,3は、それぞれパウリ演算子 I, X, Y, Z 対応。
 * @param[in] target_qubit_index_count 作用する量子ビットの数
 * @param[in] angle 回転角
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_Pauli_rotation_gate_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, double angle, CTYPE __global *state,
    ITYPE dim);

/**
 * \~english
 * Apply multi-qubit Pauli rotation operator to the quantum state with a partial
 * list.
 *
 * Apply multi-qubit Pauli rotation operator to state with a partial list of the
 * Pauli operators.
 *
 * @param[in] target_qubit_index_list list of the target qubits
 * @param[in] Pauli_operator_type_list list of {0,1,2,3} of length
 * target_qubit_index_count. {0,1,2,3} corresponds to I, X, Y, Z for the target
 * qubits.
 * @param[in] target_qubit_index_count the number of the target qubits
 * @param[in] angle rotation angle
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 複数量子ビット（部分的に指定）のパウリ演算子による回転演算を用いて状態を更新。
 *
 * 複数量子ビット（部分的に指定）のパウリ演算子 A による回転演算
 *
 * cos ( angle/2 ) I + i sin ( angle/2 ) A
 *
 * を用いて状態を更新。パウリ演算子Aは、target_qubit_index_list
 * で指定される一部の量子ビットに対するパウリ演算子のリスト
 * Pauli_operator_type_list
 * によって定義。回転角には(1/2)の因子はついていない。パウリ演算子はI, X, Y,
 * Zがそれぞれ 0,1,2,3 に対応。
 *
 * 例) ５量子ビットに対する YIZXI の場合は、target_qubit_index_list = {1,2,4},
 * Pauli_operator_type_list = {1,3,2}（量子ビットの順番は右から数えている）。
 *
 * @param[in] target_qubit_index_list 作用する量子ビットのインデックス
 * @param[in] Pauli_operator_type_list 長さ target_qubit_index_countの {0,1,2,3}
 * のリスト。0,1,2,3は、それぞれパウリ演算子 I, X, Y, Z 対応。
 * @param[in] target_qubit_index_count 作用する量子ビットの数
 * @param[in] angle 回転角
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_Pauli_rotation_gate_partial_list_single_thread(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, double angle, CTYPE __global *state,
    ITYPE dim);

/**
 * \~english
 * Apply a two-qubit arbitrary gate.
 *
 * Apply a two-qubit arbitrary gate defined by a dense matrix as a
 * one-dimentional array matrix[].
 * @param[in] target_qubit_index1 index of the first target qubit
 * @param[in] target_qubit_index2 index of the second target qubit
 * @param[in] matrix description of a multi-qubit gate as a one-dimensional
 * array with length 16
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 任意の２量子ビット演算を作用させて状態を更新。
 *
 * 任意の２量子ビット演算を作用させて状態を更新。演算は、その行列成分を１次元配列として
 * matrix[] で与える。(j,k)成分は、matrix[dim*j+k]に対応。
 *
 * @param[in] target_qubit_index1 作用する量子ビットの一つ目の添え字
 * @param[in] target_qubit_index2 作用する量子ビットの二つ目の添え字
 * @param[in] matrix 複数量子ビット演算を定義する長さ 16 の一次元配列。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void double_qubit_dense_matrix_gate_c(UINT target_qubit_index1,
                                      UINT target_qubit_index2,
                                      __constant CTYPE matrix[16],
                                      CTYPE __global *state, ITYPE dim);
void double_qubit_dense_matrix_gate_nosimd(UINT target_qubit_index1,
                                           UINT target_qubit_index2,
                                           __constant CTYPE matrix[16],
                                           CTYPE __global *state, ITYPE dim);
void double_qubit_dense_matrix_gate_simd(UINT target_qubit_index1,
                                         UINT target_qubit_index2,
                                         __constant CTYPE matrix[16],
                                         CTYPE __global *state, ITYPE dim);
void double_qubit_dense_matrix_gate_sve(UINT target_qubit_index1,
                                        UINT target_qubit_index2,
                                        __constant CTYPE matrix[16],
                                        CTYPE __global *state, ITYPE dim);

/**
 * \~english
 * Apply a multi-qubit arbitrary gate.
 *
 * Apply a multi-qubit arbitrary gate defined by a dense matrix as a
 * one-dimentional array matrix[].
 * @param[in] target_qubit_index_list list of the target qubits
 * @param[in] target_qubit_index_count the number of the target qubits
 * @param[in] matrix description of a multi-qubit gate as a one-dimensional
 * array
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 任意の複数量子ビット演算を作用させて状態を更新。
 *
 * 任意の複数量子ビット演算を作用させて状態を更新。演算は、その行列成分を１次元配列として
 * matrix[] で与える。(j,k)成分は、matrix[dim*j+k]に対応。
 *
 * 例）パウリX演算子：{0, 1, 1,
 * 0}、アダマール演算子：{1/sqrt(2),1/sqrt(2),1/sqrt(2),-1/sqrt(2)}、
 *
 * CNOT演算:
 *
 * {1,0,0,0,
 *
 *  0,1,0,0,
 *
 *  0,0,0,1,
 *
 *  0,0,1,0}
 *
 * @param[in] target_qubit_index_list 作用する量子ビットのリスト
 * @param[in] target_qubit_index_count 作用する量子ビットの数
 * @param[in] matrix 複数量子ビット演算を定義する長さ 2^(2*
 * target_qubit_index_count) の一次元配列。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_dense_matrix_gate(UINT __global *target_qubit_index_list,
                                   UINT target_qubit_index_count,
                                   __constant CTYPE *matrix,
                                   CTYPE __global *state, ITYPE dim, HEAP heap);
void multi_qubit_dense_matrix_gate_parallel(
    UINT __global *target_qubit_index_list, UINT target_qubit_index_count,
    __constant CTYPE *matrix, CTYPE __global *state, ITYPE dim, HEAP heap);

/**
 * \~english
 * Apply a single-qubit controlled multi-qubit gate.
 *
 * Apply a single-qubit controlled multi-qubit gate. The multi-qubit gate is by
 * a dense matrix as a one-dimentional array matrix[].
 * @param[in] control_qubit_index index of the control qubit
 * @param[in] control_value value of the control qubit
 * @param[in] target_qubit_index_list list of the target qubits
 * @param[in] target_qubit_index_count the number of the target qubits
 * @param[in] matrix description of a multi-qubit gate as a one-dimensional
 * array
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * １つの制御量子ビットによる任意の複数量子ビット演算の制御演算を作用させて状態を更新。
 *
 *  １つの制御量子ビットによる任意の複数量子ビット演算の制御演算を作用させて状態を更新。制御量子ビットが
 * 0 もしくは 1 のどちらの場合に作用するかは control_value
 * によって指定。複数量子ビット演算は、その行列成分を１次元配列として matrix[]
 * で与える。(j,k)成分は、matrix[dim*j+k]に対応。
 *
 * 例）パウリX演算子：{0, 1, 1,
 * 0}、アダマール演算子：{1/sqrt(2),1/sqrt(2),1/sqrt(2),-1/sqrt(2)}、
 *
 * CNOT演算:
 *
 * {1,0,0,0,
 *
 *  0,1,0,0,
 *
 *  0,0,0,1,
 *
 *  0,0,1,0}
 *
 * @param[in] control_qubit_index 制御量子ビットのインデックス
 * @param[in] control_value 制御量子ビットの値
 * @param[in] target_qubit_index_list 作用する量子ビットのリスト
 * @param[in] target_qubit_index_count 作用する量子ビットの数
 * @param[in] matrix 複数量子ビット演算を定義する長さ 2^(2*
 * target_qubit_index_count) の一次元配列。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void single_qubit_control_multi_qubit_dense_matrix_gate(
    UINT control_qubit_index, UINT control_value,
    UINT __global *target_qubit_index_list, UINT target_qubit_index_count,
    CTYPE matrix[16], CTYPE __global *state, ITYPE dim, HEAP heap);

/**
 * \~english
 * Apply a multi-qubit controlled multi-qubit gate.
 *
 * Apply a multi-qubit controlled multi-qubit gate. The multi-qubit gate is by a
 * dense matrix as a one-dimentional array matrix[].
 *
 * @param[in] control_qubit_index_list list of the indexes of the control qubits
 * @param[in] control_value_list list of the vlues of the control qubits
 * @param[in] target_qubit_index_list list of the target qubits
 * @param[in] target_qubit_index_count the number of the target qubits
 * @param[in] matrix description of a multi-qubit gate as a one-dimensional
 * array
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * 複数の制御量子ビットによる任意の複数量子ビット演算の制御演算を作用させて状態を更新。
 *
 *  複数の制御量子ビットによる任意の複数量子ビット演算の制御演算を作用させて状態を更新。制御量子ビットは、control_qubit_index_listで指定され、制御演算が作用する
 * control_qubit_index_count 個の制御量子ビットの値は、 control_value_list
 * によって指定。複数量子ビット演算は、その行列成分を１次元配列として matrix[]
 * で与える。(j,k)成分は、matrix[dim*j+k]に対応。
 *
 * 例）パウリX演算子：{0, 1, 1,
 * 0}、アダマール演算子：{1/sqrt(2),1/sqrt(2),1/sqrt(2),-1/sqrt(2)}、
 *
 * CNOT演算:
 *
 * {1,0,0,0,
 *
 *  0,1,0,0,
 *
 *  0,0,0,1,
 *
 *  0,0,1,0}
 *
 * @param[in] control_qubit_index_list 制御量子ビットのインデックスのリスト
 * @param[in] control_value_list 制御演算が作用する制御量子ビットの値のリスト
 * @param[in] control_qubit_index_count 制御量子ビットの数
 * @param[in] target_qubit_index_list ターゲット量子ビットのリスト
 * @param[in] target_qubit_index_count ターゲット量子ビットの数
 * @param[in] matrix 複数量子ビット演算を定義する長さ
 * 2^(2*target_qubit_index_count) の一次元配列。
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */
void multi_qubit_control_multi_qubit_dense_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT __global *target_qubit_index_list,
    UINT target_qubit_index_count, CTYPE __global *matrix,
    CTYPE __global *state, ITYPE dim, HEAP heap);

/**
 * Diagonal gate
 **/
void multi_qubit_diagonal_matrix_gate(UINT __global *target_qubit_index_list,
                                      UINT target_qubit_index_count,
                                      CTYPE __global *diagonal_element,
                                      CTYPE __global *state, ITYPE dim,
                                      HEAP heap);

void multi_qubit_control_multi_qubit_diagonal_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT __global *target_qubit_index_list,
    UINT target_qubit_index_count, CTYPE __global *diagonal_element,
    CTYPE __global *state, ITYPE dim, HEAP heap);

/**
 * \~english
 * Reflect state according to another given state.
 *
 * Reflect state according to another given state. When reflect quantum state
 * |a> to state |s>, unitary operator give by 2|s><s|-I is applied to |a>.
 *
 * @param[in] reflection_state quantum state to characterize reflection unitary
 * operator
 * @param[in,out] state quantum state to update
 * @param[in] dim dimension
 *
 *
 * \~japanese-en
 * reflection gateを作用する。
 *
 * 与えられた状態にreflection
 * gateを作用する。量子状態|a>を量子状態|s>に関してreflectするとは、量子状態|a>に対して2|s><s|-Iというユニタリ操作をすることに対応する。
 *
 * @param[in] reflection_state reflection gateのユニタリを特徴づける量子状態
 * @param[in,out] state quantum 更新する量子状態
 * @param[in] dim 次元
 *
 */
void reflection_gate(CTYPE __global *reflection_state, CTYPE __global *state,
                     ITYPE dim);

void state_add(CTYPE __global *state_added, CTYPE __global *state, ITYPE dim);
void state_add_with_coef(CTYPE coef, CTYPE __global *state_added,
                         CTYPE __global *state, ITYPE dim);
void state_add_with_coef_single_thread(CTYPE coef, CTYPE __global *state_added,
                                       CTYPE __global *state, ITYPE dim);
void state_multiply(CTYPE coef, CTYPE __global *state, ITYPE dim);

////////////////////////////////

/// QFT
/**
 * update quantum state with one-qubit controlled dense matrix gate
 *
 * update quantum state with one-qubit controlled dense matrix gate
 * @param[in] angle rotation angle
 * @param[in] c_bit index of control qubit
 * @param[in] t_bit index of target qubit
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 */
void CUz_gate(double angle, UINT c_bit, UINT t_bit, CTYPE __global *state,
              ITYPE dim);
/**
 * update quantum state with large dense matrix
 *
 * update quantum state with large dense matrix
 * @param[in] k ???
 * @param[in] Nbits the number of qubits
 * @param[in] doSWAP ???
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 */
void qft(UINT k, UINT Nbits, int doSWAP, CTYPE __global *state, ITYPE dim);
/**
 * update quantum state with large dense matrix
 *
 * update quantum state with large dense matrix
 * @param[in] k ???
 * @param[in] Nbits the number of qubits
 * @param[in] doSWAP ???
 * @param[in,out] state quantum state
 * @param[in] dim dimension
 */
void inverse_qft(UINT k, UINT Nbits, int doSWAP, CTYPE __global *state,
                 ITYPE dim);

/*
 rules for arguments of update functions:
   The order of arguments must be
      information_about_applying_qubits -> Pauli_operator ->
 information_about_applying_gate_matrix_elements -> a_kind_of_rotation_angle ->
 state_vector -> dimension If there is control-qubit and target-qubit,
 control-qubit is the first argument. If an array of which the size is not known
 comes, the size of that array follows.

  Definition of update function is divided to named_gates,
 single_target_qubit_gates, multiple_target_qubit_gates, QFT_gates
 */

void dm_normalize(double squared_norm, CTYPE __global *state, ITYPE dim);
void dm_single_qubit_dense_matrix_gate_non_constant(UINT target_qubit_index,
                                                    CTYPE matrix[4],
                                                    CTYPE __global *state,
                                                    ITYPE dim);
void dm_single_qubit_dense_matrix_gate(UINT target_qubit_index,
                                       __constant CTYPE matrix[4],
                                       CTYPE __global *state, ITYPE dim);
void dm_multi_qubit_control_single_qubit_dense_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT target_qubit_index,
    __constant CTYPE matrix[4], CTYPE __global *state, ITYPE dim, HEAP heap);
void dm_multi_qubit_dense_matrix_gate(UINT __global *target_qubit_index_list,
                                      UINT target_qubit_index_count,
                                      CTYPE __global *matrix,
                                      CTYPE __global *state, ITYPE dim,
                                      HEAP heap);
void dm_multi_qubit_control_multi_qubit_dense_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT __global *target_qubit_index_list,
    UINT target_qubit_index_count, CTYPE __global *matrix,
    CTYPE __global *state, ITYPE dim, HEAP heap);

void dm_X_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_Y_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_Z_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_S_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_Sdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_T_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_Tdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_sqrtX_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_sqrtXdag_gate(UINT target_qubit_index, CTYPE __global *state,
                      ITYPE dim);
void dm_sqrtY_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_sqrtYdag_gate(UINT target_qubit_index, CTYPE __global *state,
                      ITYPE dim);
void dm_H_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_CNOT_gate(UINT control_qubit_index, UINT target_qubit_index,
                  CTYPE __global *state, ITYPE dim, HEAP heap);
void dm_CZ_gate(UINT control_qubit_index, UINT target_qubit_index,
                CTYPE __global *state, ITYPE dim, HEAP heap);
void dm_SWAP_gate(UINT target_qubit_index_0, UINT target_qubit_index_1,
                  CTYPE __global *state, ITYPE dim, HEAP heap);
void dm_P0_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_P1_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim);
void dm_RX_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
                ITYPE dim);
void dm_RY_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
                ITYPE dim);
void dm_RZ_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
                ITYPE dim);
void dm_multi_qubit_Pauli_gate_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim, HEAP heap);
void dm_multi_qubit_Pauli_rotation_gate_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, double angle, CTYPE __global *state,
    ITYPE dim, HEAP heap);

void double_qubit_dense_matrix_gate(UINT target_qubit_index0,
                                    UINT target_qubit_index1,
                                    __constant CTYPE matrix[16],
                                    CTYPE __global *state, ITYPE dim);
void double_qubit_dense_matrix_gate4cd(UINT target_qubit_index0,
                                       UINT target_qubit_index1,
                                       Matrix4cd eigen_matrix,
                                       CTYPE __global *state, ITYPE dim,
                                       HEAP heap);
void double_qubit_dense_matrix_gate_eigen(UINT target_qubit_index0,
                                          UINT target_qubit_index1,
                                          Matrix4cd eigen_matrix,
                                          CTYPE __global *state, ITYPE dim,
                                          HEAP heap);

// FALTAN
/*
void multi_qubit_dense_matrix_gate_eigen(
    UINT __global *target_qubit_index_list, UINT target_qubit_index_count,
    CTYPE __global* matrix, CTYPE __global* state, ITYPE dim);
void multi_qubit_dense_matrix_gate_eigen_xXcd(
    UINT __global *target_qubit_index_list, UINT target_qubit_index_count,
     Eigen::MatrixXcd& eigen_matrix, CTYPE __global* state, ITYPE dim);
void multi_qubit_dense_matrix_gate_eigen_matrix(
    UINT __global *target_qubit_index_list, UINT target_qubit_index_count,
     Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,
        Eigen::RowMajor>& eigen_matrix,
    CTYPE __global* state, ITYPE dim);

//FALTAN
void multi_qubit_sparse_matrix_gate_eigen(
    UINT __global *target_qubit_index_list, UINT target_qubit_index_count,
     Eigen::SparseMatrix<std::complex<double>>& eigen_matrix, CTYPE __global*
state, ITYPE dim);
*/

/**
 * \~english
 * Apply reversible boolean function as a unitary gate.
 *
 * Apply reversible boolean function as a unitary gate. Boolean function is
 * given as a pointer of int -> int function.
 *
 * @param[in] target_qubit_index_list ターゲット量子ビットのリスト
 * @param[in] target_qubit_index_count ターゲット量子ビットの数
 * @param[in] matrix 添え字および対象ビットの次元を受け取ると添え字を返す関数
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 *
 * \~japanese-en
 * 可逆回路関数をユニタリゲートとして作用する
 *
 *  可逆回路関数をユニタリゲートとして作用する。可逆回路関数は添え字を与えると結果の添え字を返す関数。
 *
 * @param[in] target_qubit_index_list ターゲット量子ビットのリスト
 * @param[in] target_qubit_index_count ターゲット量子ビットの数
 * @param[in] matrix 添え字および対象ビットの次元を受け取ると添え字を返す関数
 * @param[in,out] state 量子状態
 * @param[in] dim 次元
 *
 */

// FALTAN
/*
void reversible_boolean_gate(UINT __global *target_qubit_index_list,
    UINT target_qubit_index_count,
    std::function<ITYPE(ITYPE, ITYPE)> function_ptr, CTYPE __global* state,
ITYPE dim, HEAP heap);
*/

// elementary set of gates
#define c_zero                                                                 \
  { 0., 0. }
#define c_one                                                                  \
  { 1., 0. }
#define c_minus_one                                                            \
  { -1., 0. }
#define c_i                                                                    \
  { 0., 1. }
#define c_minus_i                                                              \
  { 0., -1. }

__constant CTYPE PAULI_MATRIX[4][4] = {{c_one, c_zero, c_zero, c_one},
                                       {c_zero, c_one, c_one, c_zero},
                                       {c_zero, c_minus_i, c_i, c_zero},
                                       {c_one, c_zero, c_zero, c_minus_one}};

__constant CTYPE S_GATE_MATRIX[4] = {c_one, c_zero, c_zero, c_i};

__constant CTYPE S_DAG_GATE_MATRIX[4] = {c_one, c_zero, c_zero, c_minus_i};

__constant CTYPE T_GATE_MATRIX[4] = {
    {COSPI8, -SINPI8}, c_zero, c_zero, {COSPI8, SINPI8}};

__constant CTYPE T_DAG_GATE_MATRIX[4] = {
    {COSPI8, SINPI8}, c_zero, c_zero, {COSPI8, -SINPI8}};

__constant CTYPE HADAMARD_MATRIX[4] = {
    {1. / SQRT2, 0}, {1. / SQRT2, 0}, {1. / SQRT2, 0}, {-1. / SQRT2, 0}};

__constant CTYPE SQRT_X_GATE_MATRIX[4] = {
    {0.5, 0.5}, {0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}};

__constant CTYPE SQRT_Y_GATE_MATRIX[4] = {
    {0.5, 0.5}, {-0.5, -0.5}, {0.5, 0.5}, {0.5, 0.5}};

__constant CTYPE SQRT_X_DAG_GATE_MATRIX[4] = {
    {0.5, -0.5}, {0.5, 0.5}, {0.5, 0.5}, {0.5, -0.5}};

__constant CTYPE SQRT_Y_DAG_GATE_MATRIX[4] = {
    {0.5, -0.5}, {0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}};

__constant CTYPE PROJ_0_MATRIX[4] = {c_one, c_zero, c_zero, c_zero};

__constant CTYPE PROJ_1_MATRIX[4] = {c_zero, c_zero, c_zero, c_one};

__constant CTYPE PHASE_90ROT[4] = {c_one, c_i, c_minus_one, c_minus_i};
__constant CTYPE PHASE_M90ROT[4] = {c_one, c_minus_i, c_minus_one, c_i};

void get_Pauli_masks_partial_list(UINT __global *target_qubit_index_list,
                                  UINT *Pauli_operator_type_list,
                                  UINT target_qubit_index_count,
                                  ITYPE *bit_flip_mask, ITYPE *phase_flip_mask,
                                  UINT *global_phase_90rot_count,
                                  UINT *pivot_qubit_index) {
  (*bit_flip_mask) = 0;
  (*phase_flip_mask) = 0;
  (*global_phase_90rot_count) = 0;
  (*pivot_qubit_index) = 0;
  for (UINT cursor = 0; cursor < target_qubit_index_count; ++cursor) {
    UINT target_qubit_index = target_qubit_index_list[cursor];
    switch (Pauli_operator_type_list[cursor]) {
    case 0: // I
      break;
    case 1: // X
      (*bit_flip_mask) ^= 1ULL << target_qubit_index;
      (*pivot_qubit_index) = target_qubit_index;
      break;
    case 2: // Y
      (*bit_flip_mask) ^= 1ULL << target_qubit_index;
      (*phase_flip_mask) ^= 1ULL << target_qubit_index;
      (*global_phase_90rot_count)++;
      (*pivot_qubit_index) = target_qubit_index;
      break;
    case 3: // Z
      (*phase_flip_mask) ^= 1ULL << target_qubit_index;
      break;
    default:
      printf("ERROR: Invalid Pauli operator ID called");
    }
  }
}

void get_Pauli_masks_whole_list(UINT *Pauli_operator_type_list,
                                UINT target_qubit_index_count,
                                ITYPE *bit_flip_mask, ITYPE *phase_flip_mask,
                                UINT *global_phase_90rot_count,
                                UINT *pivot_qubit_index) {
  (*bit_flip_mask) = 0;
  (*phase_flip_mask) = 0;
  (*global_phase_90rot_count) = 0;
  (*pivot_qubit_index) = 0;
  for (UINT target_qubit_index = 0;
       target_qubit_index < target_qubit_index_count; ++target_qubit_index) {
    switch (Pauli_operator_type_list[target_qubit_index]) {
    case 0: // I
      break;
    case 1: // X
      (*bit_flip_mask) ^= 1ULL << target_qubit_index;
      (*pivot_qubit_index) = target_qubit_index;
      break;
    case 2: // Y
      (*bit_flip_mask) ^= 1ULL << target_qubit_index;
      (*phase_flip_mask) ^= 1ULL << target_qubit_index;
      (*global_phase_90rot_count)++;
      (*pivot_qubit_index) = target_qubit_index;
      break;
    case 3: // Z
      (*phase_flip_mask) ^= 1ULL << target_qubit_index;
      break;
    default:
      printf("ERROR: Invalid Pauli operator ID called");
    }
  }
}

ITYPE __global *create_matrix_mask_list(UINT __global *qubit_index_list,
                                        UINT qubit_index_count, HEAP heap) {
  ITYPE matrix_dim = 1ULL << qubit_index_count;
  ITYPE __global *mask_list =
      (ITYPE __global *)calloc(heap, (size_t)matrix_dim, sizeof(ITYPE));
  ITYPE cursor = 0;

  for (cursor = 0; cursor < matrix_dim; ++cursor) {
    for (UINT bit_cursor = 0; bit_cursor < qubit_index_count; ++bit_cursor) {
      if ((cursor >> bit_cursor) % 2) {
        UINT bit_index = qubit_index_list[bit_cursor];
        mask_list[cursor] ^= (1ULL << bit_index);
      }
    }
  }
  return mask_list;
}

ITYPE
create_control_mask(UINT __global *qubit_index_list, UINT *value_list,
                    UINT size) {
  ITYPE mask = 0;
  for (UINT cursor = 0; cursor < size; ++cursor) {
    mask ^= (1ULL << qubit_index_list[cursor]) * value_list[cursor];
  }
  return mask;
}

static int compare_ui(UINT __global *a, UINT __global *b) {
  return (*a - *b);
}

// Podría también añadir left and right para controlar en caso de 
// querer ordenar solo una parte de un array
void sort_ui(UINT __global *v, int size, HEAP heap) {
    int i, j;
    bool swapped;
    for (i = 0; i < size - 1; i++) {
        swapped = false;
        for (j = 0; j < size - i - 1; j++) {
            if (compare_ui(&v[j], &v[j + 1]) > 0) {
                swap(&v[j], &v[j + 1], 1, heap);
                swapped = true;
            }
        }
        // If no two elements were swapped by inner loop,
        // then break
        if (swapped == false)
            break;
    }
}

UINT __global *create_sorted_ui_list(UINT __global *array, size_t size,
                                     HEAP heap) {
  UINT __global *new_array = (UINT __global *)calloc(heap, size, sizeof(UINT));
  memcpy(new_array, array, size * sizeof(UINT));
  sort_ui(new_array, size, heap);
  return new_array;
}
UINT __global *create_sorted_ui_list_value(UINT __global *array, size_t size,
                                           UINT value, HEAP heap) {
  UINT __global *new_array =
      (UINT __global *)calloc(heap, size + 1, sizeof(UINT));
  memcpy(new_array, array, size * sizeof(UINT));
  new_array[size] = value;
  sort_ui(new_array, size + 1, heap);
  return new_array;
}
UINT __global *create_sorted_ui_list_list(UINT __global *array1, size_t size1,
                                          UINT __global *array2, size_t size2,
                                          HEAP heap) {
  UINT __global *new_array =
      (UINT __global *)calloc(heap, size1 + size2, sizeof(UINT));
  memcpy(new_array, array1, size1 * sizeof(UINT));
  memcpy(new_array + size1, array2, size2 * sizeof(UINT));
  sort_ui(new_array, size1 + size2, heap);
  return new_array;
}

// memory allocation
CTYPE __global *allocate_quantum_state(ITYPE dim, HEAP heap) {
  CTYPE __global *state =
      (CTYPE __global *)malloc(heap, (size_t)(sizeof(CTYPE) * dim));
  // CTYPE __global* state =
  // (CTYPE*)_aligned_malloc((size_t)(sizeof(CTYPE)*dim), 32);

  if (!state) {
    printf("ERROR: Out of memory\n");
  }
  return state;
}

void release_quantum_state(CTYPE __global *state, HEAP heap) {
  free(heap, (uintptr_t)state);
  //_aligned_free(state);
}

CTYPE __global *dm_allocate_quantum_state(ITYPE dim, HEAP heap) {
  CTYPE __global *state =
      (CTYPE __global *)malloc(heap, (size_t)(sizeof(CTYPE) * dim * dim));
  if (!state) {
    printf("ERROR: Out of memory\n");
  }
  return state;
}

void dm_initialize_quantum_state(CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  for (index = 0; index < dim * dim; ++index) {
    state[index] = double_to_complex(0, 0);
  }
  state[0] = double_to_complex(1.0, 0);
}

void dm_release_quantum_state(CTYPE __global *state, HEAP heap) {
  free(heap, (uintptr_t)state);
}

void dm_initialize_with_pure_state(CTYPE __global *state,
                                   CTYPE __global *pure_state, ITYPE dim) {
  ITYPE ind_y;
  for (ind_y = 0; ind_y < dim; ++ind_y) {
    ITYPE ind_x;
    for (ind_x = 0; ind_x < dim; ++ind_x) {
      state[ind_y * dim + ind_x] =
          mul(pure_state[ind_y], conj(pure_state[ind_x]));
    }
  }
}

// calculate norm
double state_norm_squared(CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  double norm = 0;
  for (index = 0; index < dim; ++index) {
    norm += pow(_cabs(state[index]), 2);
  }
  return norm;
}

// calculate norm
double state_norm_squared_single_thread(CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  double norm = 0;
  for (index = 0; index < dim; ++index) {
    norm += pow(_cabs(state[index]), 2);
  }
  return norm;
}

// calculate inner product of two state vector
CTYPE
state_inner_product(CTYPE __global *state_bra, CTYPE __global *state_ket,
                    ITYPE dim) {
  double real_sum = 0.;
  double imag_sum = 0.;
  ITYPE index;
  for (index = 0; index < dim; ++index) {
    CTYPE value;
    value = add(value, mul(conj(state_bra[index]), state_ket[index]));
    real_sum += _creal(value);
    imag_sum += _cimag(value);
  }
  return double_to_complex(real_sum, imag_sum);
}

void state_tensor_product(CTYPE __global *state_left, ITYPE dim_left,
                          CTYPE __global *state_right, ITYPE dim_right,
                          CTYPE *state_dst) {
  ITYPE index_left, index_right;
  for (index_left = 0; index_left < dim_left; ++index_left) {
    CTYPE val_left = state_left[index_left];
    for (index_right = 0; index_right < dim_right; ++index_right) {
      state_dst[index_left * dim_right + index_right] =
          mul(val_left, state_right[index_right]);
    }
  }
}
void state_permutate_qubit(UINT *qubit_order, CTYPE __global *state_src,
                           CTYPE *state_dst, UINT qubit_count, ITYPE dim) {
  ITYPE index;
  for (index = 0; index < dim; ++index) {
    ITYPE src_index = 0;
    for (UINT qubit_index = 0; qubit_index < qubit_count; ++qubit_index) {
      if ((index >> qubit_index) % 2) {
        src_index += 1ULL << qubit_order[qubit_index];
      }
    }
    state_dst[index] = state_src[src_index];
  }
}

void state_drop_qubits(UINT __global *target, UINT *projection,
                       UINT target_count, CTYPE __global *state_src,
                       CTYPE *state_dst, ITYPE dim, HEAP heap) {
  ITYPE dst_dim = dim >> target_count;
  UINT __global *sorted_target =
      create_sorted_ui_list(target, target_count, heap);
  ITYPE projection_mask = 0;
  for (UINT target_index = 0; target_index < target_count; ++target_index) {
    projection_mask ^= (projection[target_index] << target[target_index]);
  }

  ITYPE index;
  for (index = 0; index < dst_dim; ++index) {
    ITYPE src_index = index;
    for (UINT target_index = 0; target_index < target_count; ++target_index) {
      UINT insert_index = sorted_target[target_index];
      src_index = insert_zero_to_basis_index(src_index, 1ULL << insert_index,
                                             insert_index);
    }
    src_index ^= projection_mask;
    state_dst[index] = state_src[src_index];
  }
  free(heap, (uintptr_t)sorted_target);
}

// calculate norm
double dm_state_norm_squared(CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  double norm = 0;
  for (index = 0; index < dim; ++index) {
    norm += _creal(state[index * dim + index]);
  }
  return norm;
}

// calculate entropy of probability distribution of Z-basis measurements
double dm_measurement_distribution_entropy(CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  double ent = 0;
  double eps = 1e-15;
  for (index = 0; index < dim; ++index) {
    double prob = _creal(state[index * dim + index]);
    if (prob > eps) {
      ent += -1.0 * prob * log(prob);
    }
  }
  return ent;
}

// calculate probability with which we obtain 0 at target qubit
double dm_M0_prob(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = 1ULL << target_qubit_index;
  ITYPE state_index;
  double sum = 0.;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_0 =
        insert_zero_to_basis_index(state_index, mask, target_qubit_index);
    sum += _creal(state[basis_0 * dim + basis_0]);
  }
  return sum;
}

// calculate probability with which we obtain 1 at target qubit
double dm_M1_prob(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = 1ULL << target_qubit_index;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_1 =
        insert_zero_to_basis_index(state_index, mask, target_qubit_index) ^
        mask;
    sum += _creal(state[basis_1 * dim + basis_1]);
  }

  return sum;
}

// calculate merginal probability with which we obtain the set of values
// measured_value_list at sorted_target_qubit_index_list warning:
// sorted_target_qubit_index_list must be sorted.
double dm_marginal_prob(UINT *sorted_target_qubit_index_list,
                        UINT *measured_value_list,
                        UINT target_qubit_index_count, CTYPE __global *state,
                        ITYPE dim) {
  ITYPE loop_dim = dim >> target_qubit_index_count;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis = state_index;
    for (UINT cursor = 0; cursor < target_qubit_index_count; cursor++) {
      UINT insert_index = sorted_target_qubit_index_list[cursor];
      ITYPE mask = 1ULL << insert_index;
      basis = insert_zero_to_basis_index(basis, mask, insert_index);
      basis ^= mask * measured_value_list[cursor];
    }
    sum += _creal(state[basis * dim + basis]);
  }

  return sum;
}

void dm_state_add(CTYPE __global *state_added, CTYPE __global *state,
                  ITYPE dim) {
  ITYPE index;

  for (index = 0; index < dim * dim; ++index) {
    state[index] = add(state[index], state_added[index]);
  }
}

void dm_state_add_with_coef(CTYPE coef, CTYPE __global *state_added,
                            CTYPE __global *state, ITYPE dim) {
  ITYPE index;

  for (index = 0; index < dim * dim; ++index) {
    state[index] = add(state[index], mul(coef, state_added[index]));
  }
}

void dm_state_multiply(CTYPE coef, CTYPE __global *state, ITYPE dim) {
  ITYPE index;

  for (index = 0; index < dim * dim; ++index) {
    state[index] = mul(state[index], coef);
  }
}

double dm_expectation_value_multi_qubit_Pauli_operator_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim) {
  CTYPE sum = double_to_complex(0, 0);
  for (ITYPE state_index = 0; state_index < dim; ++state_index) {
    CTYPE coef = double_to_complex(1.0, 0);
    ITYPE state_index_sub = state_index;
    for (UINT i = 0; i < target_qubit_index_count; ++i) {
      UINT pauli_type = Pauli_operator_type_list[i];
      UINT target_qubit_index = target_qubit_index_list[i];
      if (pauli_type == 1) {
        state_index_sub ^= ((1ULL) << target_qubit_index);
      } else if (pauli_type == 2) {
        coef = mul(coef, double_to_complex(0, 1.0));
        if (state_index_sub & ((1ULL) << target_qubit_index)) {
          coef = mul(coef, double_to_complex(-1.0, 0));
        }
        state_index_sub ^= ((1ULL) << target_qubit_index);
      } else if (pauli_type == 3) {
        if (state_index_sub & ((1ULL) << target_qubit_index)) {
          coef = mul(coef, double_to_complex(-1.0, 0));
        }
      }
    }
    sum = add(sum, mul(coef, state[state_index * dim + state_index_sub]));
  }
  return _creal(sum);
}

void dm_state_tensor_product(CTYPE __global *state_left, ITYPE dim_left,
                             CTYPE __global *state_right, ITYPE dim_right,
                             CTYPE *state_dst) {
  ITYPE y_left, x_left, y_right, x_right;
  ITYPE dim_new = dim_left * dim_right;
  for (y_left = 0; y_left < dim_left; ++y_left) {
    for (x_left = 0; x_left < dim_left; ++x_left) {
      CTYPE val_left = state_left[y_left * dim_left + x_left];
      for (y_right = 0; y_right < dim_right; ++y_right) {
        for (x_right = 0; x_right < dim_right; ++x_right) {
          CTYPE val_right = state_right[y_right * dim_right + x_right];
          ITYPE x_new = x_left * dim_right + x_right;
          ITYPE y_new = y_left * dim_right + y_right;
          state_dst[y_new * dim_new + x_new] = mul(val_right, val_left);
        }
      }
    }
  }
}

void dm_state_permutate_qubit(UINT *qubit_order, CTYPE __global *state_src,
                              CTYPE *state_dst, UINT qubit_count, ITYPE dim) {
  ITYPE y, x;
  for (y = 0; y < dim; ++y) {
    for (x = 0; x < dim; ++x) {
      ITYPE src_x = 0, src_y = 0;
      for (UINT qubit_index = 0; qubit_index < qubit_count; ++qubit_index) {
        if ((x >> qubit_index) % 2) {
          src_x += 1ULL << qubit_order[qubit_index];
        }
        if ((y >> qubit_index) % 2) {
          src_y += 1ULL << qubit_order[qubit_index];
        }
      }
      state_dst[y * dim + x] = state_src[src_y * dim + src_x];
    }
  }
}

void dm_state_partial_trace_from_density_matrix(UINT __global *target,
                                                UINT target_count,
                                                CTYPE __global *state_src,
                                                CTYPE *state_dst, ITYPE dim,
                                                HEAP heap) {
  ITYPE dst_dim = dim >> target_count;
  ITYPE trace_dim = 1ULL << target_count;
  UINT __global *sorted_target =
      create_sorted_ui_list(target, target_count, heap);
  ITYPE __global *mask_list =
      create_matrix_mask_list(target, target_count, heap);

  ITYPE y, x;
  for (y = 0; y < dst_dim; ++y) {
    for (x = 0; x < dst_dim; ++x) {
      ITYPE base_x = x;
      ITYPE base_y = y;
      for (UINT target_index = 0; target_index < target_count; ++target_index) {
        UINT insert_index = sorted_target[target_index];
        base_x = insert_zero_to_basis_index(base_x, 1ULL << insert_index,
                                            insert_index);
        base_y = insert_zero_to_basis_index(base_y, 1ULL << insert_index,
                                            insert_index);
      }
      CTYPE val = double_to_complex(0.0, 0);
      for (ITYPE idx = 0; idx < trace_dim; ++idx) {
        ITYPE src_x = base_x ^ mask_list[idx];
        ITYPE src_y = base_y ^ mask_list[idx];
        val = add(val, state_src[src_y * dim + src_x]);
      }
      state_dst[y * dst_dim + x] = val;
    }
  }
  free(heap, (uintptr_t)sorted_target);
  free(heap, (uintptr_t)mask_list);
}

void dm_state_partial_trace_from_state_vector(UINT __global *target,
                                              UINT target_count,
                                              CTYPE __global *state_src,
                                              CTYPE *state_dst, ITYPE dim,
                                              HEAP heap) {
  ITYPE dst_dim = dim >> target_count;
  ITYPE trace_dim = 1ULL << target_count;
  UINT __global *sorted_target =
      create_sorted_ui_list(target, target_count, heap);
  ITYPE __global *mask_list =
      create_matrix_mask_list(target, target_count, heap);

  ITYPE y, x;
  for (y = 0; y < dst_dim; ++y) {
    for (x = 0; x < dst_dim; ++x) {
      ITYPE base_x = x;
      ITYPE base_y = y;
      for (UINT target_index = 0; target_index < target_count; ++target_index) {
        UINT insert_index = sorted_target[target_index];
        base_x = insert_zero_to_basis_index(base_x, 1ULL << insert_index,
                                            insert_index);
        base_y = insert_zero_to_basis_index(base_y, 1ULL << insert_index,
                                            insert_index);
      }
      CTYPE val = double_to_complex(0.0, 0);
      for (ITYPE idx = 0; idx < trace_dim; ++idx) {
        ITYPE src_x = base_x ^ mask_list[idx];
        ITYPE src_y = base_y ^ mask_list[idx];
        val = add(val, mul(state_src[src_y], conj(state_src[src_x])));
      }
      state_dst[y * dst_dim + x] = val;
    }
  }
  // AQUÍ HAY UN ERROR DE SEGFAULT (con los dos)
  free(heap, (uintptr_t)sorted_target);
  free(heap, (uintptr_t)mask_list);
}

void dm_normalize(double squared_norm, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim;
  double normalize_factor = 1. / squared_norm;
  ITYPE state_index_y;

  for (state_index_y = 0; state_index_y < loop_dim; ++state_index_y) {
    ITYPE state_index_x;
    for (state_index_x = 0; state_index_x < loop_dim; ++state_index_x) {
      state[state_index_y * dim + state_index_x] =
          mul(state[state_index_y * dim + state_index_x],
              double_to_complex(normalize_factor, 0));
    }
  }
}

void dm_single_qubit_dense_matrix_gate_non_constant(UINT target_qubit_index,
                                                    CTYPE matrix[4],
                                                    CTYPE __global *state,
                                                    ITYPE dim) {
  // target mask
  ITYPE target_mask = 1ULL << target_qubit_index;

  // loop variables
  ITYPE loop_dim = dim / 2;

  // create extended matrix
  CTYPE ext_matrix[16];
  for (int y = 0; y < 4; ++y) {
    int y1 = y / 2;
    int y2 = y % 2;
    for (int x = 0; x < 4; ++x) {
      int x1 = x / 2;
      int x2 = x % 2;
      ext_matrix[y * 4 + x] =
          mul(matrix[y1 * 2 + x1], conj(matrix[y2 * 2 + x2]));
    }
  }

  ITYPE state_index_x, state_index_y;

  for (state_index_y = 0; state_index_y < loop_dim; ++state_index_y) {
    // create vertical index
    ITYPE basis_0_y = insert_zero_to_basis_index(state_index_y, target_mask,
                                                 target_qubit_index);
    // flip target bit
    ITYPE basis_1_y = basis_0_y ^ target_mask;

    for (state_index_x = 0; state_index_x < loop_dim; ++state_index_x) {
      // create vertical index
      ITYPE basis_0_x = insert_zero_to_basis_index(state_index_x, target_mask,
                                                   target_qubit_index);
      // flip target bit
      ITYPE basis_1_x = basis_0_x ^ target_mask;

      ITYPE basis_00 = basis_0_y * dim + basis_0_x;
      ITYPE basis_01 = basis_0_y * dim + basis_1_x;
      ITYPE basis_10 = basis_1_y * dim + basis_0_x;
      ITYPE basis_11 = basis_1_y * dim + basis_1_x;

      // fetch values
      CTYPE cval_00 = state[basis_00];
      CTYPE cval_01 = state[basis_01];
      CTYPE cval_10 = state[basis_10];
      CTYPE cval_11 = state[basis_11];

      // set values
      state[basis_00] =
          add(add(add(mul(ext_matrix[0], cval_00), mul(ext_matrix[1], cval_01)),
                  mul(ext_matrix[2], cval_10)),
              mul(ext_matrix[3], cval_11));
      state[basis_01] =
          add(add(add(mul(ext_matrix[4], cval_00), mul(ext_matrix[5], cval_01)),
                  mul(ext_matrix[6], cval_10)),
              mul(ext_matrix[7], cval_11));
      state[basis_10] =
          add(add(add(mul(ext_matrix[8], cval_00), mul(ext_matrix[9], cval_01)),
                  mul(ext_matrix[10], cval_10)),
              mul(ext_matrix[11], cval_11));
      state[basis_11] = add(
          add(add(mul(ext_matrix[12], cval_00), mul(ext_matrix[13], cval_01)),
              mul(ext_matrix[14], cval_10)),
          mul(ext_matrix[15], cval_11));
    }
  }
}

void dm_single_qubit_dense_matrix_gate(UINT target_qubit_index,
                                       __constant CTYPE matrix[4],
                                       CTYPE __global *state, ITYPE dim) {
  // target mask
  ITYPE target_mask = 1ULL << target_qubit_index;

  // loop variables
  ITYPE loop_dim = dim / 2;

  // create extended matrix
  CTYPE ext_matrix[16];
  for (int y = 0; y < 4; ++y) {
    int y1 = y / 2;
    int y2 = y % 2;
    for (int x = 0; x < 4; ++x) {
      int x1 = x / 2;
      int x2 = x % 2;
      ext_matrix[y * 4 + x] =
          mul(matrix[y1 * 2 + x1], conj(matrix[y2 * 2 + x2]));
    }
  }

  ITYPE state_index_x, state_index_y;

  for (state_index_y = 0; state_index_y < loop_dim; ++state_index_y) {
    // create vertical index
    ITYPE basis_0_y = insert_zero_to_basis_index(state_index_y, target_mask,
                                                 target_qubit_index);
    // flip target bit
    ITYPE basis_1_y = basis_0_y ^ target_mask;

    for (state_index_x = 0; state_index_x < loop_dim; ++state_index_x) {
      // create vertical index
      ITYPE basis_0_x = insert_zero_to_basis_index(state_index_x, target_mask,
                                                   target_qubit_index);
      // flip target bit
      ITYPE basis_1_x = basis_0_x ^ target_mask;

      ITYPE basis_00 = basis_0_y * dim + basis_0_x;
      ITYPE basis_01 = basis_0_y * dim + basis_1_x;
      ITYPE basis_10 = basis_1_y * dim + basis_0_x;
      ITYPE basis_11 = basis_1_y * dim + basis_1_x;

      // fetch values
      CTYPE cval_00 = state[basis_00];
      CTYPE cval_01 = state[basis_01];
      CTYPE cval_10 = state[basis_10];
      CTYPE cval_11 = state[basis_11];

      // set values
      state[basis_00] =
          add(add(add(mul(ext_matrix[0], cval_00), mul(ext_matrix[1], cval_01)),
                  mul(ext_matrix[2], cval_10)),
              mul(ext_matrix[3], cval_11));
      state[basis_01] =
          add(add(add(mul(ext_matrix[4], cval_00), mul(ext_matrix[5], cval_01)),
                  mul(ext_matrix[6], cval_10)),
              mul(ext_matrix[7], cval_11));
      state[basis_10] =
          add(add(add(mul(ext_matrix[8], cval_00), mul(ext_matrix[9], cval_01)),
                  mul(ext_matrix[10], cval_10)),
              mul(ext_matrix[11], cval_11));
      state[basis_11] = add(
          add(add(mul(ext_matrix[12], cval_00), mul(ext_matrix[13], cval_01)),
              mul(ext_matrix[14], cval_10)),
          mul(ext_matrix[15], cval_11));
    }
  }
}

void dm_multi_qubit_control_single_qubit_dense_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT target_qubit_index,
    __constant CTYPE matrix[4], CTYPE __global *state, ITYPE dim, HEAP heap) {
  // insert index list
  UINT insert_index_list_count = control_qubit_index_count + 1;
  UINT __global *insert_index_list = create_sorted_ui_list_value(
      control_qubit_index_list, control_qubit_index_count, target_qubit_index,
      heap);

  // target mask
  ITYPE target_mask = 1ULL << target_qubit_index;

  // control mask
  ITYPE control_mask = create_control_mask(
      control_qubit_index_list, control_value_list, control_qubit_index_count);

  // loop variables
  ITYPE loop_dim = dim >> insert_index_list_count;

  CTYPE adjoint_matrix[4];
  adjoint_matrix[0] = conj(matrix[0]);
  adjoint_matrix[1] = conj(matrix[2]);
  adjoint_matrix[2] = conj(matrix[1]);
  adjoint_matrix[3] = conj(matrix[3]);

  ITYPE state_index_x, state_index_y;

  for (state_index_x = 0; state_index_x < dim; ++state_index_x) {
    for (state_index_y = 0; state_index_y < loop_dim; ++state_index_y) {
      // create base index
      ITYPE basis_c_t0_y = state_index_y;
      for (UINT cursor = 0; cursor < insert_index_list_count; ++cursor) {
        basis_c_t0_y = insert_zero_to_basis_index(
            basis_c_t0_y, 1ULL << insert_index_list[cursor],
            insert_index_list[cursor]);
      }

      // flip controls
      basis_c_t0_y ^= control_mask;

      // gather target
      ITYPE basis_c_t1_y = basis_c_t0_y ^ target_mask;

      // set index
      ITYPE basis_0 = basis_c_t0_y * dim + state_index_x;
      ITYPE basis_1 = basis_c_t1_y * dim + state_index_x;

      // fetch values
      CTYPE cval_0 = state[basis_0];
      CTYPE cval_1 = state[basis_1];

      // set values
      state[basis_0] = add(mul(matrix[0], cval_0), mul(matrix[1], cval_1));
      state[basis_1] = add(mul(matrix[2], cval_0), mul(matrix[3], cval_1));
    }
  }

  for (state_index_y = 0; state_index_y < dim; ++state_index_y) {
    for (state_index_x = 0; state_index_x < loop_dim; ++state_index_x) {
      // create base index
      ITYPE basis_c_t0_x = state_index_x;
      for (UINT cursor = 0; cursor < insert_index_list_count; ++cursor) {
        basis_c_t0_x = insert_zero_to_basis_index(
            basis_c_t0_x, 1ULL << insert_index_list[cursor],
            insert_index_list[cursor]);
      }

      // flip controls
      basis_c_t0_x ^= control_mask;

      // gather target
      ITYPE basis_c_t1_x = basis_c_t0_x ^ target_mask;

      // set index
      ITYPE basis_0 = state_index_y * dim + basis_c_t0_x;
      ITYPE basis_1 = state_index_y * dim + basis_c_t1_x;

      // fetch values
      CTYPE cval_0 = state[basis_0];
      CTYPE cval_1 = state[basis_1];

      // set values
      state[basis_0] =
          add(mul(cval_0, adjoint_matrix[0]), mul(cval_1, adjoint_matrix[2]));
      state[basis_1] =
          add(mul(cval_0, adjoint_matrix[1]), mul(cval_1, adjoint_matrix[3]));
    }
  }

  free(heap, (uintptr_t)insert_index_list);
}

void dm_multi_qubit_dense_matrix_gate(UINT __global *target_qubit_index_list,
                                      UINT target_qubit_index_count,
                                      CTYPE __global *matrix,
                                      CTYPE __global *state, ITYPE dim,
                                      HEAP heap) {
  // matrix dim, mask, buffer
  ITYPE matrix_dim = 1ULL << target_qubit_index_count;
  ITYPE __global *matrix_mask_list = create_matrix_mask_list(
      target_qubit_index_list, target_qubit_index_count, heap);

  // create extended matrix
  CTYPE __global *adjoint_matrix = (CTYPE __global *)malloc(
      heap, (size_t)(sizeof(CTYPE) * (matrix_dim * matrix_dim)));
  for (ITYPE y = 0; y < matrix_dim; ++y) {
    for (ITYPE x = 0; x < matrix_dim; ++x) {
      adjoint_matrix[y * matrix_dim + x] = conj(matrix[x * matrix_dim + y]);
    }
  }

  // insert index
  UINT __global *sorted_insert_index_list = create_sorted_ui_list(
      target_qubit_index_list, target_qubit_index_count, heap);

  // loop variables
  ITYPE loop_dim = dim >> target_qubit_index_count;

#ifndef _OPENMP
  CTYPE __global *buffer = (CTYPE __global *)malloc(
      heap, (size_t)(sizeof(CTYPE) * matrix_dim * matrix_dim));
  ITYPE state_index_y;
  for (state_index_y = 0; state_index_y < loop_dim; ++state_index_y) {
    // create base index
    ITYPE basis_0_y = state_index_y;
    for (UINT cursor = 0; cursor < target_qubit_index_count; cursor++) {
      UINT insert_index = sorted_insert_index_list[cursor];
      basis_0_y = insert_zero_to_basis_index(basis_0_y, 1ULL << insert_index,
                                             insert_index);
    }

    ITYPE state_index_x;
    for (state_index_x = 0; state_index_x < loop_dim; ++state_index_x) {
      // create base index
      ITYPE basis_0_x = state_index_x;
      for (UINT cursor = 0; cursor < target_qubit_index_count; cursor++) {
        UINT insert_index = sorted_insert_index_list[cursor];
        basis_0_x = insert_zero_to_basis_index(basis_0_x, 1ULL << insert_index,
                                               insert_index);
      }

      // compute matrix-matrix multiply
      // TODO: improve matmul
      for (ITYPE y = 0; y < matrix_dim; ++y) {
        for (ITYPE x = 0; x < matrix_dim; ++x) {
          buffer[y * matrix_dim + x] = double_to_complex(0, 0);
          for (ITYPE k = 0; k < matrix_dim; ++k) {
            ITYPE dm_index_x = basis_0_x ^ matrix_mask_list[x];
            ITYPE dm_index_k = basis_0_y ^ matrix_mask_list[k];
            buffer[y * matrix_dim + x] =
                add(buffer[y * matrix_dim + x],
                    mul(matrix[y * matrix_dim + k],
                        state[dm_index_k * dim + dm_index_x]));
          }
        }
      }

      for (ITYPE y = 0; y < matrix_dim; ++y) {
        for (ITYPE x = 0; x < matrix_dim; ++x) {
          ITYPE dm_index_x = basis_0_x ^ matrix_mask_list[x];
          ITYPE dm_index_y = basis_0_y ^ matrix_mask_list[y];
          ITYPE dm_index = dm_index_y * dim + dm_index_x;
          state[dm_index] = double_to_complex(0, 0);
          for (ITYPE k = 0; k < matrix_dim; ++k) {
            state[dm_index] =
                add(state[dm_index], mul(buffer[y * matrix_dim + k],
                                         adjoint_matrix[k * matrix_dim + x]));
          }
        }
      }
    }
  }
  free(heap, (uintptr_t)buffer);
#else
  OMPutil::get_inst().set_qulacs_num_threads(dim, 0);
  UINT thread_count = omp_get_max_threads();
  CTYPE *buffer_list = (CTYPE *)malloc(
      heap, (size_t)(sizeof(CTYPE) * matrix_dim * matrix_dim * thread_count));

  ITYPE block_size = loop_dim / thread_count;
  ITYPE residual = loop_dim % thread_count;

  {
    UINT thread_id = omp_get_thread_num();
    ITYPE start_index =
        block_size * thread_id + (residual > thread_id ? thread_id : residual);
    ITYPE end_index = block_size * (thread_id + 1) +
                      (residual > (thread_id + 1) ? (thread_id + 1) : residual);
    CTYPE *buffer = buffer_list + thread_id * matrix_dim * matrix_dim;

    ITYPE state_index_y;
    for (state_index_y = start_index; state_index_y < end_index;
         ++state_index_y) {
      // create base index
      ITYPE basis_0_y = state_index_y;
      for (UINT cursor = 0; cursor < target_qubit_index_count; cursor++) {
        UINT insert_index = sorted_insert_index_list[cursor];
        basis_0_y = insert_zero_to_basis_index(basis_0_y, 1ULL << insert_index,
                                               insert_index);
      }

      ITYPE state_index_x;
      for (state_index_x = 0; state_index_x < loop_dim; ++state_index_x) {
        // create base index
        ITYPE basis_0_x = state_index_x;
        for (UINT cursor = 0; cursor < target_qubit_index_count; cursor++) {
          UINT insert_index = sorted_insert_index_list[cursor];
          basis_0_x = insert_zero_to_basis_index(
              basis_0_x, 1ULL << insert_index, insert_index);
        }

        // compute matrix-matrix multiply
        // TODO: improve matmul
        for (ITYPE y = 0; y < matrix_dim; ++y) {
          for (ITYPE x = 0; x < matrix_dim; ++x) {
            buffer[y * matrix_dim + x] = 0;
            for (ITYPE k = 0; k < matrix_dim; ++k) {
              ITYPE dm_index_x = basis_0_x ^ matrix_mask_list[x];
              ITYPE dm_index_k = basis_0_y ^ matrix_mask_list[k];
              buffer[y * matrix_dim + x] +=
                  matrix[y * matrix_dim + k] *
                  state[dm_index_k * dim + dm_index_x];
            }
          }
        }

        for (ITYPE y = 0; y < matrix_dim; ++y) {
          for (ITYPE x = 0; x < matrix_dim; ++x) {
            ITYPE dm_index_x = basis_0_x ^ matrix_mask_list[x];
            ITYPE dm_index_y = basis_0_y ^ matrix_mask_list[y];
            ITYPE dm_index = dm_index_y * dim + dm_index_x;
            state[dm_index] = 0;
            for (ITYPE k = 0; k < matrix_dim; ++k) {
              state[dm_index] += buffer[y * matrix_dim + k] *
                                 adjoint_matrix[k * matrix_dim + x];
            }
          }
        }
      }
    }
  }
  OMPutil::get_inst().reset_qulacs_num_threads();
  free(heap, (uintptr_t)buffer_list);
#endif
  free(heap, (uintptr_t)adjoint_matrix);
  free(heap, (uintptr_t)sorted_insert_index_list);
  free(heap, (uintptr_t)matrix_mask_list);
}

void dm_multi_qubit_control_multi_qubit_dense_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT __global *target_qubit_index_list,
    UINT target_qubit_index_count, CTYPE __global *matrix,
    CTYPE __global *state, ITYPE dim, HEAP heap) {
  // matrix dim, mask, buffer
  ITYPE matrix_dim = 1ULL << target_qubit_index_count;
  ITYPE __global *matrix_mask_list = create_matrix_mask_list(
      target_qubit_index_list, target_qubit_index_count, heap);

  // insert index
  UINT insert_index_count =
      target_qubit_index_count + control_qubit_index_count;
  UINT __global *sorted_insert_index_list = create_sorted_ui_list_list(
      target_qubit_index_list, target_qubit_index_count,
      control_qubit_index_list, control_qubit_index_count, heap);

  // control mask
  ITYPE control_mask = create_control_mask(
      control_qubit_index_list, control_value_list, control_qubit_index_count);

  // loop varaibles
  ITYPE loop_dim =
      dim >> (target_qubit_index_count + control_qubit_index_count);

  CTYPE __global *adjoint_matrix = (CTYPE __global *)malloc(
      heap, (size_t)(sizeof(CTYPE) * matrix_dim * matrix_dim));
  for (ITYPE y = 0; y < matrix_dim; ++y) {
    for (ITYPE x = 0; x < matrix_dim; ++x) {
      adjoint_matrix[y * matrix_dim + x] = conj(matrix[x * matrix_dim + y]);
    }
  }

#ifndef _OPENMP
  CTYPE __global *buffer =
      (CTYPE __global *)malloc(heap, (size_t)(sizeof(CTYPE) * matrix_dim));
  ITYPE state_index_x, state_index_y;
  for (state_index_x = 0; state_index_x < dim; ++state_index_x) {
    for (state_index_y = 0; state_index_y < loop_dim; ++state_index_y) {
      // create base index
      ITYPE basis_0_y = state_index_y;
      for (UINT cursor = 0; cursor < insert_index_count; cursor++) {
        UINT insert_index = sorted_insert_index_list[cursor];
        basis_0_y = insert_zero_to_basis_index(basis_0_y, 1ULL << insert_index,
                                               insert_index);
      }

      // flip control masks
      basis_0_y ^= control_mask;

      // compute matrix vector mul
      for (ITYPE y = 0; y < matrix_dim; ++y) {
        buffer[y] = double_to_complex(0, 0);
        for (ITYPE x = 0; x < matrix_dim; ++x) {
          ITYPE dm_index_y = basis_0_y ^ matrix_mask_list[x];
          buffer[y] =
              add(buffer[y], mul(matrix[y * matrix_dim + x],
                                 state[dm_index_y * dim + state_index_x]));
        }
      }

      // set result
      for (ITYPE y = 0; y < matrix_dim; ++y) {
        ITYPE dm_index_y = basis_0_y ^ matrix_mask_list[y];
        state[dm_index_y * dim + state_index_x] = buffer[y];
      }
    }
  }
  for (state_index_y = 0; state_index_y < dim; ++state_index_y) {
    for (state_index_x = 0; state_index_x < loop_dim; ++state_index_x) {
      // create base index
      ITYPE basis_0_x = state_index_x;
      for (UINT cursor = 0; cursor < insert_index_count; cursor++) {
        UINT insert_index = sorted_insert_index_list[cursor];
        basis_0_x = insert_zero_to_basis_index(basis_0_x, 1ULL << insert_index,
                                               insert_index);
      }

      // flip control masks
      basis_0_x ^= control_mask;

      // compute matrix vector mul
      for (ITYPE y = 0; y < matrix_dim; ++y) {
        buffer[y] = double_to_complex(0, 0);
        for (ITYPE x = 0; x < matrix_dim; ++x) {
          ITYPE dm_index_x = basis_0_x ^ matrix_mask_list[x];
          buffer[y] =
              add(buffer[y], mul(state[state_index_y * dim + dm_index_x],
                                 adjoint_matrix[x * matrix_dim + y]));
        }
      }

      // set result
      for (ITYPE y = 0; y < matrix_dim; ++y) {
        ITYPE dm_index_x = basis_0_x ^ matrix_mask_list[y];
        state[state_index_y * dim + dm_index_x] = buffer[y];
      }
    }
  }
  free(heap, (uintptr_t)buffer);
#else
  OMPutil::get_inst().set_qulacs_num_threads(dim, 0);
  UINT thread_count = omp_get_max_threads();
  CTYPE *buffer_list = (CTYPE *)malloc(
      heap, (size_t)(sizeof(CTYPE) * matrix_dim * thread_count));
  ITYPE block_size = dim / thread_count;
  ITYPE residual = dim % thread_count;

  {
    UINT thread_id = omp_get_thread_num();
    ITYPE start_index =
        block_size * thread_id + (residual > thread_id ? thread_id : residual);
    ITYPE end_index = block_size * (thread_id + 1) +
                      (residual > (thread_id + 1) ? (thread_id + 1) : residual);
    CTYPE *buffer = buffer_list + thread_id * matrix_dim;

    ITYPE state_index_y, state_index_x;
    for (state_index_x = start_index; state_index_x < end_index;
         ++state_index_x) {
      for (state_index_y = 0; state_index_y < loop_dim; ++state_index_y) {
        // create base index
        ITYPE basis_0_y = state_index_y;
        for (UINT cursor = 0; cursor < insert_index_count; cursor++) {
          UINT insert_index = sorted_insert_index_list[cursor];
          basis_0_y = insert_zero_to_basis_index(
              basis_0_y, 1ULL << insert_index, insert_index);
        }

        // flip control masks
        basis_0_y ^= control_mask;

        // compute matrix vector mul
        for (ITYPE y = 0; y < matrix_dim; ++y) {
          buffer[y] = 0;
          for (ITYPE x = 0; x < matrix_dim; ++x) {
            ITYPE dm_index_y = basis_0_y ^ matrix_mask_list[x];
            buffer[y] += matrix[y * matrix_dim + x] *
                         state[dm_index_y * dim + state_index_x];
          }
        }

        // set result
        for (ITYPE y = 0; y < matrix_dim; ++y) {
          ITYPE dm_index_y = basis_0_y ^ matrix_mask_list[y];
          state[dm_index_y * dim + state_index_x] = buffer[y];
        }
      }
    }
    for (state_index_y = start_index; state_index_y < end_index;
         ++state_index_y) {
      for (state_index_x = 0; state_index_x < loop_dim; ++state_index_x) {
        // create base index
        ITYPE basis_0_x = state_index_x;
        for (UINT cursor = 0; cursor < insert_index_count; cursor++) {
          UINT insert_index = sorted_insert_index_list[cursor];
          basis_0_x = insert_zero_to_basis_index(
              basis_0_x, 1ULL << insert_index, insert_index);
        }

        // flip control masks
        basis_0_x ^= control_mask;

        // compute matrix vector mul
        for (ITYPE y = 0; y < matrix_dim; ++y) {
          buffer[y] = 0;
          for (ITYPE x = 0; x < matrix_dim; ++x) {
            ITYPE dm_index_x = basis_0_x ^ matrix_mask_list[x];
            buffer[y] += state[state_index_y * dim + dm_index_x] *
                         adjoint_matrix[x * matrix_dim + y];
          }
        }

        // set result
        for (ITYPE y = 0; y < matrix_dim; ++y) {
          ITYPE dm_index_x = basis_0_x ^ matrix_mask_list[y];
          state[state_index_y * dim + dm_index_x] = buffer[y];
        }
      }
    }
  }
  OMPutil::get_inst().reset_qulacs_num_threads();
  free(heap, (uintptr_t)buffer_list);
#endif
  free(heap, (uintptr_t)adjoint_matrix);
  free(heap, (uintptr_t)sorted_insert_index_list);
  free(heap, (uintptr_t)matrix_mask_list);
}

void dm_X_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, PAULI_MATRIX[1], state,
                                    dim);
}
void dm_Y_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, PAULI_MATRIX[2], state,
                                    dim);
}
void dm_Z_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, PAULI_MATRIX[3], state,
                                    dim);
}
void dm_S_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, S_GATE_MATRIX, state,
                                    dim);
}
void dm_Sdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, S_DAG_GATE_MATRIX,
                                    state, dim);
}
void dm_T_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, T_GATE_MATRIX, state,
                                    dim);
}
void dm_Tdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, T_DAG_GATE_MATRIX,
                                    state, dim);
}
void dm_sqrtX_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, SQRT_X_GATE_MATRIX,
                                    state, dim);
}
void dm_sqrtXdag_gate(UINT target_qubit_index, CTYPE __global *state,
                      ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, SQRT_X_DAG_GATE_MATRIX,
                                    state, dim);
}
void dm_sqrtY_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, SQRT_Y_GATE_MATRIX,
                                    state, dim);
}
void dm_sqrtYdag_gate(UINT target_qubit_index, CTYPE __global *state,
                      ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, SQRT_Y_DAG_GATE_MATRIX,
                                    state, dim);
}
void dm_H_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, HADAMARD_MATRIX, state,
                                    dim);
}
void dm_P0_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, PROJ_0_MATRIX, state,
                                    dim);
}
void dm_P1_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  dm_single_qubit_dense_matrix_gate(target_qubit_index, PROJ_1_MATRIX, state,
                                    dim);
}
void dm_CNOT_gate(UINT control_qubit_index, UINT target_qubit_index,
                  CTYPE __global *state, ITYPE dim, HEAP heap) {
  UINT __global *control_index_list =
      (UINT __global *)malloc(heap, sizeof(UINT));
  UINT control_value_list[1];
  control_index_list[0] = control_qubit_index;
  control_value_list[0] = 1;
  dm_multi_qubit_control_single_qubit_dense_matrix_gate(
      control_index_list, control_value_list, 1, target_qubit_index,
      PAULI_MATRIX[1], state, dim, heap);
}
void dm_CZ_gate(UINT control_qubit_index, UINT target_qubit_index,
                CTYPE __global *state, ITYPE dim, HEAP heap) {
  UINT __global *control_index_list =
      (UINT __global *)malloc(heap, sizeof(UINT));
  UINT control_value_list[1];
  control_index_list[0] = control_qubit_index;
  control_value_list[0] = 1;
  dm_multi_qubit_control_single_qubit_dense_matrix_gate(
      control_index_list, control_value_list, 1, target_qubit_index,
      PAULI_MATRIX[3], state, dim, heap);
}
void dm_SWAP_gate(UINT target_qubit_index_0, UINT target_qubit_index_1,
                  CTYPE __global *state, ITYPE dim, HEAP heap) {
  MatrixXcd matrix = initialize_matrixXcd(heap, 4, 4);
  matrix[0 * 4 + 0] = double_to_complex(1, 0);
  matrix[1 * 4 + 2] = double_to_complex(1, 0);
  matrix[2 * 4 + 1] = double_to_complex(1, 0);
  matrix[3 * 4 + 3] = double_to_complex(1, 0);
  UINT __global *target_index = (UINT __global *)malloc(heap, sizeof(UINT) * 2);
  target_index[0] = target_qubit_index_0;
  target_index[1] = target_qubit_index_1;
  dm_multi_qubit_dense_matrix_gate(target_index, 2, matrix, state, dim, heap);
}
void dm_RX_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
                ITYPE dim) {
  UINT i, j;
  CTYPE rotation_gate[4];
  for (i = 0; i < 2; ++i)
    for (j = 0; j < 2; ++j)
      rotation_gate[i * 2 + j] = add(
          mul(double_to_complex(cos(angle / 2), 0), PAULI_MATRIX[0][i * 2 + j]),
          mul(mul(double_to_complex(sin(angle / 2), 0),
                  double_to_complex(0, 1.0)),
              PAULI_MATRIX[1][i * 2 + j]));
  dm_single_qubit_dense_matrix_gate_non_constant(target_qubit_index,
                                                 rotation_gate, state, dim);
}
void dm_RY_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
                ITYPE dim) {
  UINT i, j;
  CTYPE rotation_gate[4];
  for (i = 0; i < 2; ++i)
    for (j = 0; j < 2; ++j)
      rotation_gate[i * 2 + j] = add(
          mul(double_to_complex(cos(angle / 2), 0), PAULI_MATRIX[0][i * 2 + j]),
          mul(mul(double_to_complex(sin(angle / 2), 0),
                  double_to_complex(0, 1.0)),
              PAULI_MATRIX[2][i * 2 + j]));
  dm_single_qubit_dense_matrix_gate_non_constant(target_qubit_index,
                                                 rotation_gate, state, dim);
}
void dm_RZ_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
                ITYPE dim) {
  UINT i, j;
  CTYPE rotation_gate[4];
  for (i = 0; i < 2; ++i)
    for (j = 0; j < 2; ++j)
      rotation_gate[i * 2 + j] = add(
          mul(double_to_complex(cos(angle / 2), 0), PAULI_MATRIX[0][i * 2 + j]),
          mul(mul(double_to_complex(sin(angle / 2), 0),
                  double_to_complex(0, 1.0)),
              PAULI_MATRIX[3][i * 2 + j]));
  dm_single_qubit_dense_matrix_gate_non_constant(target_qubit_index,
                                                 rotation_gate, state, dim);
}

void dm_multi_qubit_Pauli_gate_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim,
    HEAP heap) {
  // TODO faster impl
  ITYPE matrix_dim = 1ULL << target_qubit_index_count;
  CTYPE __global *matrix =
      (CTYPE __global *)malloc(heap, sizeof(CTYPE) * matrix_dim * matrix_dim);
  for (ITYPE y = 0; y < matrix_dim; ++y) {
    for (ITYPE x = 0; x < matrix_dim; ++x) {
      CTYPE coef = double_to_complex(1.0, 0);
      for (UINT i = 0; i < target_qubit_index_count; ++i) {
        UINT xi = (x >> i) % 2;
        UINT yi = (y >> i) % 2;
        coef =
            mul(coef, PAULI_MATRIX[Pauli_operator_type_list[i]][yi * 2 + xi]);
      }
      matrix[y * matrix_dim + x] = coef;
    }
  }
  dm_multi_qubit_dense_matrix_gate(target_qubit_index_list,
                                   target_qubit_index_count, matrix, state, dim,
                                   heap);
  free(heap, (uintptr_t)matrix);
}
void dm_multi_qubit_Pauli_rotation_gate_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, double angle, CTYPE __global *state,
    ITYPE dim, HEAP heap) {
  // TODO faster impl
  ITYPE matrix_dim = 1ULL << target_qubit_index_count;
  CTYPE __global *matrix =
      (CTYPE __global *)malloc(heap, sizeof(CTYPE) * matrix_dim * matrix_dim);
  for (ITYPE y = 0; y < matrix_dim; ++y) {
    for (ITYPE x = 0; x < matrix_dim; ++x) {
      CTYPE coef = double_to_complex(1.0, 0);
      for (UINT i = 0; i < target_qubit_index_count; ++i) {
        UINT xi = (x >> i) % 2;
        UINT yi = (y >> i) % 2;
        coef =
            mul(coef, PAULI_MATRIX[Pauli_operator_type_list[i]][yi * 2 + xi]);
      }
      if (y == x) {
        matrix[y * matrix_dim + x] =
            add(mul(double_to_complex(cos(angle / 2), 0),
                    double_to_complex(1.0, 0)),
                mul(mul(double_to_complex(0, 1.0),
                        double_to_complex(sin(angle / 2), 0)),
                    coef));
      } else {
        matrix[y * matrix_dim + x] =
            mul(mul(double_to_complex(0, 1.0),
                    double_to_complex(sin(angle / 2), 0)),
                coef);
      }
    }
  }
  dm_multi_qubit_dense_matrix_gate(target_qubit_index_list,
                                   target_qubit_index_count, matrix, state, dim,
                                   heap);
  free(heap, (uintptr_t)matrix);
}

// state initialization
void initialize_quantum_state_parallel(CTYPE __global *state, ITYPE dim);
void initialize_quantum_state(CTYPE __global *state, ITYPE dim) {
  initialize_quantum_state_parallel(state, dim);
}

void initialize_quantum_state_parallel(CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  for (index = 0; index < dim; ++index) {
    state[index] = double_to_complex(0, 0);
  }
  state[0] = double_to_complex(1.0, 0);
}

// state randomization
unsigned long xor_shift(unsigned long *state);
double random_uniform(unsigned long *state);
double random_normal(unsigned long *state);
void initialize_Haar_random_state_with_seed_single(CTYPE __global *state,
                                                   ITYPE dim, UINT seed);
void initialize_Haar_random_state_with_seed_parallel(CTYPE __global *state,
                                                     ITYPE dim, UINT seed,
                                                     HEAP heap);

// TODO: poner una seed de verdad (time no se puede usar)
void initialize_Haar_random_state(CTYPE __global *state, ITYPE dim) {
  UINT seed = 25;
  initialize_Haar_random_state_with_seed(state, dim, seed);
}
void initialize_Haar_random_state_with_seed(CTYPE __global *state, ITYPE dim,
                                            UINT seed) {
  initialize_Haar_random_state_with_seed_single(state, dim, seed);
}

int rand(int seed) // 1 <= *seed < m
{
  int a = 16807;      // ie 7**5
  int m = 2147483647; // ie 2**31-1

  seed = ((long)(seed * a)) % m;
  return seed;
}

// single thread
void initialize_Haar_random_state_with_seed_single(CTYPE __global *state,
                                                   ITYPE dim, UINT seed) {
  int ignore_first = 40;
  double norm = 0.;
  unsigned long random_state[4];
  random_state[0] = rand(seed);
  random_state[1] = rand(seed);
  random_state[2] = rand(seed);
  random_state[3] = rand(seed);
  for (int i = 0; i < ignore_first; ++i)
    xor_shift(random_state);
  for (ITYPE index = 0; index < dim; ++index) {
    double r1, r2;
    r1 = random_normal(random_state);
    r2 = random_normal(random_state);
    state[index] = add(double_to_complex(r1, 0), mul(double_to_complex(0, 1.0),
                                                     double_to_complex(r2, 0)));
    norm += r1 * r1 + r2 * r2;
  }
  norm = sqrt(norm);
  for (ITYPE index = 0; index < dim; ++index) {
    state[index] = div(state[index], double_to_complex(norm, 0));
  }
}

unsigned long xor_shift(unsigned long *state) {
  unsigned long t;
  t = (state[0] ^ (state[0] << 11));
  state[0] = state[1];
  state[1] = state[2];
  state[2] = state[3];
  return (state[3] = (state[3] ^ (state[3] >> 19)) ^ (t ^ (t >> 8)));
}
double random_uniform(unsigned long *state) {
  return xor_shift(state) / ((double)ULONG_MAX + 1);
}
double random_normal(unsigned long *state) {
  return sqrt(-1.0 * log(1 - random_uniform(state))) *
         sin(2.0 * M_PI * random_uniform(state));
}

double expectation_value_X_Pauli_operator(UINT target_qubit_index,
                                          CTYPE __global *state, ITYPE dim);
double expectation_value_Y_Pauli_operator(UINT target_qubit_index,
                                          CTYPE __global *state, ITYPE dim);
double expectation_value_Z_Pauli_operator(UINT target_qubit_index,
                                          CTYPE __global *state, ITYPE dim);
double expectation_value_multi_qubit_Pauli_operator_XZ_mask(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, CTYPE __global *state, ITYPE dim);
double expectation_value_multi_qubit_Pauli_operator_Z_mask(
    ITYPE phase_flip_mask, CTYPE __global *state, ITYPE dim);

// calculate expectation value of X on target qubit
double expectation_value_X_Pauli_operator(UINT target_qubit_index,
                                          CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = 1ULL << target_qubit_index;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_0 =
        insert_zero_to_basis_index(state_index, mask, target_qubit_index);
    ITYPE basis_1 = basis_0 ^ mask;

    sum += _creal(mul(conj(state[basis_0]), state[basis_1])) * 2;
  }

  return sum;
}

// calculate expectation value of Y on target qubit
double expectation_value_Y_Pauli_operator(UINT target_qubit_index,
                                          CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = 1ULL << target_qubit_index;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_0 =
        insert_zero_to_basis_index(state_index, mask, target_qubit_index);
    ITYPE basis_1 = basis_0 ^ mask;
    sum += _cimag(mul(conj(state[basis_0]), state[basis_1])) * 2;
  }

  return sum;
}

// calculate expectation value of Z on target qubit
double expectation_value_Z_Pauli_operator(UINT target_qubit_index,
                                          CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    int sign = 1 - 2 * ((state_index >> target_qubit_index) % 2);
    sum += _creal(mul(conj(state[state_index]), state[state_index])) * sign;
  }

  return sum;
}

// calculate expectation value for single-qubit pauli operator
double expectation_value_single_qubit_Pauli_operator(UINT target_qubit_index,
                                                     UINT Pauli_operator_type,
                                                     CTYPE __global *state,
                                                     ITYPE dim) {
  if (Pauli_operator_type == 0) {
    return state_norm_squared(state, dim);
  } else if (Pauli_operator_type == 1) {
    return expectation_value_X_Pauli_operator(target_qubit_index, state, dim);
  } else if (Pauli_operator_type == 2) {
    return expectation_value_Y_Pauli_operator(target_qubit_index, state, dim);
  } else if (Pauli_operator_type == 3) {
    return expectation_value_Z_Pauli_operator(target_qubit_index, state, dim);
  } else {
    printf("ERROR: invalid expectation value of pauli operator is called");
    return 0;
  }
}

// calculate expectation value of multi-qubit Pauli operator on qubits.
// bit-flip mask : the n-bit binary string of which the i-th element is 1 iff
// the i-th pauli operator is X or Y phase-flip mask : the n-bit binary string
// of which the i-th element is 1 iff the i-th pauli operator is Y or Z We
// assume bit-flip mask is nonzero, namely, there is at least one X or Y
// operator. the pivot qubit is any qubit index which has X or Y To generate
// bit-flip mask and phase-flip mask, see get_masks_*_list at utility.h
double expectation_value_multi_qubit_Pauli_operator_XZ_mask(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE pivot_mask = 1ULL << pivot_qubit_index;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_0 =
        insert_zero_to_basis_index(state_index, pivot_mask, pivot_qubit_index);
    ITYPE basis_1 = basis_0 ^ bit_flip_mask;
    UINT sign_0 = count_population(basis_0 & phase_flip_mask) % 2;

    sum +=
        _creal(mul(mul(state[basis_0], conj(state[basis_1])),
                   PHASE_90ROT[(global_phase_90rot_count + sign_0 * 2) % 4])) *
        2.0;
  }

  return sum;
}

double expectation_value_multi_qubit_Pauli_operator_Z_mask(
    ITYPE phase_flip_mask, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    int bit_parity = count_population(state_index & phase_flip_mask) % 2;
    int sign = 1 - 2 * bit_parity;
    sum += pow(_cabs(state[state_index]), 2) * sign;
  }

  return sum;
}

double expectation_value_multi_qubit_Pauli_operator_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim) {
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_partial_list(
      target_qubit_index_list, Pauli_operator_type_list,
      target_qubit_index_count, &bit_flip_mask, &phase_flip_mask,
      &global_phase_90rot_count, &pivot_qubit_index);
  double result;
  if (bit_flip_mask == 0) {
    result = expectation_value_multi_qubit_Pauli_operator_Z_mask(
        phase_flip_mask, state, dim);
  } else {
    result = expectation_value_multi_qubit_Pauli_operator_XZ_mask(
        bit_flip_mask, phase_flip_mask, global_phase_90rot_count,
        pivot_qubit_index, state, dim);
  }
  return result;
}

double expectation_value_multi_qubit_Pauli_operator_whole_list(
    UINT *Pauli_operator_type_list, UINT qubit_count, CTYPE __global *state,
    ITYPE dim) {
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_whole_list(Pauli_operator_type_list, qubit_count,
                             &bit_flip_mask, &phase_flip_mask,
                             &global_phase_90rot_count, &pivot_qubit_index);
  double result;
  if (bit_flip_mask == 0) {
    result = expectation_value_multi_qubit_Pauli_operator_Z_mask(
        phase_flip_mask, state, dim);
  } else {
    result = expectation_value_multi_qubit_Pauli_operator_XZ_mask(
        bit_flip_mask, phase_flip_mask, global_phase_90rot_count,
        pivot_qubit_index, state, dim);
  }
  return result;
}

/****
 * Single thread version of expectation value
 **/
// calculate expectation value of multi-qubit Pauli operator on qubits.
// bit-flip mask : the n-bit binary string of which the i-th element is 1 iff
// the i-th pauli operator is X or Y phase-flip mask : the n-bit binary string
// of which the i-th element is 1 iff the i-th pauli operator is Y or Z We
// assume bit-flip mask is nonzero, namely, there is at least one X or Y
// operator. the pivot qubit is any qubit index which has X or Y To generate
// bit-flip mask and phase-flip mask, see get_masks_*_list at utility.h
double expectation_value_multi_qubit_Pauli_operator_XZ_mask_single_thread(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE pivot_mask = 1ULL << pivot_qubit_index;
  ITYPE state_index;
  double sum = 0.;
  UINT sign_0;
  ITYPE basis_0, basis_1;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    basis_0 =
        insert_zero_to_basis_index(state_index, pivot_mask, pivot_qubit_index);
    basis_1 = basis_0 ^ bit_flip_mask;
    sign_0 = count_population(basis_0 & phase_flip_mask) % 2;
    sum +=
        _creal(mul(mul(state[basis_0], conj(state[basis_1])),
                   PHASE_90ROT[(global_phase_90rot_count + sign_0 * 2) % 4])) *
        2.0;
  }
  return sum;
}

double expectation_value_multi_qubit_Pauli_operator_Z_mask_single_thread(
    ITYPE phase_flip_mask, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim;
  ITYPE state_index;
  double sum = 0.;
  int bit_parity, sign;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    bit_parity = count_population(state_index & phase_flip_mask) % 2;
    sign = 1 - 2 * bit_parity;
    sum += pow(_cabs(state[state_index]), 2) * sign;
  }
  return sum;
}

double expectation_value_multi_qubit_Pauli_operator_partial_list_single_thread(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim) {
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_partial_list(
      target_qubit_index_list, Pauli_operator_type_list,
      target_qubit_index_count, &bit_flip_mask, &phase_flip_mask,
      &global_phase_90rot_count, &pivot_qubit_index);
  double result;
  if (bit_flip_mask == 0) {
    result = expectation_value_multi_qubit_Pauli_operator_Z_mask_single_thread(
        phase_flip_mask, state, dim);
  } else {
    result = expectation_value_multi_qubit_Pauli_operator_XZ_mask_single_thread(
        bit_flip_mask, phase_flip_mask, global_phase_90rot_count,
        pivot_qubit_index, state, dim);
  }
  return result;
}

// calculate probability with which we obtain 0 at target qubit
double M0_prob(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = 1ULL << target_qubit_index;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_0 =
        insert_zero_to_basis_index(state_index, mask, target_qubit_index);
    sum += pow(_cabs(state[basis_0]), 2);
  }

  return sum;
}

// calculate probability with which we obtain 1 at target qubit
double M1_prob(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = 1ULL << target_qubit_index;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_1 =
        insert_zero_to_basis_index(state_index, mask, target_qubit_index) ^
        mask;
    sum += pow(_cabs(state[basis_1]), 2);
  }

  return sum;
}

// calculate merginal probability with which we obtain the set of values
// measured_value_list at sorted_target_qubit_index_list warning:
// sorted_target_qubit_index_list must be sorted.
double marginal_prob( UINT* sorted_target_qubit_index_list,
     UINT __global* measured_value_list, UINT target_qubit_index_count,
     CTYPE __global* state, ITYPE dim)  {
  ITYPE loop_dim = dim >> target_qubit_index_count;
  ITYPE state_index;
  double sum = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis = state_index;
    for (UINT cursor = 0; cursor < target_qubit_index_count; cursor++) {
      UINT insert_index = sorted_target_qubit_index_list[cursor];
      ITYPE mask = 1ULL << insert_index;
      basis = insert_zero_to_basis_index(basis, mask, insert_index);
      basis ^= mask * measured_value_list[cursor];
    }
    sum += pow(_cabs(state[basis]), 2);
  }

  return sum;
}

// calculate entropy of probability distribution of Z-basis measurements
double measurement_distribution_entropy(CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  double ent = 0;
  double eps = 1e-15;

  for (index = 0; index < dim; ++index) {
    double prob = pow(_cabs(state[index]), 2);
    prob = (prob > eps) ? prob : eps;
    ent += -1.0 * prob * log(prob);
  }

  return ent;
}

CTYPE
transition_amplitude_multi_qubit_Pauli_operator_XZ_mask(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, CTYPE __global *state_bra,
    CTYPE __global *state_ket, ITYPE dim);
CTYPE
transition_amplitude_multi_qubit_Pauli_operator_Z_mask(
    ITYPE phase_flip_mask, CTYPE __global *state_bra, CTYPE __global *state_ket,
    ITYPE dim);

CTYPE
transition_amplitude_multi_qubit_Pauli_operator_XZ_mask(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, CTYPE __global *state_bra,
    CTYPE __global *state_ket, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE pivot_mask = 1ULL << pivot_qubit_index;
  ITYPE state_index;

  double sum_real = 0.;
  double sum_imag = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_0 =
        insert_zero_to_basis_index(state_index, pivot_mask, pivot_qubit_index);
    ITYPE basis_1 = basis_0 ^ bit_flip_mask;
    UINT sign_0 = count_population(basis_0 & phase_flip_mask) % 2;
    UINT sign_1 = count_population(basis_1 & phase_flip_mask) % 2;
    CTYPE val1 = mul(mul(state_ket[basis_0], conj(state_bra[basis_1])),
                     PHASE_90ROT[(global_phase_90rot_count + sign_0 * 2) % 4]);
    CTYPE val2 = mul(mul(state_ket[basis_1], conj(state_bra[basis_0])),
                     PHASE_90ROT[(global_phase_90rot_count + sign_1 * 2) % 4]);
    sum_real += _creal(val1);
    sum_imag += _cimag(val1);
    sum_real += _creal(val2);
    sum_imag += _cimag(val2);
  }
  CTYPE sum = double_to_complex(sum_real, sum_imag);

  return sum;
}

CTYPE
transition_amplitude_multi_qubit_Pauli_operator_Z_mask(
    ITYPE phase_flip_mask, CTYPE __global *state_bra, CTYPE __global *state_ket,
    ITYPE dim) {
  ITYPE loop_dim = dim;
  ITYPE state_index;

  double sum_real = 0.;
  double sum_imag = 0.;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    int bit_parity = count_population(state_index & phase_flip_mask) % 2;
    double sign = 1 - 2 * bit_parity;
    CTYPE val = mul(mul(double_to_complex(sign, 0), state_ket[state_index]),
                    conj(state_bra[state_index]));
    sum_real += _creal(val);
    sum_imag += _cimag(val);
  }
  CTYPE sum = double_to_complex(sum_real, sum_imag);

  return sum;
}

CTYPE
transition_amplitude_multi_qubit_Pauli_operator_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state_bra,
    CTYPE __global *state_ket, ITYPE dim) {
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_partial_list(
      target_qubit_index_list, Pauli_operator_type_list,
      target_qubit_index_count, &bit_flip_mask, &phase_flip_mask,
      &global_phase_90rot_count, &pivot_qubit_index);
  CTYPE result;
  if (bit_flip_mask == 0) {
    result = transition_amplitude_multi_qubit_Pauli_operator_Z_mask(
        phase_flip_mask, state_bra, state_ket, dim);
  } else {
    result = transition_amplitude_multi_qubit_Pauli_operator_XZ_mask(
        bit_flip_mask, phase_flip_mask, global_phase_90rot_count,
        pivot_qubit_index, state_bra, state_ket, dim);
  }
  return result;
}

CTYPE
transition_amplitude_multi_qubit_Pauli_operator_whole_list(
    UINT *Pauli_operator_type_list, UINT qubit_count, CTYPE __global *state_bra,
    CTYPE __global *state_ket, ITYPE dim) {
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_whole_list(Pauli_operator_type_list, qubit_count,
                             &bit_flip_mask, &phase_flip_mask,
                             &global_phase_90rot_count, &pivot_qubit_index);
  CTYPE result;
  if (bit_flip_mask == 0) {
    result = transition_amplitude_multi_qubit_Pauli_operator_Z_mask(
        phase_flip_mask, state_bra, state_ket, dim);
  } else {
    result = transition_amplitude_multi_qubit_Pauli_operator_XZ_mask(
        bit_flip_mask, phase_flip_mask, global_phase_90rot_count,
        pivot_qubit_index, state_bra, state_ket, dim);
  }
  return result;
}

void multi_qubit_control_multi_qubit_dense_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT __global *target_qubit_index_list,
    UINT target_qubit_index_count, CTYPE __global *matrix,
    CTYPE __global *state, ITYPE dim, HEAP heap) {
  // matrix dim, mask, buffer
  ITYPE matrix_dim = 1ULL << target_qubit_index_count;
  ITYPE __global *matrix_mask_list = create_matrix_mask_list(
      target_qubit_index_list, target_qubit_index_count, heap);
  CTYPE __global *buffer =
      (CTYPE __global *)malloc(heap, (size_t)(sizeof(CTYPE) * matrix_dim));

  // insert index
  UINT insert_index_count =
      target_qubit_index_count + control_qubit_index_count;
  UINT __global *sorted_insert_index_list = create_sorted_ui_list_list(
      target_qubit_index_list, target_qubit_index_count,
      control_qubit_index_list, control_qubit_index_count, heap);

  // control mask
  ITYPE control_mask = create_control_mask(
      control_qubit_index_list, control_value_list, control_qubit_index_count);

  // loop varaibles
  ITYPE loop_dim =
      dim >> (target_qubit_index_count + control_qubit_index_count);
  ITYPE state_index;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create base index
    ITYPE basis_0 = state_index;
    for (UINT cursor = 0; cursor < insert_index_count; cursor++) {
      UINT insert_index = sorted_insert_index_list[cursor];
      basis_0 = insert_zero_to_basis_index(basis_0, 1ULL << insert_index,
                                           insert_index);
    }

    // flip control masks
    basis_0 ^= control_mask;

    // compute matrix mul
    for (ITYPE y = 0; y < matrix_dim; ++y) {
      buffer[y] = double_to_complex(0, 0);
      for (ITYPE x = 0; x < matrix_dim; ++x) {
        buffer[y] = add(buffer[y], mul(matrix[y * matrix_dim + x],
                                       state[basis_0 ^ matrix_mask_list[x]]));
      }
    }

    // set result
    for (ITYPE y = 0; y < matrix_dim; ++y) {
      state[basis_0 ^ matrix_mask_list[y]] = buffer[y];
    }
  }
  free(heap, (uintptr_t)sorted_insert_index_list);
  free(heap, (uintptr_t)buffer);
  free(heap, (uintptr_t)matrix_mask_list);
}

void create_shift_mask_list_from_list_and_value_buf(UINT __global *array,
                                                    UINT count, UINT target,
                                                    UINT __global *dst_array,
                                                    ITYPE *dst_mask,
                                                    HEAP heap) {
  UINT size = count + 1;
  memcpy(dst_array, array, sizeof(UINT) * count);
  dst_array[count] = target;
  sort_ui(dst_array, size, heap);
  for (UINT i = 0; i < size; ++i) {
    dst_mask[i] = (1UL << dst_array[i]) - 1;
  }
}

void multi_qubit_control_single_qubit_dense_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT target_qubit_index, CTYPE matrix[4],
    CTYPE __global *state, ITYPE dim, HEAP heap) {
  if (control_qubit_index_count == 1) {
    single_qubit_control_single_qubit_dense_matrix_gate(
        control_qubit_index_list[0], control_value_list[0], target_qubit_index,
        matrix, state, dim);
    return;
  }

  multi_qubit_control_single_qubit_dense_matrix_gate_unroll(
      control_qubit_index_list, control_value_list, control_qubit_index_count,
      target_qubit_index, matrix, state, dim, heap);
}

void multi_qubit_control_single_qubit_dense_matrix_gate_unroll(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT target_qubit_index, CTYPE matrix[4],
    CTYPE __global *state, ITYPE dim, HEAP heap) {
  UINT __global *sort_array = (UINT __global *)malloc(heap, sizeof(UINT) * 64);
  ITYPE mask_array[64];
  create_shift_mask_list_from_list_and_value_buf(
      control_qubit_index_list, control_qubit_index_count, target_qubit_index,
      sort_array, mask_array, heap);
  ITYPE target_mask = 1ULL << target_qubit_index;
  ITYPE control_mask = create_control_mask(
      control_qubit_index_list, control_value_list, control_qubit_index_count);

  UINT insert_index_list_count = control_qubit_index_count + 1;
  ITYPE loop_dim = dim >> insert_index_list_count;

  if (target_qubit_index == 0) {
    ITYPE state_index;
    for (state_index = 0; state_index < loop_dim; ++state_index) {
      // create base index
      ITYPE basis_0 = state_index;
      for (UINT cursor = 0; cursor < insert_index_list_count; ++cursor) {
        basis_0 = (basis_0 & mask_array[cursor]) +
                  ((basis_0 & (~mask_array[cursor])) << 1);
      }
      basis_0 += control_mask;

      // fetch values
      CTYPE cval0 = state[basis_0];
      CTYPE cval1 = state[basis_0 + 1];
      // set values
      state[basis_0] = add(mul(matrix[0], cval0), mul(matrix[1], cval1));
      state[basis_0 + 1] = add(mul(matrix[2], cval0), add(matrix[3], cval1));
    }
  } else if (sort_array[0] == 0) {
    ITYPE state_index;
    for (state_index = 0; state_index < loop_dim; ++state_index) {
      // create base index
      ITYPE basis_0 = state_index;
      for (UINT cursor = 0; cursor < insert_index_list_count; ++cursor) {
        basis_0 = (basis_0 & mask_array[cursor]) +
                  ((basis_0 & (~mask_array[cursor])) << 1);
      }
      basis_0 += control_mask;
      ITYPE basis_1 = basis_0 + target_mask;

      // fetch values
      CTYPE cval0 = state[basis_0];
      CTYPE cval1 = state[basis_1];
      // set values
      state[basis_0] = add(mul(matrix[0], cval0), mul(matrix[1], cval1));
      state[basis_1] = add(mul(matrix[2], cval0), mul(matrix[3], cval1));
    }
  } else {
    ITYPE state_index;
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      // create base index
      ITYPE basis_0 = state_index;
      for (UINT cursor = 0; cursor < insert_index_list_count; ++cursor) {
        basis_0 = (basis_0 & mask_array[cursor]) +
                  ((basis_0 & (~mask_array[cursor])) << 1);
      }
      basis_0 += control_mask;
      ITYPE basis_1 = basis_0 + target_mask;

      // fetch values
      CTYPE cval0 = state[basis_0];
      CTYPE cval1 = state[basis_1];
      CTYPE cval2 = state[basis_0 + 1];
      CTYPE cval3 = state[basis_1 + 1];
      // set values
      state[basis_0] = add(mul(matrix[0], cval0), mul(matrix[1], cval1));
      state[basis_1] = add(mul(matrix[2], cval0), mul(matrix[3], cval1));
      state[basis_0 + 1] = add(mul(matrix[0], cval2), mul(matrix[1], cval3));
      state[basis_1 + 1] = add(mul(matrix[2], cval2), mul(matrix[3], cval3));
    }
  }
}

void single_qubit_control_multi_qubit_dense_matrix_gate(
    UINT control_qubit_index, UINT control_value,
    UINT __global *target_qubit_index_list, UINT target_qubit_index_count,
    CTYPE matrix[16], CTYPE __global *state, ITYPE dim, HEAP heap) {
  // matrix dim, mask, buffer
  ITYPE matrix_dim = 1ULL << target_qubit_index_count;
  ITYPE __global* matrix_mask_list = create_matrix_mask_list(
        target_qubit_index_list, target_qubit_index_count, heap);
  CTYPE __global *buffer =
      (CTYPE __global *)malloc(heap, (size_t)(sizeof(CTYPE) * matrix_dim));

  // insert list
  UINT insert_index_count = target_qubit_index_count + 1;
  UINT __global *sorted_insert_index_list = create_sorted_ui_list_value(
      target_qubit_index_list, target_qubit_index_count, control_qubit_index,
      heap);
  
  /*printf("Lista ordenada: ");
  for(int i=0; i < target_qubit_index_count + 1; i++){
    printf("%d ", sorted_insert_index_list[i]);
  }
  printf("\n\n");*/

  // control mask
  ITYPE control_mask = (1ULL << control_qubit_index) * control_value;

  // loop varaibles
  ITYPE loop_dim = dim >> insert_index_count;
  ITYPE state_index;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create base index
    ITYPE basis_0 = state_index;
    for (UINT cursor = 0; cursor < insert_index_count; cursor++) {
      UINT insert_index = sorted_insert_index_list[cursor];
      basis_0 = insert_zero_to_basis_index(basis_0, 1ULL << insert_index,
                                           insert_index);
    }

    // flip control
    basis_0 ^= control_mask;

    // compute matrix mul
    for (ITYPE y = 0; y < matrix_dim; ++y) {
      buffer[y] = double_to_complex(0, 0);
      for (ITYPE x = 0; x < matrix_dim; ++x) {
        buffer[y] = add(buffer[y], mul(matrix[y * matrix_dim + x],
                                       state[basis_0 ^ matrix_mask_list[x]]));
      }
      //printf("\nBuffer:\n");
      //printComplex(buffer[y]);
    }
    //printf("\n\n");

    // set result
    for (ITYPE y = 0; y < matrix_dim; ++y) {
      state[basis_0 ^ matrix_mask_list[y]] = buffer[y];
    }
  }
  free(heap, (uintptr_t)sorted_insert_index_list);
  free(heap, (uintptr_t)buffer);
  free(heap, (uintptr_t)matrix_mask_list);
}

void single_qubit_control_single_qubit_dense_matrix_gate(
    UINT control_qubit_index, UINT control_value, UINT target_qubit_index,
    CTYPE matrix[4], CTYPE __global *state, ITYPE dim) {

  single_qubit_control_single_qubit_dense_matrix_gate_unroll(
      control_qubit_index, control_value, target_qubit_index, matrix, state,
      dim);
}

void single_qubit_control_single_qubit_dense_matrix_gate_unroll(
    UINT control_qubit_index, UINT control_value, UINT target_qubit_index,
    CTYPE matrix[4], CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 4;

  ITYPE target_mask = 1ULL << target_qubit_index;
  ITYPE control_mask = 1ULL << control_qubit_index;

  UINT min_qubit_index = get_min_ui(control_qubit_index, target_qubit_index);
  UINT max_qubit_index = get_max_ui(control_qubit_index, target_qubit_index);
  ITYPE min_qubit_mask = 1ULL << min_qubit_index;
  ITYPE max_qubit_mask = 1ULL << (max_qubit_index - 1);
  ITYPE low_mask = min_qubit_mask - 1;
  ITYPE mid_mask = (max_qubit_mask - 1) ^ low_mask;
  ITYPE high_mask = ~(max_qubit_mask - 1);

  ITYPE state_index;
  if (target_qubit_index == 0) {
    for (state_index = 0; state_index < loop_dim; ++state_index) {
      ITYPE basis_index =
          (state_index & low_mask) + ((state_index & mid_mask) << 1) +
          ((state_index & high_mask) << 2) + control_mask * control_value;

      // fetch values
      CTYPE cval0 = state[basis_index];
      CTYPE cval1 = state[basis_index + 1];

      // set values
      state[basis_index] = add(mul(matrix[0], cval0), mul(matrix[1], cval1));
      state[basis_index + 1] =
          add(mul(matrix[2], cval0), mul(matrix[3], cval1));
    }
  } else if (control_qubit_index == 0) {
    for (state_index = 0; state_index < loop_dim; ++state_index) {
      ITYPE basis_index_0 =
          (state_index & low_mask) + ((state_index & mid_mask) << 1) +
          ((state_index & high_mask) << 2) + control_mask * control_value;
      ITYPE basis_index_1 = basis_index_0 + target_mask;

      // fetch values
      CTYPE cval0 = state[basis_index_0];
      CTYPE cval1 = state[basis_index_1];

      // set values
      state[basis_index_0] = add(mul(matrix[0], cval0), mul(matrix[1], cval1));
      state[basis_index_1] = add(mul(matrix[2], cval0), mul(matrix[3], cval1));
    }
  } else {
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis_index_0 =
          (state_index & low_mask) + ((state_index & mid_mask) << 1) +
          ((state_index & high_mask) << 2) + control_mask * control_value;
      ITYPE basis_index_1 = basis_index_0 + target_mask;

      // fetch values
      CTYPE cval0 = state[basis_index_0];
      CTYPE cval1 = state[basis_index_1];
      CTYPE cval2 = state[basis_index_0 + 1];
      CTYPE cval3 = state[basis_index_1 + 1];

      // set values
      state[basis_index_0] = add(mul(matrix[0], cval0), mul(matrix[1], cval1));
      state[basis_index_1] = add(mul(matrix[2], cval0), mul(matrix[3], cval1));
      state[basis_index_0 + 1] =
          add(mul(matrix[0], cval2), mul(matrix[1], cval3));
      state[basis_index_1 + 1] =
          add(mul(matrix[2], cval2), mul(matrix[3], cval3));
    }
  }
}

void double_qubit_dense_matrix_gate_c(UINT target_qubit_index1,
                                      UINT target_qubit_index2,
                                      __constant CTYPE matrix[16],
                                      CTYPE __global *state, ITYPE dim) {
  double_qubit_dense_matrix_gate_nosimd(
      target_qubit_index1, target_qubit_index2, matrix, state, dim);
}

void double_qubit_dense_matrix_gate_nosimd(UINT target_qubit_index1,
                                           UINT target_qubit_index2,
                                           __constant CTYPE matrix[16],
                                           CTYPE __global *state, ITYPE dim) {
  UINT min_qubit_index = get_min_ui(target_qubit_index1, target_qubit_index2);
  UINT max_qubit_index = get_max_ui(target_qubit_index1, target_qubit_index2);
  ITYPE min_qubit_mask = 1ULL << min_qubit_index;
  ITYPE max_qubit_mask = 1ULL << (max_qubit_index - 1);
  ITYPE low_mask = min_qubit_mask - 1;
  ITYPE mid_mask = (max_qubit_mask - 1) ^ low_mask;
  ITYPE high_mask = ~(max_qubit_mask - 1);

  ITYPE target_mask1 = 1ULL << target_qubit_index1;
  ITYPE target_mask2 = 1ULL << target_qubit_index2;

  // loop variables
  ITYPE loop_dim = dim / 4;
  ITYPE state_index;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create index
    ITYPE basis_0 = (state_index & low_mask) + ((state_index & mid_mask) << 1) +
                    ((state_index & high_mask) << 2);

    // gather index
    ITYPE basis_1 = basis_0 + target_mask1;
    ITYPE basis_2 = basis_0 + target_mask2;
    ITYPE basis_3 = basis_1 + target_mask2;

    // fetch values
    CTYPE cval_0 = state[basis_0];
    CTYPE cval_1 = state[basis_1];
    CTYPE cval_2 = state[basis_2];
    CTYPE cval_3 = state[basis_3];

    // set values
    state[basis_0] =
        add(add(add(mul(matrix[0], cval_0), mul(matrix[1], cval_1)),
                mul(matrix[2], cval_2)),
            mul(matrix[3], cval_3));
    state[basis_1] =
        add(add(add(mul(matrix[4], cval_0), mul(matrix[5], cval_1)),
                mul(matrix[6], cval_2)),
            mul(matrix[7], cval_3));
    state[basis_2] =
        add(add(add(mul(matrix[8], cval_0), mul(matrix[9], cval_1)),
                mul(matrix[10], cval_2)),
            mul(matrix[11], cval_3));
    state[basis_3] =
        add(add(add(mul(matrix[12], cval_0), mul(matrix[13], cval_1)),
                mul(matrix[14], cval_2)),
            mul(matrix[15], cval_3));
  }
}

#define EIGEN_DONT_PARALLELIZE

void double_qubit_dense_matrix_gate(UINT target_qubit_index1,
                                    UINT target_qubit_index2,
                                    __constant CTYPE matrix[16],
                                    CTYPE __global *state, ITYPE dim) {
  double_qubit_dense_matrix_gate_c(target_qubit_index1, target_qubit_index2,
                                   matrix, state, dim);
}

void double_qubit_dense_matrix_gate4cd(UINT target_qubit_index1,
                                       UINT target_qubit_index2,
                                       Matrix4cd eigen_matrix,
                                       CTYPE __global *state, ITYPE dim,
                                       HEAP heap) {
  double_qubit_dense_matrix_gate_eigen(target_qubit_index1, target_qubit_index2,
                                       eigen_matrix, state, dim, heap);
}

void double_qubit_dense_matrix_gate_eigen(UINT target_qubit_index1,
                                          UINT target_qubit_index2,
                                          Matrix4cd eigen_matrix,
                                          CTYPE __global *state, ITYPE dim,
                                          HEAP heap) {
  // target mask

  UINT min_qubit_index = get_min_ui(target_qubit_index1, target_qubit_index2);
  UINT max_qubit_index = get_max_ui(target_qubit_index1, target_qubit_index2);
  ITYPE min_qubit_mask = 1ULL << min_qubit_index;
  ITYPE max_qubit_mask = 1ULL << (max_qubit_index - 1);
  ITYPE low_mask = min_qubit_mask - 1;
  ITYPE mid_mask = (max_qubit_mask - 1) ^ low_mask;
  ITYPE high_mask = ~(max_qubit_mask - 1);

  ITYPE target_mask1 = 1ULL << target_qubit_index1;
  ITYPE target_mask2 = 1ULL << target_qubit_index2;
  CTYPE __global *eigen_state = state;

  // loop variables
  ITYPE loop_dim = dim / 4.0;
  ITYPE state_index;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create index
    ITYPE basis_0 = (state_index & low_mask) + ((state_index & mid_mask) << 1) +
                    ((state_index & high_mask) << 2);

    // gather index
    ITYPE basis_1 = basis_0 + target_mask1;
    ITYPE basis_2 = basis_0 + target_mask2;
    ITYPE basis_3 = basis_1 + target_mask2;

    // fetch values
    Vector4cd vec, res;
    initialise_vector4cd(res);
    vec[0] = state[basis_0];
    vec[1] = state[basis_1];
    vec[2] = state[basis_2];
    vec[3] = state[basis_3];

    mul4cd(eigen_matrix, vec, res);
    eigen_state[basis_0] = res[0];
    eigen_state[basis_1] = res[1];
    eigen_state[basis_2] = res[2];
    eigen_state[basis_3] = res[3];
  }
}

void create_shift_mask_list_from_list_buf(UINT __global *array, UINT count,
                                          UINT __global *dst_array,
                                          ITYPE *dst_mask, HEAP heap);

void multi_qubit_dense_matrix_gate(UINT __global *target_qubit_index_list,
                                   UINT target_qubit_index_count,
                                   __constant CTYPE *matrix,
                                   CTYPE __global *state, ITYPE dim,
                                   HEAP heap) {
  if (target_qubit_index_count == 1) {
    single_qubit_dense_matrix_gate(target_qubit_index_list[0], matrix, state,
                                   dim);
  } else if (target_qubit_index_count == 2) {
    double_qubit_dense_matrix_gate_c(target_qubit_index_list[0],
                                     target_qubit_index_list[1], matrix, state,
                                     dim);
  } else {
    multi_qubit_dense_matrix_gate_parallel(target_qubit_index_list,
                                           target_qubit_index_count, matrix,
                                           state, dim, heap);
  }
}

void create_shift_mask_list_from_list_buf(UINT __global *array, UINT count,
                                          UINT __global *dst_array,
                                          ITYPE *dst_mask, HEAP heap) {
  memcpy(dst_array, array, sizeof(UINT) * count);
  sort_ui(dst_array, count, heap);
  for (UINT i = 0; i < count; ++i) {
    dst_mask[i] = (1UL << dst_array[i]) - 1;
  }
}

void multi_qubit_dense_matrix_gate_parallel(
    UINT __global *target_qubit_index_list, UINT target_qubit_index_count,
    __constant CTYPE *matrix, CTYPE __global *state, ITYPE dim, HEAP heap) {
  UINT __global *sort_array = (UINT __global *)malloc(heap, sizeof(UINT) * 64);
  ITYPE mask_array[64];
  create_shift_mask_list_from_list_buf(target_qubit_index_list,
                                       target_qubit_index_count, sort_array,
                                       mask_array, heap);

  // matrix dim, mask, buffer
  ITYPE matrix_dim = 1ULL << target_qubit_index_count;
  ITYPE __global *matrix_mask_list = create_matrix_mask_list(
      target_qubit_index_list, target_qubit_index_count, heap);
  // loop variables
  ITYPE loop_dim = dim >> target_qubit_index_count;

  UINT thread_count = 1;
  CTYPE __global *buffer_list = (CTYPE __global *)malloc(
      heap, (size_t)(sizeof(CTYPE) * matrix_dim * thread_count));

  ITYPE block_size = loop_dim / thread_count;
  ITYPE residual = loop_dim % thread_count;

  {
    UINT thread_id = 0;
    ITYPE start_index =
        block_size * thread_id + (residual > thread_id ? thread_id : residual);
    ITYPE end_index = block_size * (thread_id + 1) +
                      (residual > (thread_id + 1) ? (thread_id + 1) : residual);
    CTYPE __global *buffer = buffer_list + thread_id * matrix_dim;

    ITYPE state_index;
    for (state_index = start_index; state_index < end_index; ++state_index) {
      // create base index
      ITYPE basis_0 = state_index;
      for (UINT cursor = 0; cursor < target_qubit_index_count; ++cursor) {
        basis_0 = (basis_0 & mask_array[cursor]) +
                  ((basis_0 & (~mask_array[cursor])) << 1);
      }

      // compute matrix-vector multiply
      for (ITYPE y = 0; y < matrix_dim; ++y) {
        buffer[y] = double_to_complex(0, 0);
        for (ITYPE x = 0; x < matrix_dim; ++x) {
          buffer[y] = add(buffer[y], mul(matrix[y * matrix_dim + x],
                                         state[basis_0 ^ matrix_mask_list[x]]));
        }
      }

      // set result
      for (ITYPE y = 0; y < matrix_dim; ++y) {
        state[basis_0 ^ matrix_mask_list[y]] = buffer[y];
      }
    }
  }
  free(heap, (uintptr_t)buffer_list);
  free(heap, (uintptr_t)matrix_mask_list);
}

void single_qubit_dense_matrix_gate(UINT target_qubit_index,
                                    __constant CTYPE matrix[4],
                                    CTYPE __global *state, ITYPE dim) {
  single_qubit_dense_matrix_gate_parallel(target_qubit_index, matrix, state,
                                          dim);
}

void single_qubit_dense_matrix_gate_non_constant(UINT target_qubit_index,
                                                 CTYPE matrix[4],
                                                 CTYPE __global *state,
                                                 ITYPE dim) {
  single_qubit_dense_matrix_gate_parallel_no_constant(target_qubit_index,
                                                      matrix, state, dim);
}

void single_qubit_dense_matrix_gate_parallel_no_constant(
    UINT target_qubit_index, CTYPE matrix[4], CTYPE __global *state,
    ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = (1ULL << target_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;

  ITYPE state_index = 0;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_0 = (state_index & mask_low) + ((state_index & mask_high) << 1);
    ITYPE basis_1 = basis_0 + mask;

    // fetch values
    CTYPE cval_0 = state[basis_0];
    CTYPE cval_1 = state[basis_1];

    // set values
    state[basis_0] = add(mul(matrix[0], cval_0), mul(matrix[1], cval_1));
    state[basis_1] = add(mul(matrix[2], cval_0), mul(matrix[3], cval_1));
  }
}

void single_qubit_dense_matrix_gate_parallel(UINT target_qubit_index,
                                             __constant CTYPE matrix[4],
                                             CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = (1ULL << target_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;

  ITYPE state_index = 0;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE basis_0 = (state_index & mask_low) + ((state_index & mask_high) << 1);
    ITYPE basis_1 = basis_0 + mask;

    // fetch values
    CTYPE cval_0 = state[basis_0];
    CTYPE cval_1 = state[basis_1];

    // set values
    state[basis_0] = add(mul(matrix[0], cval_0), mul(matrix[1], cval_1));
    state[basis_1] = add(mul(matrix[2], cval_0), mul(matrix[3], cval_1));
  }
}

void single_qubit_dense_matrix_gate_parallel_unroll(UINT target_qubit_index,
                                                    CTYPE matrix[4],
                                                    CTYPE __global *state,
                                                    ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = (1ULL << target_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;

  if (target_qubit_index == 0) {
    ITYPE basis = 0;
    for (basis = 0; basis < dim; basis += 2) {
      CTYPE val0a = state[basis];
      CTYPE val1a = state[basis + 1];
      CTYPE res0a = add(mul(val0a, matrix[0]), mul(val1a, matrix[1]));
      CTYPE res1a = add(mul(val0a, matrix[2]), mul(val1a, matrix[3]));
      state[basis] = res0a;
      state[basis + 1] = res1a;
    }
  } else {
    ITYPE state_index = 0;
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis_0 =
          (state_index & mask_low) + ((state_index & mask_high) << 1);
      ITYPE basis_1 = basis_0 + mask;
      CTYPE val0a = state[basis_0];
      CTYPE val0b = state[basis_0 + 1];
      CTYPE val1a = state[basis_1];
      CTYPE val1b = state[basis_1 + 1];

      CTYPE res0a = add(mul(val0a, matrix[0]), mul(val1a, matrix[1]));
      CTYPE res1b = add(mul(val0b, matrix[2]), mul(val1b, matrix[3]));
      CTYPE res1a = add(mul(val0a, matrix[2]), mul(val1a, matrix[3]));
      CTYPE res0b = add(mul(val0b, matrix[0]), mul(val1b, matrix[1]));

      state[basis_0] = res0a;
      state[basis_0 + 1] = res0b;
      state[basis_1] = res1a;
      state[basis_1 + 1] = res1b;
    }
  }
}

void multi_qubit_diagonal_matrix_gate(UINT __global *target_qubit_index_list,
                                      UINT target_qubit_index_count,
                                      CTYPE __global *diagonal_element,
                                      CTYPE __global *state, ITYPE dim,
                                      HEAP heap) {
  // matrix dim, mask, buffer
  ITYPE matrix_dim = 1ULL << target_qubit_index_count;
  ITYPE __global *matrix_mask_list = create_matrix_mask_list(
      target_qubit_index_list, target_qubit_index_count, heap);
  // insert index
  UINT __global *sorted_insert_index_list = create_sorted_ui_list(
      target_qubit_index_list, target_qubit_index_count, heap);
  // loop variables
  ITYPE loop_dim = dim >> target_qubit_index_count;
  ITYPE state_index;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create base index
    ITYPE basis_0 = state_index;
    for (UINT cursor = 0; cursor < target_qubit_index_count; cursor++) {
      UINT insert_index = sorted_insert_index_list[cursor];
      basis_0 = insert_zero_to_basis_index(basis_0, 1ULL << insert_index,
                                           insert_index);
    }

    // compute matrix-vector multiply
    for (ITYPE y = 0; y < matrix_dim; ++y) {
      state[basis_0 ^ matrix_mask_list[y]] =
          mul(state[basis_0 ^ matrix_mask_list[y]], diagonal_element[y]);
    }
  }

  free(heap, (uintptr_t)sorted_insert_index_list);
  free(heap, (uintptr_t)matrix_mask_list);
}

void multi_qubit_control_multi_qubit_diagonal_matrix_gate(
    UINT __global *control_qubit_index_list, UINT *control_value_list,
    UINT control_qubit_index_count, UINT __global *target_qubit_index_list,
    UINT target_qubit_index_count, CTYPE __global *diagonal_element,
    CTYPE __global *state, ITYPE dim, HEAP heap) {
  // matrix dim, mask, buffer
  ITYPE matrix_dim = 1ULL << target_qubit_index_count;
  ITYPE __global *matrix_mask_list = create_matrix_mask_list(
      target_qubit_index_list, target_qubit_index_count, heap);

  // insert index
  UINT insert_index_count =
      target_qubit_index_count + control_qubit_index_count;
  UINT __global *sorted_insert_index_list = create_sorted_ui_list_list(
      target_qubit_index_list, target_qubit_index_count,
      control_qubit_index_list, control_qubit_index_count, heap);

  // control mask
  ITYPE control_mask = create_control_mask(
      control_qubit_index_list, control_value_list, control_qubit_index_count);

  // loop varaibles
  ITYPE loop_dim =
      dim >> (target_qubit_index_count + control_qubit_index_count);
  ITYPE state_index;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create base index
    ITYPE basis_0 = state_index;
    for (UINT cursor = 0; cursor < insert_index_count; cursor++) {
      UINT insert_index = sorted_insert_index_list[cursor];
      basis_0 = insert_zero_to_basis_index(basis_0, 1ULL << insert_index,
                                           insert_index);
    }

    // flip control masks
    basis_0 ^= control_mask;

    // compute matrix mul
    for (ITYPE y = 0; y < matrix_dim; ++y) {
      state[basis_0 ^ matrix_mask_list[y]] =
          mul(state[basis_0 ^ matrix_mask_list[y]], diagonal_element[y]);
    }
  }

  free(heap, (uintptr_t)sorted_insert_index_list);
  free(heap, (uintptr_t)matrix_mask_list);
}

void single_qubit_diagonal_matrix_gate(UINT target_qubit_index,
                                       CTYPE diagonal_matrix[2],
                                       CTYPE __global *state, ITYPE dim) {
  single_qubit_diagonal_matrix_gate_parallel_unroll(
      target_qubit_index, diagonal_matrix, state, dim);
}

void single_qubit_diagonal_matrix_gate_parallel_unroll(UINT target_qubit_index,
                                                       CTYPE diagonal_matrix[2],
                                                       CTYPE __global *state,
                                                       ITYPE dim) {
  // loop variables
  ITYPE loop_dim = dim;
  ITYPE state_index;
  if (target_qubit_index == 0) {
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      state[state_index] = mul(state[state_index], diagonal_matrix[0]);
      state[state_index + 1] = mul(state[state_index + 1], diagonal_matrix[1]);
    }
  } else {
    ITYPE mask = 1ULL << target_qubit_index;

    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      int bitval = ((state_index & mask) != 0);
      state[state_index] = mul(state[state_index], diagonal_matrix[bitval]);
      state[state_index + 1] =
          mul(state[state_index + 1], diagonal_matrix[bitval]);
    }
  }
}

void single_qubit_phase_gate(UINT target_qubit_index, CTYPE phase,
                             CTYPE __global *state, ITYPE dim) {
  single_qubit_phase_gate_parallel_unroll(target_qubit_index, phase, state,
                                          dim);
}

void single_qubit_phase_gate_parallel_unroll(UINT target_qubit_index,
                                             CTYPE phase, CTYPE __global *state,
                                             ITYPE dim) {
  // target tmask
  ITYPE mask = 1ULL << target_qubit_index;
  ITYPE low_mask = mask - 1;
  ITYPE high_mask = ~low_mask;

  // loop varaibles
  ITYPE loop_dim = dim / 2;
  if (target_qubit_index == 0) {
    ITYPE state_index;
    for (state_index = 1; state_index < dim; state_index += 2) {
      state[state_index] = mul(state[state_index], phase);
    }
  } else {
    ITYPE state_index;
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis =
          (state_index & low_mask) + ((state_index & high_mask) << 1) + mask;
      state[basis] = mul(state[basis], phase);
      state[basis + 1] = mul(state[basis + 1], phase);
    }
  }
}

void single_qubit_phase_gate(UINT target_qubit_index, CTYPE phase,
                             CTYPE __global *state, ITYPE dim);

void S_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  single_qubit_phase_gate(target_qubit_index, double_to_complex(0, 1.0), state,
                          dim);
}
void Sdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  single_qubit_phase_gate(target_qubit_index, double_to_complex(0, -1.0), state,
                          dim);
}
void T_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  single_qubit_phase_gate(target_qubit_index,
                          double_to_complex(1.0 / sqrt(2.), 1.0 / sqrt(2.)),
                          state, dim);
}
void Tdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  single_qubit_phase_gate(target_qubit_index,
                          double_to_complex(1.0 / sqrt(2.), -1.0 / sqrt(2.)),
                          state, dim);
}
void sqrtX_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  single_qubit_dense_matrix_gate(target_qubit_index, SQRT_X_GATE_MATRIX, state,
                                 dim);
}
void sqrtXdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  single_qubit_dense_matrix_gate(target_qubit_index, SQRT_X_DAG_GATE_MATRIX,
                                 state, dim);
}
void sqrtY_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  single_qubit_dense_matrix_gate(target_qubit_index, SQRT_Y_GATE_MATRIX, state,
                                 dim);
}
void sqrtYdag_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  single_qubit_dense_matrix_gate(target_qubit_index, SQRT_Y_DAG_GATE_MATRIX,
                                 state, dim);
}

void CNOT_gate(UINT control_qubit_index, UINT target_qubit_index,
               CTYPE __global *state, ITYPE dim) {
  CNOT_gate_parallel_unroll(control_qubit_index, target_qubit_index, state,
                            dim);
}

void CNOT_gate_parallel_unroll(UINT control_qubit_index,
                               UINT target_qubit_index, CTYPE __global *state,
                               ITYPE dim) {
  ITYPE loop_dim = dim / 4;

  ITYPE target_mask = 1ULL << target_qubit_index;
  ITYPE control_mask = 1ULL << control_qubit_index;

  UINT min_qubit_index = get_min_ui(control_qubit_index, target_qubit_index);
  UINT max_qubit_index = get_max_ui(control_qubit_index, target_qubit_index);
  ITYPE min_qubit_mask = 1ULL << min_qubit_index;
  ITYPE max_qubit_mask = 1ULL << (max_qubit_index - 1);
  ITYPE low_mask = min_qubit_mask - 1;
  ITYPE mid_mask = (max_qubit_mask - 1) ^ low_mask;
  ITYPE high_mask = ~(max_qubit_mask - 1);

  ITYPE state_index = 0;
  if (target_qubit_index == 0) {
    // swap neighboring two basis
    for (state_index = 0; state_index < loop_dim; ++state_index) {
      ITYPE basis_index = ((state_index & mid_mask) << 1) +
                          ((state_index & high_mask) << 2) + control_mask;
      CTYPE temp = state[basis_index];
      state[basis_index] = state[basis_index + 1];
      state[basis_index + 1] = temp;
    }
  } else if (control_qubit_index == 0) {
    // no neighboring swap
    for (state_index = 0; state_index < loop_dim; ++state_index) {
      ITYPE basis_index_0 = (state_index & low_mask) +
                            ((state_index & mid_mask) << 1) +
                            ((state_index & high_mask) << 2) + control_mask;
      ITYPE basis_index_1 = basis_index_0 + target_mask;
      CTYPE temp = state[basis_index_0];
      state[basis_index_0] = state[basis_index_1];
      state[basis_index_1] = temp;
    }
  } else {
    // a,a+1 is swapped to a^m, a^m+1, respectively
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis_index_0 = (state_index & low_mask) +
                            ((state_index & mid_mask) << 1) +
                            ((state_index & high_mask) << 2) + control_mask;
      ITYPE basis_index_1 = basis_index_0 + target_mask;
      CTYPE temp0 = state[basis_index_0];
      CTYPE temp1 = state[basis_index_0 + 1];
      state[basis_index_0] = state[basis_index_1];
      state[basis_index_0 + 1] = state[basis_index_1 + 1];
      state[basis_index_1] = temp0;
      state[basis_index_1 + 1] = temp1;
    }
  }
}

void CZ_gate(UINT control_qubit_index, UINT target_qubit_index,
             CTYPE __global *state, ITYPE dim) {
  CZ_gate_parallel_unroll(control_qubit_index, target_qubit_index, state, dim);
}

void CZ_gate_parallel_unroll(UINT control_qubit_index, UINT target_qubit_index,
                             CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim / 4;

  ITYPE target_mask = 1ULL << target_qubit_index;
  ITYPE control_mask = 1ULL << control_qubit_index;

  UINT min_qubit_index = get_min_ui(control_qubit_index, target_qubit_index);
  UINT max_qubit_index = get_max_ui(control_qubit_index, target_qubit_index);
  ITYPE min_qubit_mask = 1ULL << min_qubit_index;
  ITYPE max_qubit_mask = 1ULL << (max_qubit_index - 1);
  ITYPE low_mask = min_qubit_mask - 1;
  ITYPE mid_mask = (max_qubit_mask - 1) ^ low_mask;
  ITYPE high_mask = ~(max_qubit_mask - 1);

  ITYPE mask = target_mask + control_mask;
  ITYPE state_index = 0;
  if (target_qubit_index == 0 || control_qubit_index == 0) {
    for (state_index = 0; state_index < loop_dim; ++state_index) {
      ITYPE basis_index = (state_index & low_mask) +
                          ((state_index & mid_mask) << 1) +
                          ((state_index & high_mask) << 2) + mask;
      state[basis_index] = mul(state[basis_index], double_to_complex(-1, 0));
    }
  } else {
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis_index = (state_index & low_mask) +
                          ((state_index & mid_mask) << 1) +
                          ((state_index & high_mask) << 2) + mask;
      state[basis_index] = mul(state[basis_index], double_to_complex(-1, 0));
      state[basis_index + 1] =
          mul(state[basis_index + 1], double_to_complex(-1, 0));
    }
  }
}

void H_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  H_gate_parallel_unroll(target_qubit_index, state, dim);
}

void H_gate_parallel_unroll(UINT target_qubit_index, CTYPE __global *state,
                            ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = (1ULL << target_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;
  double sqrt2inv = 1. / sqrt(2.);
  ITYPE state_index = 0;
  if (target_qubit_index == 0) {
    ITYPE basis_index;
    for (basis_index = 0; basis_index < dim; basis_index += 2) {
      CTYPE temp0 = state[basis_index];
      CTYPE temp1 = state[basis_index + 1];
      state[basis_index] =
          mul(add(temp0, temp1), double_to_complex(sqrt2inv, 0));
      state[basis_index + 1] =
          mul(sub(temp0, temp1), double_to_complex(sqrt2inv, 0));
    }
  } else {
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis_index_0 =
          (state_index & mask_low) + ((state_index & mask_high) << 1);
      ITYPE basis_index_1 = basis_index_0 + mask;
      CTYPE temp_a0 = state[basis_index_0];
      CTYPE temp_a1 = state[basis_index_1];
      CTYPE temp_b0 = state[basis_index_0 + 1];
      CTYPE temp_b1 = state[basis_index_1 + 1];
      state[basis_index_0] =
          mul(add(temp_a0, temp_a1), double_to_complex(sqrt2inv, 0));
      state[basis_index_1] =
          mul(sub(temp_a0, temp_a1), double_to_complex(sqrt2inv, 0));
      state[basis_index_0 + 1] =
          mul(add(temp_b0, temp_b1), double_to_complex(sqrt2inv, 0));
      state[basis_index_1 + 1] =
          mul(sub(temp_b0, temp_b1), double_to_complex(sqrt2inv, 0));
    }
  }
}

void SWAP_gate(UINT target_qubit_index_0, UINT target_qubit_index_1,
               CTYPE __global *state, ITYPE dim) {
  SWAP_gate_parallel_unroll(target_qubit_index_0, target_qubit_index_1, state,
                            dim);
}

void SWAP_gate_parallel_unroll(UINT target_qubit_index_0,
                               UINT target_qubit_index_1, CTYPE __global *state,
                               ITYPE dim) {
  ITYPE loop_dim = dim / 4;

  ITYPE mask_0 = 1ULL << target_qubit_index_0;
  ITYPE mask_1 = 1ULL << target_qubit_index_1;
  ITYPE mask = mask_0 + mask_1;

  UINT min_qubit_index = get_min_ui(target_qubit_index_0, target_qubit_index_1);
  UINT max_qubit_index = get_max_ui(target_qubit_index_0, target_qubit_index_1);
  ITYPE min_qubit_mask = 1ULL << min_qubit_index;
  ITYPE max_qubit_mask = 1ULL << (max_qubit_index - 1);
  ITYPE low_mask = min_qubit_mask - 1;
  ITYPE mid_mask = (max_qubit_mask - 1) ^ low_mask;
  ITYPE high_mask = ~(max_qubit_mask - 1);

  ITYPE state_index = 0;
  if (target_qubit_index_0 == 0 || target_qubit_index_1 == 0) {
    for (state_index = 0; state_index < loop_dim; ++state_index) {
      ITYPE basis_index_0 = (state_index & low_mask) +
                            ((state_index & mid_mask) << 1) +
                            ((state_index & high_mask) << 2) + mask_0;
      ITYPE basis_index_1 = basis_index_0 ^ mask;
      CTYPE temp = state[basis_index_0];
      state[basis_index_0] = state[basis_index_1];
      state[basis_index_1] = temp;
    }
  } else {
    // a,a+1 is swapped to a^m, a^m+1, respectively
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis_index_0 = (state_index & low_mask) +
                            ((state_index & mid_mask) << 1) +
                            ((state_index & high_mask) << 2) + mask_0;
      ITYPE basis_index_1 = basis_index_0 ^ mask;
      CTYPE temp0 = state[basis_index_0];
      CTYPE temp1 = state[basis_index_0 + 1];
      state[basis_index_0] = state[basis_index_1];
      state[basis_index_0 + 1] = state[basis_index_1 + 1];
      state[basis_index_1] = temp0;
      state[basis_index_1 + 1] = temp1;
    }
  }
}

void X_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  X_gate_parallel_unroll(target_qubit_index, state, dim);
}

void X_gate_parallel_unroll(UINT target_qubit_index, CTYPE __global *state,
                            ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = (1ULL << target_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;
  ITYPE state_index = 0;
  if (target_qubit_index == 0) {
    ITYPE basis_index = 0;
    for (basis_index = 0; basis_index < dim; basis_index += 2) {
      CTYPE temp = state[basis_index];
      state[basis_index] = state[basis_index + 1];
      state[basis_index + 1] = temp;
    }
  } else {
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis_index_0 =
          (state_index & mask_low) + ((state_index & mask_high) << 1);
      ITYPE basis_index_1 = basis_index_0 + mask;
      CTYPE temp0 = state[basis_index_0];
      CTYPE temp1 = state[basis_index_0 + 1];
      state[basis_index_0] = state[basis_index_1];
      state[basis_index_0 + 1] = state[basis_index_1 + 1];
      state[basis_index_1] = temp0;
      state[basis_index_1 + 1] = temp1;
    }
  }
}

void Y_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  Y_gate_parallel_unroll(target_qubit_index, state, dim);
}

void Y_gate_parallel_unroll(UINT target_qubit_index, CTYPE __global *state,
                            ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = (1ULL << target_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;
  ITYPE state_index = 0;
  CTYPE imag = double_to_complex(0, 1.0);
  CTYPE neg_imag = double_to_complex(0, -1.0);
  if (target_qubit_index == 0) {
    ITYPE basis_index;
    for (basis_index = 0; basis_index < dim; basis_index += 2) {
      CTYPE temp0 = state[basis_index];
      state[basis_index] = mul(neg_imag, state[basis_index + 1]);
      state[basis_index + 1] = mul(imag, temp0);
    }
  } else {
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis_index_0 =
          (state_index & mask_low) + ((state_index & mask_high) << 1);
      ITYPE basis_index_1 = basis_index_0 + mask;
      CTYPE temp0 = state[basis_index_0];
      CTYPE temp1 = state[basis_index_0 + 1];
      state[basis_index_0] = mul(neg_imag, state[basis_index_1]);
      state[basis_index_0 + 1] = mul(neg_imag, state[basis_index_1 + 1]);
      state[basis_index_1] = mul(imag, temp0);
      state[basis_index_1 + 1] = mul(imag, temp1);
    }
  }
}

void Z_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  Z_gate_parallel_unroll(target_qubit_index, state, dim);
}

void Z_gate_parallel_unroll(UINT target_qubit_index, CTYPE __global *state,
                            ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = (1ULL << target_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;
  ITYPE state_index = 0;
  if (target_qubit_index == 0) {
    for (state_index = 1; state_index < dim; state_index += 2) {
      state[state_index] = mul(state[state_index], double_to_complex(-1, 0));
    }
  } else {
    for (state_index = 0; state_index < loop_dim; state_index += 2) {
      ITYPE basis_index =
          (state_index & mask_low) + ((state_index & mask_high) << 1) + mask;
      state[basis_index] = mul(state[basis_index], double_to_complex(-1, 0));
      state[basis_index + 1] =
          mul(state[basis_index + 1], double_to_complex(-1, 0));
    }
  }
}

void P0_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  P0_gate_parallel(target_qubit_index, state, dim);
}

void P1_gate(UINT target_qubit_index, CTYPE __global *state, ITYPE dim) {
  P1_gate_parallel(target_qubit_index, state, dim);
}

void P0_gate_parallel(UINT target_qubit_index, CTYPE __global *state,
                      ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = (1ULL << target_qubit_index);
  ITYPE low_mask = mask - 1;
  ITYPE high_mask = ~low_mask;

  ITYPE state_index;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE temp_index =
        (state_index & low_mask) + ((state_index & high_mask) << 1) + mask;
    state[temp_index] = double_to_complex(0, 0);
  }
}

void P1_gate_parallel(UINT target_qubit_index, CTYPE __global *state,
                      ITYPE dim) {
  ITYPE loop_dim = dim / 2;
  ITYPE mask = (1ULL << target_qubit_index);
  ITYPE low_mask = mask - 1;
  ITYPE high_mask = ~low_mask;

  ITYPE state_index;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    ITYPE temp_index =
        (state_index & low_mask) + ((state_index & high_mask) << 1);
    state[temp_index] = double_to_complex(0, 0);
  }
}

void normalize(double squared_norm, CTYPE __global *state, ITYPE dim) {
  ITYPE loop_dim = dim;
  double normalize_factor = sqrt(1. / squared_norm);
  ITYPE state_index;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    state[state_index] =
        mul(state[state_index], double_to_complex(normalize_factor, 0));
  }
}

void normalize_single_thread(double squared_norm, CTYPE __global *state,
                             ITYPE dim) {
  ITYPE loop_dim = dim;
  double normalize_factor = sqrt(1. / squared_norm);
  ITYPE state_index;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    state[state_index] =
        mul(state[state_index], double_to_complex(normalize_factor, 0));
  }
}

void state_add(CTYPE __global *state_added, CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  for (index = 0; index < dim; ++index) {
    state[index] = add(state[index], state_added[index]);
  }
}

void state_add_with_coef(CTYPE coef, CTYPE __global *state_added,
                         CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  for (index = 0; index < dim; ++index) {
    state[index] = add(state[index], mul(coef, state_added[index]));
  }
}

void state_add_with_coef_single_thread(CTYPE coef, CTYPE __global *state_added,
                                       CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  for (index = 0; index < dim; ++index) {
    state[index] = add(state[index], mul(coef, state_added[index]));
  }
}

void state_multiply(CTYPE coef, CTYPE __global *state, ITYPE dim) {
  ITYPE index;
  for (index = 0; index < dim; ++index) {
    state[index] = mul(state[index], coef);
  }
}

/**
 * perform multi_qubit_Pauli_gate with XZ mask.
 *
 * This function assumes bit_flip_mask is not 0, i.e., at least one bit is
 * flipped. If no bit is flipped, use multi_qubit_Pauli_gate_Z_mask. This
 * function update the quantum state with Pauli operation. bit_flip_mask,
 * phase_flip_mask, global_phase_90rot_count, and pivot_qubit_index must be
 * computed before calling this function. See get_masks_from_*_list for the
 * above four arguemnts.
 */
void multi_qubit_Pauli_gate_XZ_mask(ITYPE bit_flip_mask, ITYPE phase_flip_mask,
                                    UINT global_phase_90rot_count,
                                    UINT pivot_qubit_index,
                                    CTYPE __global *state, ITYPE dim);
void multi_qubit_Pauli_rotation_gate_XZ_mask(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, double angle, CTYPE __global *state, ITYPE dim);
void multi_qubit_Pauli_gate_Z_mask(ITYPE phase_flip_mask, CTYPE __global *state,
                                   ITYPE dim);
void multi_qubit_Pauli_rotation_gate_Z_mask(ITYPE phase_flip_mask, double angle,
                                            CTYPE __global *state, ITYPE dim);

void multi_qubit_Pauli_gate_XZ_mask_single_thread(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, CTYPE __global *state, ITYPE dim);
void multi_qubit_Pauli_rotation_gate_XZ_mask_single_thread(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, double angle, CTYPE __global *state, ITYPE dim);
void multi_qubit_Pauli_gate_Z_mask_single_thread(ITYPE phase_flip_mask,
                                                 CTYPE __global *state,
                                                 ITYPE dim);
void multi_qubit_Pauli_rotation_gate_Z_mask_single_thread(ITYPE phase_flip_mask,
                                                          double angle,
                                                          CTYPE __global *state,
                                                          ITYPE dim);

void multi_qubit_Pauli_gate_XZ_mask(ITYPE bit_flip_mask, ITYPE phase_flip_mask,
                                    UINT global_phase_90rot_count,
                                    UINT pivot_qubit_index,
                                    CTYPE __global *state, ITYPE dim) {
  // loop varaibles
  ITYPE loop_dim = dim / 2;
  ITYPE state_index;

  ITYPE mask = (1ULL << pivot_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create base index
    ITYPE basis_0 = (state_index & mask_low) + ((state_index & mask_high) << 1);

    // gather index
    ITYPE basis_1 = basis_0 ^ bit_flip_mask;

    // determine sign
    UINT sign_0 = count_population(basis_0 & phase_flip_mask) % 2;
    UINT sign_1 = count_population(basis_1 & phase_flip_mask) % 2;

    // fetch values
    CTYPE cval_0 = state[basis_0];
    CTYPE cval_1 = state[basis_1];

    // set values
    state[basis_0] =
        mul(cval_1, PHASE_M90ROT[(global_phase_90rot_count + sign_0 * 2) % 4]);
    state[basis_1] =
        mul(cval_0, PHASE_M90ROT[(global_phase_90rot_count + sign_1 * 2) % 4]);
  }
}
void multi_qubit_Pauli_rotation_gate_XZ_mask(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, double angle, CTYPE __global *state, ITYPE dim) {
  // loop varaibles
  ITYPE loop_dim = dim / 2;
  ITYPE state_index;

  ITYPE mask = (1ULL << pivot_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;

  // coefs
  double cosval = cos(angle / 2);
  double sinval = sin(angle / 2);

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create base index
    ITYPE basis_0 = (state_index & mask_low) + ((state_index & mask_high) << 1);

    // gather index
    ITYPE basis_1 = basis_0 ^ bit_flip_mask;

    // determine parity
    int bit_parity_0 = count_population(basis_0 & phase_flip_mask) % 2;
    int bit_parity_1 = count_population(basis_1 & phase_flip_mask) % 2;

    // fetch values
    CTYPE cval_0 = state[basis_0];
    CTYPE cval_1 = state[basis_1];

    // set values
    state[basis_0] = add(
        mul(double_to_complex(cosval, 0), cval_0),
        mul(mul(mul(double_to_complex(0, 1.0), double_to_complex(sinval, 0)),
                cval_1),
            PHASE_M90ROT[(global_phase_90rot_count + bit_parity_0 * 2) % 4]));
    state[basis_1] = add(
        mul(double_to_complex(cosval, 0), cval_1),
        mul(mul(mul(double_to_complex(0, 1.0), double_to_complex(sinval, 0)),
                cval_0),
            PHASE_M90ROT[(global_phase_90rot_count + bit_parity_1 * 2) % 4]));
  }
}

void multi_qubit_Pauli_gate_Z_mask(ITYPE phase_flip_mask, CTYPE __global *state,
                                   ITYPE dim) {
  // loop varaibles
  ITYPE loop_dim = dim;
  ITYPE state_index;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // determine parity
    int bit_parity = count_population(state_index & phase_flip_mask) % 2;

    // set values
    if (bit_parity % 2 == 1) {
      state[state_index] = mul(state[state_index], double_to_complex(-1.0, 0));
    }
  }
}

void multi_qubit_Pauli_rotation_gate_Z_mask(ITYPE phase_flip_mask, double angle,
                                            CTYPE __global *state, ITYPE dim) {
  // loop variables
  ITYPE loop_dim = dim;
  ITYPE state_index;

  // coefs
  double cosval = cos(angle / 2);
  double sinval = sin(angle / 2);

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // determine sign
    int bit_parity = count_population(state_index & phase_flip_mask) % 2;
    int sign = 1 - 2 * bit_parity;

    // set value
    state[state_index] =
        mul(state[state_index],
            add(double_to_complex(cosval, 0),
                mul(mul(double_to_complex(sign, 0), double_to_complex(0, 1.0)),
                    double_to_complex(sinval, 0))));
  }
}

void multi_qubit_Pauli_gate_partial_list(UINT __global *target_qubit_index_list,
                                         UINT *Pauli_operator_type_list,
                                         UINT target_qubit_index_count,
                                         CTYPE __global *state, ITYPE dim) {
  // create pauli mask and call function
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_partial_list(
      target_qubit_index_list, Pauli_operator_type_list,
      target_qubit_index_count, &bit_flip_mask, &phase_flip_mask,
      &global_phase_90rot_count, &pivot_qubit_index);
  if (bit_flip_mask == 0) {
    multi_qubit_Pauli_gate_Z_mask(phase_flip_mask, state, dim);
  } else {
    multi_qubit_Pauli_gate_XZ_mask(bit_flip_mask, phase_flip_mask,
                                   global_phase_90rot_count, pivot_qubit_index,
                                   state, dim);
  }
}

void multi_qubit_Pauli_gate_whole_list(UINT *Pauli_operator_type_list,
                                       UINT qubit_count, CTYPE __global *state,
                                       ITYPE dim) {
  // create pauli mask and call function
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_whole_list(Pauli_operator_type_list, qubit_count,
                             &bit_flip_mask, &phase_flip_mask,
                             &global_phase_90rot_count, &pivot_qubit_index);
  if (bit_flip_mask == 0) {
    multi_qubit_Pauli_gate_Z_mask(phase_flip_mask, state, dim);
  } else {
    multi_qubit_Pauli_gate_XZ_mask(bit_flip_mask, phase_flip_mask,
                                   global_phase_90rot_count, pivot_qubit_index,
                                   state, dim);
  }
}

void multi_qubit_Pauli_rotation_gate_partial_list(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, double angle, CTYPE __global *state,
    ITYPE dim) {
  // create pauli mask and call function
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_partial_list(
      target_qubit_index_list, Pauli_operator_type_list,
      target_qubit_index_count, &bit_flip_mask, &phase_flip_mask,
      &global_phase_90rot_count, &pivot_qubit_index);
  if (bit_flip_mask == 0) {
    multi_qubit_Pauli_rotation_gate_Z_mask(phase_flip_mask, angle, state, dim);
  } else {
    multi_qubit_Pauli_rotation_gate_XZ_mask(
        bit_flip_mask, phase_flip_mask, global_phase_90rot_count,
        pivot_qubit_index, angle, state, dim);
  }
}

void multi_qubit_Pauli_rotation_gate_whole_list(UINT *Pauli_operator_type_list,
                                                UINT qubit_count, double angle,
                                                CTYPE __global *state,
                                                ITYPE dim) {
  // create pauli mask and call function
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_whole_list(Pauli_operator_type_list, qubit_count,
                             &bit_flip_mask, &phase_flip_mask,
                             &global_phase_90rot_count, &pivot_qubit_index);
  if (bit_flip_mask == 0) {
    multi_qubit_Pauli_rotation_gate_Z_mask(phase_flip_mask, angle, state, dim);
  } else {
    multi_qubit_Pauli_rotation_gate_XZ_mask(
        bit_flip_mask, phase_flip_mask, global_phase_90rot_count,
        pivot_qubit_index, angle, state, dim);
  }
}

void multi_qubit_Pauli_gate_XZ_mask_single_thread(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, CTYPE __global *state, ITYPE dim) {
  // loop varaibles
  ITYPE loop_dim = dim / 2;
  ITYPE state_index;

  ITYPE mask = (1ULL << pivot_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;

  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create base index
    ITYPE basis_0 = (state_index & mask_low) + ((state_index & mask_high) << 1);

    // gather index
    ITYPE basis_1 = basis_0 ^ bit_flip_mask;

    // determine sign
    UINT sign_0 = count_population(basis_0 & phase_flip_mask) % 2;
    UINT sign_1 = count_population(basis_1 & phase_flip_mask) % 2;

    // fetch values
    CTYPE cval_0 = state[basis_0];
    CTYPE cval_1 = state[basis_1];

    // set values
    state[basis_0] =
        mul(cval_1, PHASE_M90ROT[(global_phase_90rot_count + sign_0 * 2) % 4]);
    state[basis_1] =
        mul(cval_0, PHASE_M90ROT[(global_phase_90rot_count + sign_1 * 2) % 4]);
  }
}

void multi_qubit_Pauli_rotation_gate_XZ_mask_single_thread(
    ITYPE bit_flip_mask, ITYPE phase_flip_mask, UINT global_phase_90rot_count,
    UINT pivot_qubit_index, double angle, CTYPE __global *state, ITYPE dim) {
  // loop varaibles
  ITYPE loop_dim = dim / 2;
  ITYPE state_index;

  ITYPE mask = (1ULL << pivot_qubit_index);
  ITYPE mask_low = mask - 1;
  ITYPE mask_high = ~mask_low;

  // coefs
  double cosval = cos(angle / 2);
  double sinval = sin(angle / 2);
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // create base index
    ITYPE basis_0 = (state_index & mask_low) + ((state_index & mask_high) << 1);

    // gather index
    ITYPE basis_1 = basis_0 ^ bit_flip_mask;

    // determine parity
    int bit_parity_0 = count_population(basis_0 & phase_flip_mask) % 2;
    int bit_parity_1 = count_population(basis_1 & phase_flip_mask) % 2;

    // fetch values
    CTYPE cval_0 = state[basis_0];
    CTYPE cval_1 = state[basis_1];

    // set values
    state[basis_0] = add(
        mul(double_to_complex(cosval, 0), cval_0),
        mul(mul(mul(double_to_complex(0, 1.0), double_to_complex(sinval, 0)),
                cval_1),
            PHASE_M90ROT[(global_phase_90rot_count + bit_parity_0 * 2) % 4]));
    state[basis_1] = add(
        mul(double_to_complex(cosval, 0), cval_1),
        mul(mul(mul(double_to_complex(0, 1.0), double_to_complex(sinval, 0)),
                cval_0),
            PHASE_M90ROT[(global_phase_90rot_count + bit_parity_1 * 2) % 4]));
  }
}

void multi_qubit_Pauli_gate_Z_mask_single_thread(ITYPE phase_flip_mask,
                                                 CTYPE __global *state,
                                                 ITYPE dim) {
  // loop varaibles
  ITYPE loop_dim = dim;
  ITYPE state_index;
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // determine parity
    int bit_parity = count_population(state_index & phase_flip_mask) % 2;

    // set values
    if (bit_parity % 2 == 1) {
      state[state_index] = mul(state[state_index], double_to_complex(-1.0, 0));
    }
  }
}

void multi_qubit_Pauli_rotation_gate_Z_mask_single_thread(ITYPE phase_flip_mask,
                                                          double angle,
                                                          CTYPE __global *state,
                                                          ITYPE dim) {
  // loop variables
  ITYPE loop_dim = dim;
  ITYPE state_index;

  // coefs
  double cosval = cos(angle / 2);
  double sinval = sin(angle / 2);
  for (state_index = 0; state_index < loop_dim; ++state_index) {
    // determine sign
    int bit_parity = count_population(state_index & phase_flip_mask) % 2;
    int sign = 1 - 2 * bit_parity;

    // set value
    state[state_index] =
        mul(state[state_index],
            add(double_to_complex(cosval, 0),
                mul(mul(double_to_complex(sign, 0), double_to_complex(0, 1.0)),
                    double_to_complex(sinval, 0))));
  }
}

void multi_qubit_Pauli_gate_partial_list_single_thread(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, CTYPE __global *state, ITYPE dim) {
  // create pauli mask and call function
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_partial_list(
      target_qubit_index_list, Pauli_operator_type_list,
      target_qubit_index_count, &bit_flip_mask, &phase_flip_mask,
      &global_phase_90rot_count, &pivot_qubit_index);
  if (bit_flip_mask == 0) {
    multi_qubit_Pauli_gate_Z_mask_single_thread(phase_flip_mask, state, dim);
  } else {
    multi_qubit_Pauli_gate_XZ_mask_single_thread(bit_flip_mask, phase_flip_mask,
                                                 global_phase_90rot_count,
                                                 pivot_qubit_index, state, dim);
  }
}

void multi_qubit_Pauli_gate_whole_list_single_thread(
    UINT *Pauli_operator_type_list, UINT qubit_count, CTYPE __global *state,
    ITYPE dim) {
  // create pauli mask and call function
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_whole_list(Pauli_operator_type_list, qubit_count,
                             &bit_flip_mask, &phase_flip_mask,
                             &global_phase_90rot_count, &pivot_qubit_index);
  if (bit_flip_mask == 0) {
    multi_qubit_Pauli_gate_Z_mask_single_thread(phase_flip_mask, state, dim);
  } else {
    multi_qubit_Pauli_gate_XZ_mask_single_thread(bit_flip_mask, phase_flip_mask,
                                                 global_phase_90rot_count,
                                                 pivot_qubit_index, state, dim);
  }
}

void multi_qubit_Pauli_rotation_gate_partial_list_single_thread(
    UINT __global *target_qubit_index_list, UINT *Pauli_operator_type_list,
    UINT target_qubit_index_count, double angle, CTYPE __global *state,
    ITYPE dim) {
  // create pauli mask and call function
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_partial_list(
      target_qubit_index_list, Pauli_operator_type_list,
      target_qubit_index_count, &bit_flip_mask, &phase_flip_mask,
      &global_phase_90rot_count, &pivot_qubit_index);
  if (bit_flip_mask == 0) {
    multi_qubit_Pauli_rotation_gate_Z_mask_single_thread(phase_flip_mask, angle,
                                                         state, dim);
  } else {
    multi_qubit_Pauli_rotation_gate_XZ_mask_single_thread(
        bit_flip_mask, phase_flip_mask, global_phase_90rot_count,
        pivot_qubit_index, angle, state, dim);
  }
}

void multi_qubit_Pauli_rotation_gate_whole_list_single_thread(
    UINT *Pauli_operator_type_list, UINT qubit_count, double angle,
    CTYPE __global *state, ITYPE dim) {
  // create pauli mask and call function
  ITYPE bit_flip_mask = 0;
  ITYPE phase_flip_mask = 0;
  UINT global_phase_90rot_count = 0;
  UINT pivot_qubit_index = 0;
  get_Pauli_masks_whole_list(Pauli_operator_type_list, qubit_count,
                             &bit_flip_mask, &phase_flip_mask,
                             &global_phase_90rot_count, &pivot_qubit_index);
  if (bit_flip_mask == 0) {
    multi_qubit_Pauli_rotation_gate_Z_mask_single_thread(phase_flip_mask, angle,
                                                         state, dim);
  } else {
    multi_qubit_Pauli_rotation_gate_XZ_mask_single_thread(
        bit_flip_mask, phase_flip_mask, global_phase_90rot_count,
        pivot_qubit_index, angle, state, dim);
  }
}

void single_qubit_Pauli_gate(UINT target_qubit_index, UINT Pauli_operator_type,
                             CTYPE __global *state, ITYPE dim) {
  switch (Pauli_operator_type) {
  case 0:
    break;
  case 1:
    X_gate(target_qubit_index, state, dim);
    break;
  case 2:
    Y_gate(target_qubit_index, state, dim);
    break;
  case 3:
    Z_gate(target_qubit_index, state, dim);
    break;
  default:
    printf("ERROR: invalid Pauli operation is called");
  }
}

void single_qubit_Pauli_rotation_gate(UINT target_qubit_index,
                                      UINT Pauli_operator_index, double angle,
                                      CTYPE __global *state, ITYPE dim) {
  switch (Pauli_operator_index) {
  case 0:
    break;
  case 1:
    RX_gate(target_qubit_index, angle, state, dim);
    break;
  case 2:
    RY_gate(target_qubit_index, angle, state, dim);
    break;
  case 3:
    RZ_gate(target_qubit_index, angle, state, dim);
    break;
  default:
    printf("ERROR: Invalid Pauli operation is called");
  }
}

void RX_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
             ITYPE dim) {
  CTYPE matrix[4];
  double c, s;
  c = cos(angle / 2);
  s = sin(angle / 2);
  matrix[0] = double_to_complex(c, 0);
  matrix[1] = mul(double_to_complex(0, 1.0), double_to_complex(s, 0));
  matrix[2] = mul(double_to_complex(0, 1.0), double_to_complex(s, 0));
  matrix[3] = double_to_complex(c, 0);
  single_qubit_dense_matrix_gate_non_constant(target_qubit_index, matrix, state,
                                              dim);
}

void RY_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
             ITYPE dim) {
  CTYPE matrix[4];
  double c, s;
  c = cos(angle / 2);
  s = sin(angle / 2);
  matrix[0] = double_to_complex(c, 0);
  matrix[1] = double_to_complex(s, 0);
  matrix[2] = double_to_complex(-s, 0);
  matrix[3] = double_to_complex(c, 0);
  single_qubit_dense_matrix_gate_non_constant(target_qubit_index, matrix, state,
                                              dim);
}

void RZ_gate(UINT target_qubit_index, double angle, CTYPE __global *state,
             ITYPE dim) {
  CTYPE diagonal_matrix[2];
  double c, s;
  c = cos(angle / 2);
  s = sin(angle / 2);
  diagonal_matrix[0] =
      add(double_to_complex(c, 0),
          mul(double_to_complex(0, 1.0), double_to_complex(s, 0)));
  diagonal_matrix[1] =
      sub(double_to_complex(c, 0),
          mul(double_to_complex(0, 1.0), double_to_complex(s, 0)));
  single_qubit_diagonal_matrix_gate(target_qubit_index, diagonal_matrix, state,
                                    dim);
}

void CUz_gate(double angle, UINT c_bit, UINT t_bit, CTYPE __global *psi,
              ITYPE dim) {
  // ITYPE i;
  ITYPE j;
  ITYPE head, body, tail;
  ITYPE basis00, basis11;
  ITYPE left, right;
  // CTYPE U_gate[2][2];
  /*
  for(i = 0; i < 2; ++i)
      for(j = 0; j < 2; ++j)
          U_gate[i][j] = cos(angle) * PAULI_MATRIX[0][i][j] + sin(angle)
  * double_to_complex(0, 1.0)
  * PAULI_MATRIX[3][i][j];
  */

  if (t_bit > c_bit) {
    left = t_bit;
    right = c_bit;
  } else {
    left = c_bit;
    right = t_bit;
  }

  for (j = 0; j < dim / 4; ++j) {
    head = j >> (left - 1);
    body = (j & ((1ULL << (left - 1)) - 1)) >> right; // (j % 2^(k-1)) >> i
    tail = j & ((1ULL << right) - 1);                 // j%(2^i)

    basis00 = (head << (left + 1)) + (body << (right + 1)) + tail;
    basis11 = basis00 + (1ULL << c_bit) + (1ULL << t_bit);

    psi[basis11] = add(
        mul(double_to_complex(cos(angle), 0), psi[basis11]),
        mul(mul(double_to_complex(sin(angle), 0), double_to_complex(0, 1.0)),
            psi[basis11]));
  }
}

// qft: apply qft from k th to k+Nbits th
void qft(UINT k, UINT Nbits, int doSWAP, CTYPE __global *psi, ITYPE dim) {
  UINT i, j;
  for (i = 1; i < Nbits; ++i) {
    single_qubit_dense_matrix_gate(Nbits - i + k, HADAMARD_MATRIX, psi, dim);
    // printf("Had %d\n",Nbits-i+k);
    for (j = 0; j < i; ++j) {
      // printf("CUZ %d %d %.5f*PI\n", Nbits-i-1+k,
      // Nbits-j-1+k, 1.0/(1ULL<<(i-j)));
      // CUz_gate(1.0*PI/(1ULL<<(i-j)), Nbits-i-1+k, Nbits-j-1+k, psi);
      CUz_gate(1.0 * PI / (1ULL << (i - j)), Nbits - i - 1 + k,
               Nbits - j - 1 + k, psi, dim);
    }
  }
  single_qubit_dense_matrix_gate(k, HADAMARD_MATRIX, psi, dim);
  // printf("Had %d\n",k);
  if (doSWAP) {
    for (i = 0; i < Nbits / 2; ++i) {
      SWAP_gate(i + k, Nbits - 1 - i + k, psi, dim);
      // printf("SWAP %d %d\n",i+k,Nbits-1-i+k);
    }
  }
}

// inverse_qft: apply inverse_qft from k th to k+Nbits th
void inverse_qft(UINT k, UINT Nbits, int doSWAP, CTYPE __global *psi,
                 ITYPE dim) {
  UINT i, j;
  if (doSWAP) {
    for (i = 0; i < Nbits / 2; ++i) {
      SWAP_gate(i + k, Nbits - 1 - i + k, psi, dim);
      // printf("SWAP %d %d\n",i+k,Nbits-1-i+k);
    }
  }
  for (i = 0; i < Nbits; ++i) {
    single_qubit_dense_matrix_gate(i + k, HADAMARD_MATRIX, psi, dim);
    // printf("Had %d\n",i+k);
    for (j = i + 1; j < Nbits; ++j) {
      // printf("CUZ %d %d %.5f*PI\n",i+k,j+k,-1.0/(1ULL<<(j-i)));
      // CUz_gate(-1.0*PI/(1ULL<<(j-i)), i+k, j+k, psi);
      CUz_gate(-1.0 * PI / (1ULL << (j - i)), i + k, j + k, psi, dim);
    }
  }
}

void reflection_gate(CTYPE __global *reflection_state, CTYPE __global *state,
                     ITYPE dim) {
  CTYPE coef = state_inner_product(reflection_state, state, dim);

  for (ITYPE state_index = 0; state_index < dim; ++state_index) {
    state[state_index] = sub(mul(mul(double_to_complex(2.0, 0), coef),
                                 reflection_state[state_index]),
                             state[state_index]);
  }
}

// FALTAN
/*
void reversible_boolean_gate(UINT __global *target_qubit_index_list,
    UINT target_qubit_index_count,
    ITYPE (*function_ptr)(ITYPE, ITYPE), CTYPE __global* state, ITYPE dim, HEAP
heap) {
    // matrix dim, mask, buffer
     ITYPE matrix_dim = 1ULL << target_qubit_index_count;
    ITYPE __global* matrix_mask_list = create_matrix_mask_list(
        target_qubit_index_list, target_qubit_index_count, heap);

    // insert index
    UINT __global* sorted_insert_index_list = create_sorted_ui_list(
        target_qubit_index_list, target_qubit_index_count, heap);

    // loop variables
     ITYPE loop_dim = dim >> target_qubit_index_count;

    CTYPE* buffer = (CTYPE*)malloc(heap, (size_t)(sizeof(CTYPE) * matrix_dim));
    ITYPE state_index;
    for (state_index = 0; state_index < loop_dim; ++state_index) {
        // create base index
        ITYPE basis_0 = state_index;
        for (UINT cursor = 0; cursor < target_qubit_index_count; cursor++) {
            UINT insert_index = sorted_insert_index_list[cursor];
            basis_0 = insert_zero_to_basis_index(
                basis_0, 1ULL << insert_index, insert_index);
        }

        // compute matrix-vector multiply

        for (ITYPE x = 0; x < matrix_dim; ++x) {
            ITYPE y = function_ptr(x, matrix_dim);
            buffer[y] = state[basis_0 ^ matrix_mask_list[x]];
        }

        // set result
        for (ITYPE y = 0; y < matrix_dim; ++y) {
            state[basis_0 ^ matrix_mask_list[y]] = buffer[y];
        }
    }
    free(heap, buffer);
    free(heap, (UINT*)sorted_insert_index_list);
    free(heap, (ITYPE*)matrix_mask_list);
}*/