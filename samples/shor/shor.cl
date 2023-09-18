void CSWAP_gate(UINT control_qubit_index, UINT target_qubit_index1,
                UINT target_qubit_index2, CTYPE __global *state, ITYPE dim,
                HEAP heap) {
    
    UINT __global *target_qubit_index_list =
        (UINT __global *)malloc(heap, sizeof(UINT) * 2);
    CTYPE swap_matrix[16] = {
      {1, 0}, {0, 0}, {0, 0}, {0, 0},
      {0, 0}, {0, 0}, {1, 0}, {0, 0},
      {0, 0}, {1, 0}, {0, 0}, {0, 0},
      {0, 0}, {0, 0}, {0, 0}, {1, 0}
    };

    target_qubit_index_list[0] = target_qubit_index1;
    target_qubit_index_list[1] = target_qubit_index2;

    UINT control_value = 1;

    single_qubit_control_multi_qubit_dense_matrix_gate(
      control_qubit_index, control_value,
      target_qubit_index_list, 2,
      swap_matrix, state, dim, heap);
    
    free(heap, (uintptr_t)target_qubit_index_list);
}

void c_amod15(UINT control_qubit_index, CTYPE __global *state, int a, int power,
            int n_counts, ITYPE dim, HEAP heap) {
    // if a not in [2,4,7,8,11,13]:
    //     raise ValueError("'a' must be 2,4,7,8,11 or 13")

    // Cambiar
    for(int i = 0; i < power; i++){
        if (a == 2 || a == 13){
            CSWAP_gate(control_qubit_index, n_counts + 2, n_counts + 3, state, dim, heap);
            CSWAP_gate(control_qubit_index, n_counts + 1, n_counts + 2, state, dim, heap);
            CSWAP_gate(control_qubit_index, n_counts + 0, n_counts + 1, state, dim, heap);
        } if (a == 7 || a == 8){
            CSWAP_gate(control_qubit_index, n_counts + 0, n_counts + 1, state, dim, heap);
            CSWAP_gate(control_qubit_index, n_counts + 1, n_counts + 2, state, dim, heap);
            CSWAP_gate(control_qubit_index, n_counts + 2, n_counts + 3, state, dim, heap);
        } if (a == 4 || a == 11){
            CSWAP_gate(control_qubit_index, n_counts + 1, n_counts + 3, state, dim, heap);
            CSWAP_gate(control_qubit_index, n_counts + 0, n_counts + 2, state, dim, heap);
        } if (a == 7 || a == 11 || a == 13){
            for(int j = 0; j < 4; j++) {
                CNOT_gate(control_qubit_index, n_counts + j, state, dim);
            }
        }
    }
}

__kernel void shor(HEAP heap, COUNT __global* counts, size_t __global * n_counts) {
    int a = 7, N = 15, i;
    UINT dim = pow(2, *n_counts + 4);

    CTYPE __global *state = allocate_quantum_state(dim, heap);
    initialize_quantum_state(state, dim);

    for(i = 0; i < *n_counts; i++){
        H_gate(i, state, dim);
    }

    X_gate(*n_counts, state, dim);
    
    for(i = 0; i < *n_counts; i++){
        c_amod15(i, state, a, pow(2, i), *n_counts, dim, heap);
    }
    inverse_qft(0, *n_counts, 1, state, dim);

    UINT index_list[3] = {0, 1, 2};
    get_probabilities(index_list, 3, state, dim, counts);
}