__kernel void qpe(HEAP heap, COUNT __global* counts, size_t __global * n_counts) {
    int dim = pow(2, *n_counts + 1);
    CTYPE __global *state = allocate_quantum_state(dim, heap);
    initialize_quantum_state(state, dim);
    int i, j, exp;

    double phase_input = 1./2;

    CTYPE matrix[4] = {{1, 0}, {0, 0}, {0, 0}, {0, 0}};

    for(i = 0; i<3; i++) {
        H_gate(i, state, dim);
    }
    X_gate(3, state, dim);

    for(i = 0; i<3; i++) {
        exp = 1 << (3-i-1);
        matrix[3] = double_to_complex(cos(2*PI*phase_input*exp), sin(2*PI*phase_input*exp));
        single_qubit_control_single_qubit_dense_matrix_gate(
            i, 1, 3, matrix, state, dim);
    }
    inverse_qft(0, 3, 0, state, dim);

    UINT index_list[3] = {0, 1, 2};

    get_probabilities(index_list, 3, state, dim, counts);
}