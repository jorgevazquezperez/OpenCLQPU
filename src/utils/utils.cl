UINT pow(int x, int n)
{
    // Initialize result to 1
    UINT pow = 1;
 
    // Multiply x for n times
    for (int i = 0; i < n; i++) {
        pow = pow * x;
    }
 
    return pow;
}

void get_binary(UINT value, UINT __global* qubit_value){
    int i;

    for (i = 0; i <= 3; i++){
        qubit_value[i] = 0;
    }

    for(i = 0; value > 0; i++){    
        qubit_value[i] = value % 2;    
        value = value / 2; 
    }   
}

typedef struct{
    UINT qubit_value[3];
    double prob;
} COUNT;

COUNT __global* get_probabilities(UINT* sorted_target_qubit_index_list, UINT target_qubit_index_count,
     CTYPE __global* state, ITYPE dim, COUNT __global* counts){

        int i, j;
        UINT sub_dim = pow(2, target_qubit_index_count);

        for(i = 0; i < sub_dim; i++){
            get_binary(i, counts[i].qubit_value);
            
            counts[i].prob = marginal_prob(sorted_target_qubit_index_list, counts[i].qubit_value, 
                target_qubit_index_count, state, dim);
        }
        return counts;
     }

void print_counts(COUNT __global* counts, UINT target_qubit_index_count){
    
    printf("Probabilities:\n");
    for(int i = 0; i < pow(2, target_qubit_index_count); i++){
        for(int j = target_qubit_index_count - 1; j >= 0; j--){
            printf(" %d", counts[i].qubit_value[j]);
        }
        printf(": %f\n", counts[i].prob);
    }

}

void print_state(CTYPE __global *state, UINT dim){
  for(int i=0; i<dim; i++){
        if(i % 3 == 0){
          printf("\n");
        }
        printComplex(state[i]);
        printf("\t");
    }
    printf("\n");
}