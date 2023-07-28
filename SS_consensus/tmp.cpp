
int find_max_row_elem_index(float *v_row, int n_elems){ // incorrect function

    float curr_max=0;
    int max_index;

    for(int i=0; i<n_elems; i++)
	if(v_row[i] > curr_max){
	    curr_max = v_row[i];
	    max_index = i;
	}

    return max_index;
}

