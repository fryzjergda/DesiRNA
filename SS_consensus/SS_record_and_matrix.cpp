// created on 09.03.2019 by Michal Boniecki
// program inputs set of secondary structures in dot and bracket format (in one file)
// and calculates base-pairing consensus

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "SS_record_and_matrix.h"

//#define DISPLAY_MATRICES

int   n_bracket_pairs=0; //value will be established later, by reading v_bracket_paris
const char *v_bracket_pairs[] = {"()", "[]", "{}", "<>", "Aa", "Bb", "Cc", "Dd", "Ee", "Ff", "Gg", "Hh", "Ii", "Jj", "Kk", "Ll", "Mm", "Nn", "Oo", "Pp", "Qq", "Rr", "Vv", "Ww", "Xx", "Yy", "Zz", NULL}; //if you want to service additional pair of brackets, just modify this list, NULL must be at the end
const int n_ascii_codes=256;
int v_bracket_to_index[n_ascii_codes];

T_SS_record_and_matrix::T_SS_record_and_matrix(){

    v_line = NULL;
    weight = 1.0;

    n_ss_lines = 1;

    matrix_of_contacts = NULL;
    matrix_of_dividers = NULL;
}


T_SS_record_and_matrix::~T_SS_record_and_matrix(){
/*
    int seq_length = strlen(v_line);

    if(matrix_of_contacts!=NULL){
	for(int i=0; i<seq_length; i++)
	    delete [] matrix_of_contacts[i];
	delete [] matrix_of_contacts;
    }

    if(matrix_of_dividers!=NULL){
	for(int i=0; i<seq_length; i++)
	    delete [] matrix_of_dividers[i];
	delete [] matrix_of_dividers;
    }

    if(v_line!=NULL) delete [] v_line;
*/
}

vector<T_pair> T_SS_record_and_matrix::calc_pairs(){

    if(v_pairs.size()!=0){
	fprintf(stderr, "vectors v_paris is not empty, contains: %d elements\n", (int)v_pairs.size());
	exit(EXIT_FAILURE);
    }

    int i,j,k;

    int seq_length = strlen(v_line);

    n_bracket_pairs=0;
    for(i=0; v_bracket_pairs[i]!=NULL; i++)
	if(strlen(v_bracket_pairs[i]) != 2){
	    fprintf(stderr, "incorrect definition of: v_bracket_pairs\n");
	    fprintf(stderr, "each element is a pair of brackets (two character strigs), entire array must be terminated by NULL\n");
	    fprintf(stderr, "example definition:\nconst char *v_bracket_pairs[] = {\"()\", \"[]\", \"{}\", \"<>\", NULL};\n");
	    exit(EXIT_FAILURE);
	}
	else
	    n_bracket_pairs++;

//    fprintf(stderr, "n_bracket_pairs: %d\n", n_bracket_pairs);

    for(i=0; i<n_ascii_codes; i++)
	v_bracket_to_index[i] = 0;

    for(i=0; i<n_bracket_pairs; i++){
	v_bracket_to_index[v_bracket_pairs[i][0]] = i+1; //indexing opening bracket
	v_bracket_to_index[v_bracket_pairs[i][1]] = -(i+1); //indexing closing bracket
    }

//    for(i=0; i<n_ascii_codes; i++)
//	fprintf(stderr, "%d %d\n", i, v_bracket_to_index[i]);
//finally array v_bracket_to_index is:
//v_bracket_to_index['('] = 1
//v_bracket_to_index[')'] = -1
//v_bracket_to_index['['] = 2
//v_bracket_to_index[']'] = -2
//...
//...
//and so on (for non declared brackets all values must be 0)

    int v_stack_tops[n_bracket_pairs];
    int *v_stacks[n_bracket_pairs];
    for(i=0; i<n_bracket_pairs; i++){
        v_stacks[i] = new int[seq_length];
        v_stack_tops[i] = 0;
    }

    for(k=0; k<n_ss_lines; k++){

	for(i=0; i<seq_length; i++){

//	    unsigned char curr_char = (unsigned char) v_ss_lines[k][i];
	    unsigned char curr_char = (unsigned char) v_line[i];

	    if(v_bracket_to_index[curr_char] == 0 && curr_char != '.' && curr_char != ',' && curr_char != 'x' && curr_char != '-' && curr_char != '~' && curr_char != '?' && curr_char != '_'){
//fprintf(stderr, "%d %d\n", curr_char, v_bracket_pairs[curr_char]);
		fprintf(stderr, "illegal character: %c in input string:\n%s\n", v_line[i], v_line);
		fprintf(stderr, "the only allowed symbols are:\n.\n&\n");
		for(j=0; j< n_bracket_pairs; j++)
		    fprintf(stderr, "%c %c\n", v_bracket_pairs[j][0], v_bracket_pairs[j][1]);
		exit(EXIT_FAILURE);
	    }

	    int elem_index = v_bracket_to_index[curr_char];
	    int array_index = abs(elem_index)-1;
	    int i_opening;

	    if(elem_index > 0){ //positive index means opening bracket, abs(elem_index)-1 is index of the stack
		v_stacks[array_index][v_stack_tops[array_index]] = i;
		v_stack_tops[array_index]++;
	    }
	    else if(elem_index < 0){ //negative index means closing bracket, abs(elem_index)-1 is index of the stack
		v_stack_tops[array_index]--;
		if(v_stack_tops[array_index]<0){
		    fprintf(stderr, "Too many closing brackets of type: %c in:\n%s\n", v_bracket_pairs[array_index][1], v_line);
		    exit(EXIT_FAILURE);
		}
		i_opening = v_stacks[array_index][v_stack_tops[array_index]];

		T_pair curr_pair;
		curr_pair.nucl_index_1 = i_opening;
		curr_pair.nucl_index_2 = i;

		v_pairs.push_back(curr_pair);

//		matrix_of_contacts[i_opening][i] = weight; //here i is i_closing //old code --- loading to matrix
//		matrix_of_contacts[i][i_opening] = weight; //and again to make it symmetric --- loading to matrix
//	    fprintf(stderr, "%d, %d\n", i_opening, i);
	    }
	}

	for(i=0; i<n_bracket_pairs; i++)
	    if(v_stack_tops[i] > 0){
		fprintf(stderr, "Too many opening brackets of type: %c in:\n%s\n", v_bracket_pairs[i][0], v_line);
		exit(EXIT_FAILURE);
	}
    }

    for(i=0; i<n_bracket_pairs; i++)
        delete [] v_stacks[i];

    return v_pairs;
}

vector<T_pair> T_SS_record_and_matrix::get_pairs(){

    return v_pairs;
}


void T_SS_record_and_matrix::create_matrix_of_contacts(){

    if(matrix_of_contacts){
	fprintf(stderr, "improper call of void T_SS_record_and_matrix::create_matrix_of_contacts(), matrix_of_contacts seems to be alredy initialized\n");
	exit(EXIT_FAILURE);
    }

    int i,j,k;

    int seq_length = strlen(v_line);

    matrix_of_contacts = new float*[seq_length];
    for(j=0; j<seq_length; j++)
	matrix_of_contacts[j] = new float[seq_length];

    for(j=0; j<seq_length; j++)
	for(i=0; i<seq_length; i++)
	    matrix_of_contacts[j][i] = 0.0;

    calc_pairs();

    vector<T_pair>::iterator it;

    for(it=v_pairs.begin(); it!=v_pairs.end(); ++it){
//fprintf(stderr, "pair_to_matrix: %3d %3d\n", it->nucl_index_1, it->nucl_index_2);
	matrix_of_contacts[it->nucl_index_1][it->nucl_index_2] = weight;
	matrix_of_contacts[it->nucl_index_2][it->nucl_index_1] = weight;
    }
fprintf(stderr, "\n");

//    fprintf(stderr, "ss_str: %s\n", ss_str);
#ifdef DISPLAY_MATRICES
    for(j=0; j<seq_length; j++){
	for(i=0; i<seq_length; i++){
	    if(matrix_of_contacts[j][i] == 0)
		fprintf(stderr, ".");
	    else
		fprintf(stderr, "x");
	}
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
#endif
}

void T_SS_record_and_matrix::create_matrix_of_dividers(){

    if(matrix_of_dividers){
	fprintf(stderr, "improper call of void T_SS_record_and_matrix::create_matrix_of_dividers(), matrix_of_dividers seems to be alredy initialized\n");
	exit(EXIT_FAILURE);
    }

    int seq_length = strlen(v_line);

    int i,j;

    matrix_of_dividers = new float*[seq_length];
    for(j=0; j<seq_length; j++)
	matrix_of_dividers[j] = new float[seq_length];

    for(j=0; j<seq_length; j++)
	for(i=0; i<seq_length; i++)
	    matrix_of_dividers[j][i] = 0.0;

    for(j=0; j<seq_length; j++)
	if(v_line[j]!='-' && v_line[j]!='~' && v_line[j]!='?' && v_line[j]!='_')
	    for(i=0; i<seq_length; i++)
		if(v_line[i]!='-' && v_line[i]!='~' && v_line[i]!='?' && v_line[i]!='_')
			    matrix_of_dividers[j][i] = weight;

#ifdef DISPLAY_MATRICES
    for(j=0; j<seq_length; j++){
	for(i=0; i<seq_length; i++){
	    if(matrix_of_dividers[j][i] == 0)
		fprintf(stderr, ".");
	    else
		fprintf(stderr, "x");
	}
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
#endif

}
