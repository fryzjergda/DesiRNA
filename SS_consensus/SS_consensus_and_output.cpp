// created on 06.02.2018 by Michal Boniecki
// program inputs set of secondary structures in dot and bracket format (in one file)
// and calculates base-pairing consensus

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"SS_record_and_matrix.h"
#include"SS_consensus_and_output.h"

#define MAX_SS_OUTPUT_LINES 16 // this is actually numner of psknots, levels of nesting, 16 is quite much

#define EPSILON 0.0000000001

SS_consensus_and_output::SS_consensus_and_output(){

    seq_str=NULL;
    seq_name=NULL;
    seq_length=0;

    consensus_matrix=NULL;
    dividers_matrix=NULL;
}

SS_consensus_and_output::~SS_consensus_and_output(){
/*
    if(consensus_matrix!=NULL){
	for(int i=0; i<seq_length; i++)
	    delete [] consensus_matrix[i];
	delete [] consensus_matrix;
    }

    if(seq_str!=NULL) delete [] seq_str;
    if(seq_length) seq_length=0;
    if(seq_name!=NULL) delete [] seq_name;
*/
}

void SS_consensus_and_output::init(T_SS_record_and_matrix *given_item){

    if(seq_name!=NULL){
	fprintf(stderr, "error in: void SS_consensus_and_output::init(SS_lines_and_matrix *given_item): seq_name is already initialized\n");
	exit(EXIT_FAILURE);
    }

    int curr_strlen;

//    curr_strlen = strlen(given_item->seq_name);
//    seq_name=new char[curr_strlen+1];
//    strcpy(seq_name, given_item->seq_name);

    if(seq_str!=NULL || seq_length!=0){
	fprintf(stderr, "error in: void SS_consensus_and_output::init(SS_lines_and_matrix *given_item): seq_str is already initialized\n");
	exit(EXIT_FAILURE);
    }

    seq_length = strlen(given_item->v_line);
//    seq_str = new char[seq_length+1];
//    strcpy(seq_str, given_item->seq_str);
}

void SS_consensus_and_output::create_consensus_matrix(T_SS_record_and_matrix *v_ss_records, int n_ss_records){

    if(consensus_matrix){
	fprintf(stderr, "improper call of void SS_consensus_and_output::create_consensus_matrix(SS_lines_and_matrix *v_ss_records, int n_ss_records), consensus_matrix seems to be alredy initialized\n");
	exit(EXIT_FAILURE);
    }

    int i,j,k;

    consensus_matrix = new float*[seq_length];
    for(j=0; j<seq_length; j++){
	consensus_matrix[j] = new float[seq_length];
    }
    for(j=0; j<seq_length; j++)
	for(i=0; i<seq_length; i++)
	    consensus_matrix[j][i] = 0;

    dividers_matrix = new float*[seq_length];
    for(j=0; j<seq_length; j++){
	dividers_matrix[j] = new float[seq_length];
    }
    for(j=0; j<seq_length; j++)
	for(i=0; i<seq_length; i++)
	    dividers_matrix[j][i] = 0;

    for(k=0; k<n_ss_records; k++){

	vector<T_pair> v_pairs = v_ss_records[k].get_pairs();
	vector<T_pair>::iterator it;

	for(it=v_pairs.begin(); it!=v_pairs.end(); ++it){

	    consensus_matrix[it->nucl_index_1][it->nucl_index_2] += v_ss_records[k].weight;
	    consensus_matrix[it->nucl_index_2][it->nucl_index_1] += v_ss_records[k].weight;
	}
/* --- old code by summing up contant matrix --- instead, now is reading just vectors of pairs
	for(j=0; j<seq_length; j++)
	    for(i=0; i<seq_length; i++)
//		if(v_ss_records[k].matrix_of_contacts[j][i] != 0)
		    consensus_matrix[j][i] += v_ss_records[k].matrix_of_contacts[j][i];
*/
    }

    for(j=0; j<seq_length; j++)
	for(i=0; i<seq_length; i++){

	    double divider=0;

	    for(k=0; k<n_ss_records; k++){

		char *v_line = v_ss_records[k].v_line;

		if(v_line[j]!='-' && v_line[j]!='~' && v_line[j]!='?' && v_line[j]!='_') // - ~ ? _ indicate lack of data, so this position (in matrix) is ommited for this record
		    if(v_line[i]!='-' && v_line[i]!='~' && v_line[i]!='?' && v_line[i]!='_') // - ~ ? _ indicate lack of data, so this position (in matrix) is ommited for this record
			dividers_matrix[j][i] += v_ss_records[k].weight;

	    }

	}

    for(j=0; j<seq_length; j++)
	for(i=0; i<seq_length; i++){

	    if(dividers_matrix[j][i] > EPSILON) // in fact should be 0.000 here, but am adding small epsilon, avoiding division by zero
		consensus_matrix[j][i] /= dividers_matrix[j][i];
//	    else
//		consensus_matrix[j][i] = 0.0;
	    else
		consensus_matrix[j][i] = -1; // indicated that there is no data at this position
	}

#ifdef DISPLAY_MATRICES
    for(j=0; j<seq_length; j++){
	for(i=0; i<seq_length; i++){
	    fprintf(stderr, "%d", int(9*consensus_matrix[j][i]));
	}
	fprintf(stderr, "\n");
    }
#endif
}

int SS_consensus_and_output::consensus_output(double cons_thrs){
//fprintf(stderr, "inside: int SS_consensus_and_output::consensus_output()\n");
    if(consensus_matrix==NULL){
	fprintf(stderr, "improper call of void SS_consensus_and_output::consensus_output(), consensus_matrix is not created\n");
	exit(EXIT_FAILURE);
    }
//fprintf(stderr, "after consensus_matrix NULL check\n");
    int *v_ss_lines[MAX_SS_OUTPUT_LINES]; //<--- in case of psedoknots next line will be allocated, that's why is a table of tables here
    int n_ss_lines;

    int *v_ss_line_non_specific;

    int i, i2, j, j2, k;
    int open_bracket_index, close_bracket_index;
    int partial_sum;

    float row_sum, row_max_elem_value;
    int row_max_elem_index, row_n_max_elems;

    float col_sum, col_max_elem_value;
    int col_max_elem_index, col_n_max_elems;

    float *curr_column = new float[seq_length];

    int *v_consensus_values = new int[seq_length];
    int *v_consensus_values_non_specific = new int[seq_length]; // this is just any pair on a given position
    int *v_empty_rows = new int[seq_length];

    int *v_occupied = new int[seq_length];
//    int *v_

    for(j=0; j<MAX_SS_OUTPUT_LINES; j++){
	v_ss_lines[j] = new int[seq_length];
	for(i=0; i<seq_length; i++)
	    v_ss_lines[j][i] = 0;
    }
    n_ss_lines=1;
//fprintf(stderr, "after zeroing v_ss_lines[][]\n");
    v_ss_line_non_specific = new int[seq_length];
    for(i=0; i<seq_length; i++)
	v_ss_line_non_specific[i] = 0;

    for(i=0; i<seq_length; i++){
	v_consensus_values[i] = 0;
	v_consensus_values_non_specific[i] = 0;
	v_occupied[i] = 0;
    }

    for(j=0; j<seq_length; j++){

	double curr_sum = 0;
	for(i=0; i<seq_length; i++)
	    curr_sum += consensus_matrix[j][i];

	if(curr_sum <= seq_length*(-1)+EPSILON && curr_sum >=seq_length*(-1)-EPSILON) // (-1) indicates no data
	    v_empty_rows[j]=1;
	else
	    v_empty_rows[j]=0;
//fprintf(stderr, "%d", v_empty_rows[j]);
    }
//fprintf(stderr, "\n");

    if(cons_thrs>0.0){
fprintf(stderr, "	applying threshold: %lf\n", cons_thrs);
	for(j=0; j<seq_length; j++)
	    for(i=0; i<seq_length; i++)
		if(consensus_matrix[j][i]<cons_thrs && consensus_matrix[j][i]>-1)
		    consensus_matrix[j][i] = 0;
fprintf(stderr, "	[DONE]\n");
    }


    for(j=0; j<seq_length; j++){
//fprintf(stderr, "main_loop: j: %3d\n",j);
	row_sum = calc_row_sum(consensus_matrix[j], seq_length);
	row_max_elem_index = find_max_row_elem_index(consensus_matrix[j], seq_length);
//fprintf(stderr, "1: row_max_elem_index: %d\n", row_max_elem_index);
	row_max_elem_value = consensus_matrix[j][row_max_elem_index];
//fprintf(stderr, "2: row_max_elem_index: %d\n", row_max_elem_index);
	row_n_max_elems = count_max_elems(consensus_matrix[j], row_max_elem_value, seq_length);
//fprintf(stderr, "m\n");
//fprintf(stderr, "row_sum: %f, max_elem_index: %2d, max_elem_value: %f, n_max_elems: %d\n", row_sum, max_elem_index, max_elem, n_max_elems);

//WHAT IF row_max_elem_value is ZERO
        if(row_n_max_elems >= 1){
//fprintf(stderr, "if\n");
	    get_column_as_row(curr_column, consensus_matrix, row_max_elem_index, seq_length);
	    col_sum = calc_row_sum(curr_column, seq_length);
	    col_max_elem_index = find_max_row_elem_index(curr_column, seq_length);
	    col_max_elem_value = consensus_matrix[col_max_elem_index][row_max_elem_index]; // the same as curr_column[col_max_elem_index]
	    col_n_max_elems = count_max_elems(curr_column, col_max_elem_value, seq_length);
//check_whether max_elem_value and max_elem_value_partner is the same
//		if(row_max_elem_value != col_max_elem_value){
//		    fprintf(stderr, "error in: int SS_consensus_and_output::consensus_output(float cons_thrs): row_max_elem_value!=col_max_elem_value --- %f!=%f\n", row_max_elem_value, col_max_elem_value);
//		    exit(EXIT_FAILURE);
//		}
//		n_max_elems_partner = count_max_elems(consensus_matrix[max_elem_index], max_elem_value, seq_length);
//fprintf(stderr, "m\n");
	    if(col_n_max_elems >= 1 && row_max_elem_value==col_max_elem_value){
//fprintf(stderr, "m\n");

		open_bracket_index = j;
		close_bracket_index = row_max_elem_index;
/*
		if(v_occupied[open_bracket_index] || v_occupied[close_bracket_index])
		    break;
		else{
		    v_occupied[open_bracket_index] = 1;
		    v_occupied[close_bracket_index] = 1;
		}
*/
		if(open_bracket_index < close_bracket_index)
		for(j2=0; j2<MAX_SS_OUTPUT_LINES; j2++){
//		    if(j2==0) fprintf(stderr, "m\n");
		    partial_sum = 0; //partial_sum indicates what's going on between opening and closing braket, during summation partial_sum cannot be positive value --- it means that there are closing brackets first
		    if(v_ss_lines[j2][open_bracket_index]==0 && v_ss_lines[j2][close_bracket_index]==0){
			for(k=open_bracket_index; k<=close_bracket_index; k++){ //also when partial sum is not equal 0 it means that there is a pseudoknot and program should go farther to the next lines
			    partial_sum += v_ss_lines[j2][k];
			    if(partial_sum != 0){
				partial_sum += 100000; //skipping
			        break;
			    }
			}
		    }
		    else partial_sum = -1;
//fprintf(stderr,"j = %d, partial_sum = %d\n",j,partial_sum);
			
		    if(partial_sum == 0){
			if(v_ss_line_non_specific[open_bracket_index] == 0 && v_ss_line_non_specific[close_bracket_index] == 0){

			    v_ss_lines[j2][open_bracket_index]= -1;
			    v_ss_lines[j2][close_bracket_index]= 1;
			    v_consensus_values[open_bracket_index] = int(9*row_max_elem_value);
			    v_consensus_values[close_bracket_index] = int(9*row_max_elem_value);

			    v_ss_line_non_specific[open_bracket_index] = 1;
			    v_ss_line_non_specific[close_bracket_index] = 1;
			    v_consensus_values_non_specific[open_bracket_index] = int(9*row_max_elem_value);
			    v_consensus_values_non_specific[close_bracket_index] = int(9*row_max_elem_value);
			}
			break;
		    }
		    
		    if(n_ss_lines==j2+1) n_ss_lines++;
//fprintf(stderr,"n_ss_lines = %d\n",n_ss_lines);
		}

		if(j2 == MAX_SS_OUTPUT_LINES){
		    fprintf(stderr,"error: int SS_consensus_and_output::consensus_output(): too many pseudoknots detected, increase MAX_SS_OUTPUT_LINES and recompile\n");
		    exit(EXIT_FAILURE);
		}
	    }
	    else{
		if(v_consensus_values[j]==0)
		    v_consensus_values[j] = int(9*(1.0 - 0.5*(row_sum+col_sum))); //here is a "dot consensus"  //just averave row_sum and col_sum

		if(v_consensus_values_non_specific[j]==0){
		    v_consensus_values_non_specific[j] = int(9*0.5*(row_sum+col_sum)); //just 1 - averave row_sum and col_sum
		    v_ss_line_non_specific[j] = 1;
		}
	    }
	}
	else{ //pairing is not unique so, but it should be indicated by "|" in the output
//fprintf(stderr, "else\n");
	    if(v_consensus_values[j]==0)
		v_consensus_values[j] = int(9*(1.0-row_sum)); //here is a "dot consensus"

	        if(v_consensus_values_non_specific[j]==0){
		    v_consensus_values_non_specific[j] = int(9*row_sum);
//		    v_ss_line_non_specific[j] = 1;
		}
	}
/*
	else{
	    int is_bracket_there = 0;
	    for(k=0; k<n_ss_lines; k++)
		if(v_ss_lines[k][j] != 0)
		    is_bracket_there = 1;
	    if(is_bracket_there == 0){
		v_consensus_values[j] = 9; //there is no pairing in this row so "dot consensus" is maximal
		v_consensus_values_non_specific[j] = 9;
	    }
	}
*/
//	else
//	    if(v_ss_lines[j]==0)
//		v_consensus_values[j] = int(9*(1.0-row_sum)); //here is a "dot consensus"
    }
//--------------------------------
    for(i2=0; i2<seq_length; i2++)
	if(v_consensus_values[i2]<0)
	    v_consensus_values[i2]=0;
    for(i2=0; i2<seq_length; i2++)
	if(v_consensus_values_non_specific[i2]>9)
	    v_consensus_values_non_specific[i2]=9;
//--------------------------------

    FILE *outfile=stdout;

    fprintf(outfile,"\n\n");
//    fprintf(outfile, ">consensus_pairs_unique multi line\n");
    for(i=0; i<n_ss_lines; i++){
	fprintf(outfile, ">consensus_pairs_unique_line_%02d\n",i);
	for(i2=0; i2<seq_length; i2++){
	    if(v_ss_lines[i][i2] == 0){
		if(v_empty_rows[i2] == 0)
		    fprintf(outfile,".");
		else
		    fprintf(outfile,"-");
//		if(is_chain_break(i2))
//		    fprintf(outfile," ");
	    }
	    else if(v_ss_lines[i][i2] == -1)
		fprintf(outfile,"(");
	    else if(v_ss_lines[i][i2] == 1)
		fprintf(outfile,")");
	    else{
		fprintf(stderr,"error processing data in int SS_consensus_and_output::consensus_output()\nv_ss_lines[][] can contain only values 0, -1 or 1, but in v_ss_lines[%d][%d] is %d\n",i,i2,v_ss_lines[i][i2]);
		exit(EXIT_FAILURE);
	    }
	}
	fprintf(outfile,"\n");
    }
    fprintf(outfile, ">consensus_pairs_unique_score\n");
    for(i2=0; i2<seq_length; i2++){
	if(v_empty_rows[i2]==0)
	    fprintf(outfile,"%1d",v_consensus_values[i2]);
	else
	    fprintf(outfile,"-");
    }
    fprintf(outfile,"\n\n");

//    fprintf(outfile, ">consensus_pairs_unique single line with pseudoknots\n");
    fprintf(outfile, ">consensus_pairs_unique_single_line\n");
    for(i2=0; i2<seq_length; i2++){
	int non_zero_found = 0; // this will inform whether in any line is something, it must be only one non zero element -> so it's psknot at a gicen level, or all elements may ne zero -> it's just a dot (no pairing)
	int which_level, opening_or_closing; //which_level will indicate level of psknotting (thus index of type of brackets), opening_or_closing will indicate whether bracket is opening (which is index 0 in v_bracket_pairs) or closing (which is index 1 in v_bracket_pairs)
	for(i=0; i<n_ss_lines; i++){

	    if(v_ss_lines[i][i2] != 0){
		if(non_zero_found){ //it means non zero is already found before, this is error, because it must one only one non zero element in column (pairing bracket) or non (just a dot)
		    fprintf(stderr, "error in processing pseudoknots for output in: int SS_consensus_and_output::consensus_output(), at array position %d there are at least two brackets, in fact should be 1 or 0\n", i2);
		    exit(EXIT_FAILURE);
		}
		else{
		    non_zero_found = 1;
		    which_level = i;
		    opening_or_closing = v_ss_lines[i][i2]; // here opening is -1 whereas closing is 1, later during printing -1 will become index 0, 1 will become 1 (so no change :-)
		}
	    }
	}
	if(non_zero_found == 0){
	    if(v_empty_rows[i2]==0)
		fprintf(outfile,".");
	    else
		fprintf(outfile,"-");
	}
	else{
	    if(opening_or_closing == -1)
		fprintf(outfile, "%c", v_bracket_pairs[which_level][0]);
	    else if(opening_or_closing == 1)
		fprintf(outfile, "%c", v_bracket_pairs[which_level][1]);
	    else{
		fprintf(stderr,"error processing data in int SS_consensus_and_output::consensus_output()\nv_ss_lines[][] can contain only values 0, -1 or 1, but in v_ss_lines[%d][%d] is %d\n",i,i2,v_ss_lines[i][i2]);
		exit(EXIT_FAILURE);
	    }
	}
    }
    fprintf(outfile,"\n");

    fprintf(outfile, ">consensus_pairs_unique_score\n");
    for(i2=0; i2<seq_length; i2++){
	if(v_empty_rows[i2]==0)
	    fprintf(outfile,"%1d",v_consensus_values[i2]);
	else
	    fprintf(outfile,"-");
    }
    fprintf(outfile,"\n");

    fprintf(outfile,"\n");
//    fprintf(outfile, ">consensus_flat\n");
    fprintf(outfile, ">consensus_flat\n");
    for(i2=0; i2<seq_length; i2++)
	if(v_ss_line_non_specific[i2] == 0){
	    if(v_empty_rows[i2]==0)
		fprintf(outfile, ".");
	    else
		fprintf(outfile, "-");
	}
	else
	    fprintf(outfile, "|");
    fprintf(outfile,"\n");

    fprintf(outfile, ">consensus_flat_score\n");
    for(i2=0; i2<seq_length; i2++){
	if(v_empty_rows[i2]==0)
	    fprintf(outfile,"%1d",v_consensus_values_non_specific[i2]);
	else
	    fprintf(outfile,"-");
    }
    fprintf(outfile,"\n\n");

    delete [] curr_column;

    delete [] v_consensus_values;
    for(j=0; j<MAX_SS_OUTPUT_LINES; j++)
	delete [] v_ss_lines[j];

    delete [] v_consensus_values_non_specific;
    delete [] v_ss_line_non_specific;
    delete [] v_occupied;
    delete [] v_empty_rows;
}

float SS_consensus_and_output::calc_row_sum(float *v_row, int n_elems){
//fprintf(stderr, "inside: calc_row_sum()\n");
    float curr_sum=0;

    for(int i=0; i<n_elems; i++)
	if(v_row[i]>=0) //excluding -1 values, which are indicates no data
	    curr_sum += v_row[i];

    return curr_sum;
}

int SS_consensus_and_output::find_max_row_elem_index(float *v_row, int n_elems){ //when there are more then one max elem the function returns index of the first one
//fprintf(stderr, "inside: find_max_row_elem_index()\n");
    float curr_max=-1e32;
    int max_index=0; //just in case, not to return random value

    for(int i=0; i<n_elems; i++)
	if(v_row[i] > curr_max){
	    curr_max = v_row[i];
	    max_index = i;
	}
    return max_index;
}

int SS_consensus_and_output::count_max_elems(float *v_row, float max_elem, int n_elems){
//fprintf(stderr, "inside: count_max_elems()\n");
    int occur_counter=0;

    for(int i=0; i<n_elems; i++)
	if(v_row[i] == max_elem)
	    occur_counter++;

    return occur_counter;
}

void SS_consensus_and_output::get_column_as_row(float *v_result, float **v2_matrix, int row_index, int n_elems){
//fprintf(stderr, "inside: get_column_as_row(): n_elems: %d\n", n_elems);
    for(int i=0; i<n_elems; i++)
	v_result[i] = v2_matrix[i][row_index];
}
