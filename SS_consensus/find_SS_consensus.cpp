// created on 24.01.2018 by Michal Boniecki, further development on 03.09.2019
//- implementation of new idea - allowing for gapped SS lines and allowing for weights declared for specific lines

// program inputs set of secondary structures in dot and bracket format (in one file) (they may contain gaps, and optional weight at the end of each line)
// and calculates base-pairing consensus

//comparing to the initial version (from January of 2018) pretty much is changed.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions.h"
#include "SS_record_and_matrix.h"
#include "SS_consensus_and_output.h"

#define MAX_N_SS_RECORDS 1000000

//int   n_bracket_pairs=0; //value will be established later, by reading v_bracket_paris
//const char *v_bracket_pairs[] = {"()", "[]", "{}", "<>", NULL}; //if you want to service additional pair of brackets, just modify this list, NULL must be at the end
//const int n_ascii_codes=256;
//int v_bracket_to_index[n_ascii_codes];

/*
typedef struct{
    int begin_index, n_lines;
} T_record_index;
*

/*
typedef struct{

    char *v_line;
    double weight;

} T_SS_record;
*/
/*
class T_SS_record{

public:

    char *v_line;
    double weight;
};
*/


int main(int argc, char *argv[]){


    if(argc < 2){
        fprintf(stderr, "usage: find_SS_consensus file_with_SS_dot-bracket_lines.txt [consensus_thrs]\n");
	fprintf(stderr, "see attached example input files !\n");
	exit(EXIT_FAILURE);
    }

    FILE *fp;
    char *filename_str = argv[1];

    double consensus_thrs=0.0;

    if(argc > 2){
	if(sscanf(argv[2], "%lf", &consensus_thrs) != 1){
	    fprintf(stderr, "error converting: %s (second parameter, which should be consensus_thrs value) into float\n", argv[2]);
	    exit(EXIT_FAILURE);
	}
    }

fprintf(stderr, "reading input\n");
    if((fp=fopen(filename_str, "r"))==NULL){
	fprintf(stderr, "error opening file: %s for reading\n", filename_str);
	exit(EXIT_FAILURE);
    }

    char *curr_line = NULL;
    size_t line_len = 0;
    size_t ss_len = 0;
    int n_line_cnt = 1;

    while((getline(&curr_line, &line_len, fp)) != -1) // pre-read to establish how many lines is in the file, some lines may be empty, but I don't care at this moment, count all
	n_line_cnt++;
fprintf(stderr, "file: %s contains: %d lines\n", filename_str, n_line_cnt);

    T_SS_record_and_matrix *v_SS_records = new T_SS_record_and_matrix[n_line_cnt];
    int n_SS_records = 0;

    char *splitted_line[MAX_N_ITEMS];

    rewind(fp);
    n_line_cnt = 1;

    while((getline(&curr_line, &line_len, fp)) != -1){

//	fprintf(stderr, "1_%s_\n", curr_line);
//fprintf(stderr, "%d\n", strlen(curr_line));
	if(strlen(curr_line) == 0){
	    n_line_cnt++;
	    continue;
	}
//	fprintf(stderr, "2_%s_\n", curr_line);

	if(curr_line[strlen(curr_line)-1] == '\n')
	    curr_line[strlen(curr_line)-1] = 0; //deleting new line at the end of the string
	if(curr_line[strlen(curr_line)-1] == 0xD)
	    curr_line[strlen(curr_line)-1] = 0; //just in case, if DOS new line is provided

	int n_items = split_string(curr_line, " ", splitted_line, MAX_N_ITEMS);
//fprintf(stderr, "n_items: %d\n", n_items);

	if(n_items > 0){
	    if(ss_len == 0){ //it means this is first line, establishing ss length
		ss_len = strlen(splitted_line[0]);
//fprintf(stderr, "first ss_len: %d\n", (int)ss_len);
	    }
	    else{ //not first line, ss length is already established
		int curr_ss_len = strlen(splitted_line[0]);
		if(curr_ss_len != ss_len){
		    fprintf(stderr, "incorrect SS line length\n");
		    fprintf(stderr, "size of previous SS lines is: %d, but line nr: %d is: %d chars long\n", (int)ss_len, n_line_cnt, (int)curr_ss_len);
		    exit(EXIT_FAILURE);
		}
	    }

	    v_SS_records[n_SS_records].v_line = new char[ss_len+1];
	    strcpy(v_SS_records[n_SS_records].v_line, splitted_line[0]);
	    v_SS_records[n_SS_records].weight = 1.0; //defaut value

//	    fprintf(stderr, "line nr: %d\n", n_line_cnt);
//	    fprintf(stderr, "%s\n", v_SS_records[n_SS_records].v_line);

	    if(n_items > 1){ //second item will be interpreted as weight
		if(sscanf(splitted_line[1], "%lf", &v_SS_records[n_SS_records].weight) != 1){
		    fprintf(stderr, "error converting %s into float\n", splitted_line[1]);
		    fprintf(stderr, "file: %s, line: %d\n", filename_str, n_line_cnt);
		    fprintf(stderr, "second item in line, just after SS line is interpreted as weight, should vary in range (0.0:1.0>\n");
		    fprintf(stderr, "default value is: 1.0\n");
		    exit(EXIT_FAILURE);
		}
		if(v_SS_records[n_SS_records].weight > 1.0){
		    fprintf(stderr, "declared weight of SS line nr: %d is: %lf\nthe the value is set to 1.0, default value\n", n_line_cnt, v_SS_records[n_SS_records].weight);
		    fprintf(stderr, "second item in line, just after SS line is interpreted as weight, should vary in range (0.0:1.0>\n");
		    v_SS_records[n_SS_records].weight = 1.0; //defaut value
		}
//		fprintf(stderr, "line nr: %d\n", n_line_cnt);
//		fprintf(stderr, "%s\n", v_SS_records[n_SS_records].v_line);
		if(v_SS_records[n_SS_records].weight > 0.0 && v_SS_records[n_SS_records].weight <= 1.0)
		fprintf(stderr, "is loaded with weight: %lf\n\n", v_SS_records[n_SS_records].weight);

		if(v_SS_records[n_SS_records].weight <= 0.0){
		    fprintf(stderr, "declared weight of SS line nr: %d is: %lf\nthe line is skipped\n", n_line_cnt, v_SS_records[n_SS_records].weight);
		    fprintf(stderr, "second item in line, just after SS line is interpreted as weight, should vary in range (0.0:1.0>\n\n");
		    delete [] v_SS_records[n_SS_records].v_line;
		    n_SS_records--;
		}
		
	    }

	    n_SS_records++;
	
//	    if(n_SS_records>=MAX_N_SS_RECORDS){
//		fprintf(stderr, "too many SS_records, increase MAX_N_SS_RECORDS and recompile\n");
//		exit(EXIT_FAILURE);
//	    }
	}
	n_line_cnt++;
    }

    fclose(fp);
fprintf(stderr, "[DONE]\n");

    if(curr_line)
	free(curr_line);
fprintf(stderr, "deriving lists of pairs\n");
    for(int i=0; i<n_SS_records; i++){
//	fprintf(stderr, "record: %d create contact_matrix\n", i+1);
//	v_SS_records[i].create_matrix_of_contacts(); //--- old code, contact_matrix is created only as a consensus matrix, not for every SS case
	v_SS_records[i].calc_pairs(); //instread of contact matrix just list of pairs is calculated
    }
fprintf(stderr, "[DONE]\n");

//fprintf(stderr, "before matrix_of_diveders calculation\n"); //--- old code
//    for(int i=0; i<n_SS_records; i++){
//	fprintf(stderr, "record: %d create contact_matrix\n", i+1);
//	v_SS_records[i].create_matrix_of_dividers();
//    }


//fprintf(stderr, "before SS_consensus_and_output declaration\n");
    SS_consensus_and_output a_consensus_and_output;
//fprintf(stderr, "before .init(v_SS_records)\n");
    a_consensus_and_output.init(v_SS_records);
//fprintf(stderr, "before .create_consensus_matrix\n");
fprintf(stderr, "calculation of consensus matrix\n");
    a_consensus_and_output.create_consensus_matrix(v_SS_records, n_SS_records);
fprintf(stderr, "[DONE]\n");
//    fprintf(stdout, "----------------------- consensus results ----------------------------------------\n");
//fprintf(stderr, "before .consensus_output()\n");
fprintf(stderr, "outputing\n");
    a_consensus_and_output.consensus_output(consensus_thrs);
fprintf(stderr, "[DONE]\n");

    delete [] v_SS_records;

    return 0;
}
