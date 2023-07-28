// created on 06.02.2018 by Michal Boniecki, developed on 03.09.2019
// program inputs set of secondary structures in dot and bracket format (in one file)
// and calculates base-pairing consensus

#ifndef _SS_CONSENSUS_AND_OUTPUT_
#define _SS_CONSENSUS_AND_OUTPUT_

class SS_consensus_and_output{

private:

    float calc_row_sum(float *v_row, int n_elems);
    int find_max_row_elem_index(float *v_row, int n_elems);
    int count_max_elems(float *v_row, float max_elem, int n_elems);

    void get_column_as_row(float *v_result, float **v2_matrix, int row_index, int n_elems);
//    float *get_column_as_row(float *v_result, float **v2_matrix, int row_index, int n_elems);

public:

    char *seq_str;
    char *seq_name;
    int seq_length;

    float **consensus_matrix;
    float **dividers_matrix;

    void init(T_SS_record_and_matrix *given_item);
    void create_consensus_matrix(T_SS_record_and_matrix *v_ss_records, int n_ss_records);
    int consensus_output(double cons_thrs);

    SS_consensus_and_output();
    ~SS_consensus_and_output();
};

#endif
