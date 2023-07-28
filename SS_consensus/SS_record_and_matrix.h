// created on 09.03.2019 by Michal Boniecki
// program inputs set of secondary structures in dot and bracket format (in one file)
// and calculates base-pairing consensus

#ifndef _SS_RECORD_AND_MATRIX_
#define _SS_RECORD_AND_MATRIX_

#include <vector>

using namespace std;

extern const char *v_bracket_pairs[];

typedef struct{
    int nucl_index_1, nucl_index_2;
} T_pair;

class T_SS_record_and_matrix{

public:

    char *v_line;
    double weight;

    int n_ss_lines;

    vector<T_pair> v_pairs;

    vector<T_pair> calc_pairs();
    vector<T_pair> get_pairs();

    float **matrix_of_contacts;
    float **matrix_of_dividers;

    void create_matrix_of_contacts();
    void create_matrix_of_dividers();

    T_SS_record_and_matrix();
    ~T_SS_record_and_matrix();
};

#endif //_SS_RECORD_AND_MATRIX_
