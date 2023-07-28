
#include "functions.h"

//#define MAX_N_ITEMS 16

int split_string(char *string2parse, const char *delim_string, char **output_list, int max_n_items)
{
    char *curr_token,*saveptr;
    int n_items=0;

    curr_token=strtok_r(string2parse,delim_string,&saveptr);
    while(curr_token!=NULL)
    {
	output_list[n_items]=curr_token;
	n_items++;
	curr_token=strtok_r(NULL,delim_string,&saveptr);
	if(n_items>=MAX_N_ITEMS)
	{
	    fprintf(stderr,"too many items in output list from parsing in function split_string(char*,const char*,char**,int)\n");
	    fprintf(stderr,"increase constant MAX_N_ITEMS and recompile\n");
	}
    }
    return n_items;
}
