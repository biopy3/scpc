#include "newick.h"

double get_gap( Tree *t, int *parent, double *sta_d, int **stat, int *bins, double *stat_max, double *stat_min);
void get_d_index( int tip1, int tip2, int num_seqs, int *index_edge);