#include "write_species.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern int cinnodes;
extern int cterminals;
extern int *cpy_stat_terminals;
extern int **terminals;

/*
void fixing_tree( Tree *t, double gap, double *sta_d, int *parent)
{  
    
	int i,j,temp;
    for( i = 0; i < cinnodes; i++)
    {
        if( !isnan( sta_d[i]))
        {
            if( sta_d[i] < gap)
            {
                for( j = 0; j < cterminals + cinnodes - 1 ; j++)
                {
                    if( (t->edge1[j] - cterminals - 1) == i)
                    {
                        temp = t->edge1[j];
                        t->edge1[j] = parent[ temp - 1];
                        
                        parent[ t->edge2[j] - 1] = parent[ temp - 1];
                    }
                    if( t->edge2[j] - cterminals - 1 == i)
                    {
                        t->edge1[j] = 0;
                        t->edge2[j] = 0;
                    }
                }
            }
        }
    }
}
*/

void write_species( char *seq_name[], double gap, double *sta_d, int *parent, char *output)
{ 
    int switch_terminals[cterminals+1];
    for( int i = 0; i < cterminals; i++)
    {
        switch_terminals[i] = 1;
    }

    FILE *f;
    f = fopen(output, "w");
    if( f == NULL)
    {
        printf("Can't create the species table file!\n");
        exit(1);
    }
    
	int i,j,k,temp_parent;
    for( i= 1; i < cterminals+1; i++)
    {
        if( switch_terminals[i] == 1)
        {
            j = i;
            PPP:
            temp_parent = parent[ j - 1];
            if( sta_d[ temp_parent - cterminals - 1]< gap)
            {
                j = temp_parent;
                goto PPP;
            }
            else
            {
                if( j<= cterminals)
                {
                    fprintf(f, "%s\n", seq_name[ j - 1]);
                }
                else
                {
                    for( k = 0; k < cpy_stat_terminals[ j - cterminals - 1] -1 ;k++)
                    {
                        fprintf( f,"%s,",seq_name[ terminals[ j - cterminals - 1][k] -1]);
                        switch_terminals[terminals[ j - cterminals - 1][k]] = 0;
                    }
                    fprintf(f, "%s\n",seq_name[ terminals[ j - cterminals - 1][cpy_stat_terminals[ j - cterminals - 1] - 1] -1]);
                    switch_terminals[terminals[ j - cterminals - 1][cpy_stat_terminals[ j - cterminals - 1] - 1] ] = 0;
                }   
            }
        }
    }
    fclose(f);
}