#include "fixing_tree.h"
#include <math.h>
#include <stdio.h>

extern int cinnodes;
extern int cterminals;

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