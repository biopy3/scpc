#include <stdio.h>
#include <stdlib.h>
#include "statistic.h"
#include <string.h>
#include <math.h>
#include "root.h"

#define BIN_LEN 0.005
extern int cinnodes;
extern int cterminals;
extern double *d;

int index_min(double *p,int n)
{
    int i,max_i;
    for(i=0; i<n; i++)
    {
        if( isnormal(*(p+i)))
        {
            max_i = i;
            break;
        }
    }
    for(; i + 1<n; i++)
    {
        if( isnormal(*(p+i)))
        {
            if(*(p+max_i) > *(p+i))
            {    
                max_i=i;
            }
        }
    }
    return max_i;
}

void get_d_index( int tip1, int tip2, int num_seqs, int *index_edge)
{
    int temp;
    int init_num_seqs = num_seqs;
    *index_edge = -1;
    if( tip1 < tip2)
    {
        temp = tip1;
        tip1 = tip2;
        tip2 = temp;
    }else
    {
        temp = tip2;
    }
    
    while( 1 < tip2 )
    {
       (*index_edge) += (num_seqs - 1);
        num_seqs--;
        tip2--;
    }
    
    (*index_edge) += ( tip1 - temp);
}

void get_stat_terminals( Tree *t, int *stat_terminals, int innode,int edge2)
{
  for( int j = 0; j < cterminals + cinnodes - 1; j++)
    {
      if( t->edge1[j] == edge2)
      {
        if( t->edge2[j] <= cterminals)
        {
          stat_terminals[ innode]++;
        }
        else
        {
          get_stat_terminals( t, stat_terminals, innode, t->edge2[j]);
        }
      }
    }
}

void get_terminals( Tree *t, int *stat_terminals, int **terminals,int innode,int edge2)
{
  for( int j = 0; j < cterminals + cinnodes - 1; j++)
    {
      if( t->edge1[j] == edge2)
      {
        if( t->edge2[j] <= cterminals)
        {
          terminals[ innode][ stat_terminals[ innode ] - 1] = t->edge2[j];
          stat_terminals[ innode ]--;
        }
        else
        {
          get_terminals( t, stat_terminals, terminals, innode, t->edge2[j]);
        }
      }
    }
}

double get_gap( Tree *t, int *parent, double *sta_d, int **stat, int *bins, double *stat_max, double *stat_min)
{
  int stat_children[cinnodes];
  int stat_terminals[cinnodes];
  int cpy_stat_children[cinnodes];
  int cpy_stat_terminals[cinnodes];
  int temp[cinnodes];

  for( int i = 0; i < cinnodes;i++)
  {
    stat_children[i] = 0;
    stat_terminals[i] = 0;
    temp[i] = 0;
  }
  
  /* get statistic */
  /* for children */
  for( int i = 0; i < cterminals + cinnodes - 1; i++)
  {
    stat_children[ t->edge1[i] - cterminals - 1]++;
  }

  /* for terminals */
  for( int i = 0; i < cinnodes; i++)
  {
    for( int j = 0; j < cterminals + cinnodes - 1; j++)
    {
      if( t->edge1[j] - cterminals - 1 == i)
      {
        if( t->edge2[j] <= cterminals)
        {
          stat_terminals[ i]++;
        }
        else
        {
          get_stat_terminals( t, stat_terminals, i, t->edge2[j]);
        }
      } 
    }
  }

  /* copy stat_children and stat_terminals */
  memcpy( cpy_stat_children, stat_children, sizeof( stat_children));
  memcpy( cpy_stat_terminals, stat_terminals, sizeof( stat_terminals));
  
  int *children[cinnodes];
  int *terminals[cinnodes];
  /* allocate memory for children and terminals */
  for( int i = 0; i < cinnodes; i++)
  {
    children[i] = ( int *)malloc( stat_children[i]*sizeof(int));
    terminals[i] = ( int *)malloc( stat_terminals[i]*sizeof(int));
  }

  /* get children */
  for( int i = 0; i < cterminals + cinnodes - 1; i++)
  {
    children[ t->edge1[i] - cterminals - 1][ stat_children[ t->edge1[i] - cterminals - 1] - 1] = t->edge2[i];
    stat_children[ t->edge1[i] - cterminals - 1]--;
  }
  /* get parent */
  for( int i = 0; i < cterminals + cinnodes - 1; i++)
  {
      parent[ t->edge2[i] - 1] = t->edge1[i];
  }
  parent[ cterminals] = cterminals +1;

  /* get terminals */
  for( int i = 0; i < cinnodes; i++)
  {
    for( int j = 0; j < cterminals + cinnodes - 1; j++)
    {
      if( t->edge1[j] - cterminals - 1 == i)
      {
        if( t->edge2[j] <= cterminals)
        {
          terminals[ i][ stat_terminals[ i]- 1] = t->edge2[j];
           stat_terminals[ i]--;
        }
        else
        {
          get_terminals( t, stat_terminals, terminals,i, t->edge2[j]);
        }
      } 
    }
  }

  /* print the children and terminals */
  // for( int i = 0; i < cinnodes; i++)
  // {
  //   for( int j = 0; j < cpy_stat_children[i]; j++)
  //   {
  //     printf("children[%d]:%d\n", i, children[i][j]);
  //   }
  //   printf("\n");
  // } 

  // for( int i = 0; i < cinnodes; i++)
  // {
  //   for( int j = 0; j < cpy_stat_terminals[i]; j++)
  //   {
  //     printf("terminals[%d]:%d\n", i, terminals[i][j]);
  //   }
  //   printf("\n");
  // } 

    /* print the parents */
  // for( int i = 0; i < cinnodes + cterminals; i++)
  // {
  //   printf("parent[%d]:%d\n", i, parent[i]);
  // } 



  /* compute the distances of internal nodes */
  int count1,count2;
  double sum1 = 0.0;
  double sum2 = 0.0;
  int index_d;

  for( int i = 0; i < cinnodes; i++)
  {
    count1 = 0;
    sum1 = 0.0;
    for( int j = 0; j < cpy_stat_children[i]; j++)
    {
      for( int k = j + 1; k < cpy_stat_children[i]; k++)
      {
        count1++;
        sum2 = 0.0;

        if( children[i][j] <= cterminals && children[i][k] <= cterminals)
        {
            get_d_index( children[i][j], children[i][k], cterminals, &index_d);
            sum2 += d[index_d];
            sum1 += sum2;
        }
        else if( children[i][j] > cterminals && children[i][k] <= cterminals)
        {
            count2 = 0;
            for( int s = 0; s < cpy_stat_terminals[ children[i][j] - cterminals - 1]; s++)
            {
                get_d_index( terminals[ children[i][j] - cterminals - 1][s], children[i][k], cterminals, &index_d);
                sum2 += d[index_d];
                count2++;
            }
            sum1 += (sum2 / count2);
        }
        else if( children[i][j] <= cterminals && children[i][k] > cterminals)
        {
            count2 = 0;
            for( int q = 0; q < cpy_stat_terminals[ children[i][k] - cterminals - 1]; q++)
            {
                get_d_index( children[i][j], terminals[ children[i][k] - cterminals - 1][q], cterminals, &index_d);
                sum2 += d[index_d];
                count2++;
            }
            sum1 += (sum2 / count2);
        }
        else
        {
            count2 = 0;
            for( int s = 0; s < cpy_stat_terminals[ children[i][j] - cterminals - 1]; s++)
            {
                for( int q = 0; q < cpy_stat_terminals[ children[i][k] - cterminals - 1]; q++)
                {
                    get_d_index( terminals[ children[i][j] - cterminals - 1][s], terminals[ children[i][k] - cterminals - 1][q], cterminals, &index_d);
                    sum2 += d[index_d];
                    count2++;
                }
            }
            sum1 += (sum2 / count2);
        }
      }
    }
    sum1 /= count1;
    sta_d[i] = sum1;
  }

  for( int i = 0; i < cinnodes; i++)
  {
    printf("%lf ", sta_d[i]);
  }
  printf("\n");
  *stat_max = sta_d[ index_max(sta_d, cinnodes)];
  printf("stat_max:%lf\n", *stat_max);
  *stat_min = sta_d[ index_min(sta_d, cinnodes)];
  printf("stat_min:%lf\n", *stat_min);

  *bins = ceil( (*stat_max - *stat_min) / BIN_LEN);
  *stat = ( int *)calloc(*bins, sizeof(int));

  for( int i = 0; i < cinnodes ; i++)
  {
    for( int j = 1; j <= *bins; j++)
    {  
      if( sta_d[i] <= *stat_min + j*BIN_LEN)
      {
        (*stat)[j - 1]++;
        break;
      }    
    }    
  }

  // for( int i = 0; i < *bins; i++)
  // {
  //     printf("stat[%d]:%d\n", i, (*stat)[i]);
  // }

  int gap_index = (*stat)[0];
  for( int i = 0; i < *bins-1; i++)
  {
    if( (*stat)[i] == 0)
    {
      if( (*stat)[i+1] == 0) 
      {
        ;
      }
      else
      {
        gap_index = i;
        break;
      }
    }
  }

  // printf("gap_index:%d\n", gap_index);
  printf("gap:%lf\n", *stat_min + gap_index * BIN_LEN);

  /* free children and terminals */
  for( int i = 0; i < cinnodes; i++)
  {
    free(children[i]);
    free(terminals[i]);
  }
  return *stat_min + (gap_index + 1) * BIN_LEN;
}