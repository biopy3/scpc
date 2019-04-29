#include "scpc.h"
#include <math.h>

int index_max(double *p,int n)
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
            if(*(p+max_i) < *(p+i))
            {    
                max_i=i;
            }
        }
    }
    return max_i;
}

static int iii;

void bar_reorder(int node, int n, int m, int *e1, int *e2, int *neworder, int *L, int *pos)
{
	int i = node - n - 1, j, k;

	for (j = pos[i] - 1; j >= 0; j--)
		neworder[iii--] = L[i + m * j] + 1;

	for (j = 0; j < pos[i]; j++) {
		k = e2[L[i + m * j]];
		if (k > n)
			bar_reorder(k, n, m, e1, e2, neworder, L, pos);
	}
}

int *postorder_reorder_index(int *e1, int *e2, int num_seqs)
{
    /* n: nb of tips
    m: nb of nodes
    N: nb of edges */
    int i, j, k, *L, *pos, *neworder, N = 2*num_seqs - 3, n = num_seqs;
    int m = N - n + 1, degrmax = n - m + 1;

    printf("N:%d\n", N);
    printf("n:%d\n", n);
    printf("m:%d\n", m);
    L = (int*)malloc((m * degrmax)*sizeof(int));
    pos = (int*)malloc(m*sizeof(int));
    neworder = (int*)malloc(N*sizeof(int));
	memset(pos, 0, m * sizeof(int));

    for (i = 0; i < N; i++)
    {
		k = e1[i] - n - 1; /* k is the 'row' index in L corresponding to node e1[i] */
		j = pos[k]; /* the current 'column' position corresponding to k */
		pos[k]++; /* increment in case the same node is found in another row of the edge matrix */
		L[k + m * j] = i;
	}

    iii = N - 1;
    bar_reorder(n + 1, n, m, e1, e2, neworder, L, pos);
    free(pos);
    free(L);
    return neworder;
}

void postorder_reorder(int *e1, int *e2, double *edge_length,int *neworder, int len_edge)
{
    int temp1[len_edge], temp2[len_edge];
    double temp_edge_length[len_edge];
    for( int i = 0; i<len_edge; i++)
    {
        temp1[i] = e1[i];
        temp2[i] = e2[i];
        temp_edge_length[i] = edge_length[i];
    }
    for( int i = 0; i<len_edge; i++)
    {
        e1[i] = temp1[neworder[i] - 1];
        e2[i] = temp2[neworder[i] - 1];
        edge_length[i] = temp_edge_length[neworder[i] - 1];
    }
}

int *statis_nodes(int *edge1 ,int *edge2, int nodes)
{
    int *res;
    res = (int *)malloc(nodes*sizeof(int));
    /* init res to zero */
    for(int i = 0; i<nodes; i++)
        res[i] = 0;
    for(int i = 0; i<nodes; i++)
    {
        res[ edge1[i]-1 ]++;
        res[ edge2[i]-1 ]++;
    }
    return res;
}

void flounder_sub(int *tip1, int *tip2, int num_seqs, int index_edge)
{
    int init_num_seqs = num_seqs;
    *tip1 = 1, *tip2 = 0;
    num_seqs--;
    index_edge++;
    while(index_edge > num_seqs)
    {
        (*tip1)++;
        index_edge -= num_seqs;
        num_seqs--;
    }
    *tip2 = index_edge + (init_num_seqs - num_seqs);

}

int is_all_degree_than_1_except_tip1andtip2( int tip1, int tip2, int *stat_nodes , int nodes)
{
    for( int i = 0; i< nodes; i++)
    {
        if( i == tip1 - 1 || i == tip2 - 1 )
            continue;
        if(stat_nodes[i] == 1)
            //printf("%d\n", i);
            return 0;
    }
    return 1;
}

void f( int internal_node, int child, int *e1, int *e2,
            int *er1, int *er2, int num_seqs)
{
    extern int count;
    int nodes = 2*num_seqs -2, edges = nodes - 1;
    
    er1[count] = internal_node + 1;
    if( child <= num_seqs)
        er2[count++] = child;
    else
        er2[count++] = child + 1;
    
    if( child <= num_seqs)
    {
        ;
    }else
    {
        for( int i = 0; i < edges; i++)
        {
            if( e1[i] == child && e2[i] != internal_node)
            {
                f( child, e2[i], e1, e2, er1, er2, num_seqs);
            }else if( e2[i] == child && e1[i] != internal_node)
            {
                f( child, e1[i], e1, e2, er1, er2, num_seqs);
            }else
            {
                ;
            }
        }                
    }
}



int count = 0;   
void root(int *edge1, int *edge2, double *edge_length, double *d, int num_seqs, int **edge1_r, int **edge2_r)
{
    int *stat_nodes, index_max_edge, tip1, tip2;
    int nodes = 2*num_seqs - 2;
    int edges = 2*num_seqs - 3;
    int temp_chain1, temp_chain2;
    int temp1, temp2, temp_r;
    double halfway, temp_path_len = 0.0;

    int cpy_edge1[edges], cpy_edge2[edges];
    memcpy(cpy_edge1, edge1, sizeof(cpy_edge1));
    memcpy(cpy_edge2, edge2, sizeof(cpy_edge2));

    stat_nodes = statis_nodes(edge1, edge2, nodes);

    index_max_edge = index_max(d, num_seqs*(num_seqs - 1)/2 );
    flounder_sub(&tip1, &tip2, num_seqs, index_max_edge);
    
    while( !is_all_degree_than_1_except_tip1andtip2( tip1, tip2, stat_nodes , nodes) )
    {
        for( int i = 0; i<edges; i++)
        {
            if( cpy_edge1[i] == 0)
                continue;
            if( cpy_edge2[i] == tip1 || cpy_edge2[i] == tip2 )
                continue;
            temp1 = stat_nodes[ cpy_edge1[i] - 1 ];
            temp2 = stat_nodes[ cpy_edge2[i] - 1 ];

            if( temp1 == 1 || temp2 == 1)
            {
                stat_nodes[cpy_edge1[i] - 1]--;
                stat_nodes[cpy_edge2[i] - 1]--;
                cpy_edge1[i] = 0;
                cpy_edge2[i] = 0;
            }
        }
    }
    
    printf("tip1:%d, tip2:%d\n", tip1, tip2);

    halfway = d[index_max_edge] / 2.0;
    
    for( int i=0; i < num_seqs*(num_seqs - 1)/2; i++)
    {
        printf("%lf ", d[i]);
    }
    
    for( int i = 0; i<nodes; i++)
    {
        printf("cpy_e1:%d,cpy_e2:%d\n", cpy_edge1[i], cpy_edge2[i]);
    }

    printf("halfway:%lf\n", halfway);
    temp_chain1 = tip1;
    int DO = 1;
    while( DO )
    {
        for( int i = 0; i<edges; i++)
        {
            if( cpy_edge2[i] == 0)
                continue;
            if( cpy_edge1[i] == temp_chain1)
            {
                temp_path_len += edge_length[i];
                temp_chain1 = cpy_edge2[i];
                if( temp_path_len > halfway)
                {
                    temp_chain2 = cpy_edge1[i];
                    printf("Ceshi0\n");
                    DO = 0;
                    break;
                }
            }else if( cpy_edge2[i] == temp_chain1)
            {
                temp_path_len += edge_length[i];
                temp_chain1 = cpy_edge1[i];
                if( temp_path_len > halfway)
                {
                    temp_chain2 = cpy_edge2[i];
                    printf("Ceshi1\n");
                    DO = 0;
                    break;
                }
            }else
            {
                ;
            }
        }
    }

    printf("temp_path_len:%lf\n", temp_path_len);

    *edge1_r = (int *)malloc(nodes*sizeof(int));
    *edge2_r = (int *)malloc(nodes*sizeof(int));

    if( !edge1 || !edge2 || !edge_length)
    {
        printf("Malloc fail!\n");
        exit(1);
    }

    printf("num_seqs:%d\n", num_seqs);
    printf("cha1:%d, cha2:%d\n", temp_chain1, temp_chain2);

    int in1[2] = {0};
    int i_1 = 0;
    int in2[2] = {0};
    int i_2 = 0;

    for( int i = 0; i<edges; i++)
    {
        if( edge1[i] == temp_chain1 && edge2[i] != temp_chain2)
        {
            in1[i_1++] = edge2[i];
        }else if(edge2[i] == temp_chain1 && edge1[i] != temp_chain2)
        {
            in1[i_1++] = edge1[i];
        }else
        {
            ;
        }

        if( edge1[i] == temp_chain2 && edge2[i] != temp_chain1)
        {
            in2[i_2++] = edge2[i];
        }else if(edge2[i] == temp_chain2 && edge1[i] != temp_chain1)
        {
            in2[i_2++] = edge1[i];
        }else
        {
            ;
        }
    }

    (*edge1_r)[count] = num_seqs + 1;
    if( temp_chain1 <= num_seqs)
        (*edge2_r)[count++] = temp_chain1;
    else
        (*edge2_r)[count++] = temp_chain1 + 1;
    
    (*edge1_r)[count] = num_seqs + 1;
    if( temp_chain2 <= num_seqs)
        (*edge2_r)[count++] = temp_chain2;
    else
        (*edge2_r)[count++] = temp_chain2 + 1;

    if( in1[0] != 0)
        f( temp_chain1, in1[0], edge1, edge2, *edge1_r, *edge2_r, num_seqs);
    if( in1[1] != 0)
        f( temp_chain1, in1[1], edge1, edge2, *edge1_r, *edge2_r, num_seqs);
    if( in2[0] != 0)
        f( temp_chain2, in2[0], edge1, edge2, *edge1_r, *edge2_r, num_seqs);
    if( in2[1] != 0)
        f( temp_chain2, in2[1], edge1, edge2, *edge1_r, *edge2_r, num_seqs);

    free(stat_nodes);

    for( int i = 0; i<nodes; i++)
    {
        printf("er1:%d,er2:%d\n", (*edge1_r)[i], (*edge2_r)[i]);
    }
    printf("DDDD\n");
}