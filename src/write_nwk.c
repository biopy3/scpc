#include "scpc.h"
#include "write_nwk.h"

extern int cinnodes;
extern int cterminals;

void write( int *f_er1, int *f_er2, FILE *f, int ROOT, char *seq_name[])
{
	int i = 0;
    fwrite("(" , sizeof(char), 1, f);
	for( i = 0; i < cinnodes + cterminals -1; i++)
    {
        if( f_er1[i] == ROOT)
        {
            if( 1 <= f_er2[i] && f_er2[i] <= cterminals )
            {
                fprintf(f, "%s", seq_name[ f_er2[i] - 1]);
                fwrite("," , sizeof(char), 1, f);
            }
            else if( f_er2[i] > cterminals)
            {
                write( f_er1, f_er2, f, f_er2[i], seq_name);
            }
            else
            {
                ;
            }
        }
    }
    fseek( f, -1L, 1);
    fwrite(")" , sizeof(char), 1, f);
    fwrite("," , sizeof(char), 1, f);
}


void write_nwk( int *f_er1, int *f_er2, char *seq_name[], char *output)
{  
    int ROOT = cterminals + 1;
    FILE *f_nwk;
    f_nwk = fopen(output,"w");
    if( f_nwk == NULL)
    {
        printf("Can't create the species newick file!\n");
        exit(1);
    }
    write( f_er1, f_er2, f_nwk, ROOT, seq_name);
    fseek( f_nwk, -1L, 1);
    fwrite(";" , sizeof(char), 1, f_nwk);
    fclose(f_nwk);
}

void write_species_table(int *e1, int *e2, char *seq_name[], char *output)
{
    FILE *f;
    f = fopen(output, "w");
    int i, j;
    if( f == NULL)
    {
        printf("Can't create the species table file!\n");
        exit(1);
    }

    for( i = 0; i < cinnodes; i++)
    {
        for( j = 0; j < cinnodes + cterminals - 1; j++)
        {
            if( e1[j] - cterminals - 1 == i && e2[j] <= cterminals)
            {
                fprintf(f, "%s,", seq_name[ e2[j] - 1]);
            }
        }
        fseek( f, -1L, 1);
        fprintf(f, "\n");
    }
    fclose(f);
}