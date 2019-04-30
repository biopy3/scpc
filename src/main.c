#include "scpc.h"
#include "newick.h"
#include "parse_fasta2bin.h"
#include "opra_matrix.h"
#include "correct_tip_number.h"
#include "dist_dna.h"
#include "statistic.h"
#include "write_nwk.h"
#include "fixing_tree.h"

#define _FILE_NAME_LEN 260
#define _MAX_FILES 99
#define MODEL_LEN 5
#define VERSION "0.2"
#define MAX_SEQ_NUM 999

double *d;
int main( int argc, char *argv[])
{
  char fasta_file[FILENAME_MAX] = { 0 };
  char fas_file_name[_FILE_NAME_LEN] = {0};
  char newick_file[FILENAME_MAX] = { 0 };
  char output_dir[FILENAME_MAX] = { 0 };
  char model[MODEL_LEN] = { 0 };
  int output_dir_len = 0;
  int fas_filename_len = 0;
  int count_arg = 1;
  int i, j, k;

  char prev_arg = 0;
  int cout_for_fas_filename = -1;

  if( argc == 1)
  {

    printf("Usage: scpc -s <alignment> -t <newick_tree> -m <model_name>\n\
explain:\n\
  \n\
  -s | --mseqs <alignment> Input alignment in FASTA format \n\
  -t | --tree <newick_tree> Input tree in NEWICK format \n\
  -m | --model <string> Indicate the substituted model \n\
  -o | --outdir <path> Indicate where the directory of results to put \n\
  -v | --version Print product version and exit \n");
  }

  while( count_arg < argc)
  {
    if( argv[count_arg][0] == '-' && argv[count_arg][1] != '-' &&
        argv[count_arg][2] == 0)
    {
      switch ( argv[count_arg][1])
      {
        case 's':
          prev_arg = 's';
          break;
        case 't':
          prev_arg = 't';
          break;
        case 'm':
          prev_arg = 'm';
          break;
        case 'o':
          prev_arg = 'o';
          break;
        case 'v':
          printf("The current scpc version is %s.\n", VERSION);
          break;
        default:
          printf("The command isn't supply.\n");
          return EXIT_FAILURE;
      }
    }
    else if( strcmp(argv[count_arg], "--mseqs") == 0)
    {
      prev_arg = 's';
    }
    else if( strcmp(argv[count_arg], "--model") == 0)
    {
      prev_arg = 'm';
    }
    else if( strcmp(argv[count_arg], "--outdir") == 0)
    {
      prev_arg = 'o';
    }
    else if( strcmp(argv[count_arg], "--tree") == 0)
    {
      prev_arg = 't';
    }
    else if( strcmp(argv[count_arg], "--version") == 0)
    {
      printf("The current scpc version is %s.\n", VERSION);
    }
    else
    {
      switch(prev_arg)
      {
        case 's':
          strcpy( fasta_file, argv[count_arg]);
          fas_filename_len = strlen( fasta_file);
          printf("fas_filename_len:%d\n",fas_filename_len);
          for( i = fas_filename_len; -1 < i; i--)
          {
            if( fasta_file[i] == '.')
            {
              while( !(fasta_file[i] == PATH_SPLIT || i < 0))
              {
                cout_for_fas_filename++;
                i--;
              }
              memcpy( fas_file_name, fasta_file + i + 1, cout_for_fas_filename*sizeof(char));
              cout_for_fas_filename = -1;
              break;
            }
          }
          
          prev_arg = '\0';
          break;
        case 't':
          strcpy( newick_file, argv[count_arg]);
          prev_arg = '\0';
          break;
        case 'm':
          strcpy( model, argv[count_arg]);
          prev_arg = '\0';
          break;
        case 'o':
          strcpy( output_dir, argv[count_arg]);
          output_dir_len = strlen( output_dir);
          if( output_dir[ output_dir_len - 1] != PATH_SPLIT)
          {
            output_dir[ output_dir_len] = PATH_SPLIT;
            output_dir[ output_dir_len + 1] = '\0';
            output_dir_len++;
          }
          prev_arg = '\0';
          break;
        case '\0':
          printf("The given parameters is invaild, please check then try again.\n");
          return EXIT_FAILURE;
        default:
          printf("The command isn't supply.\n");
          return EXIT_FAILURE;
      }
    }
    count_arg++;
  }


  Tree *t;
  t = parse_newick_file( newick_file);
  printf( "The cteminals: %d\n", t->cterminals);
  // for( int i = 0; i < 2*(t->cterminals) - 2; i++)
  // {
  //   printf("edge1:%d\n", t->edge1[i]);
  //   printf("edge2:%d\n", t->edge2[i]);
  //   printf("branch:%s\n", t->branch[i]);
  //   printf("label:%s\n\n", t->label[i]);
  // }

  extern int cterminals;
  // printf("cterminals:%d\n", cterminals);
  int num_distances = 0;
  int len_seq[cterminals], len_name[cterminals], konw_bases = 0;

  int sta = chec_file( fasta_file, cterminals, len_seq, len_name, &konw_bases);
	if(sta != 0)
  {
      printf("check file fail!\n");
      return EXIT_FAILURE;
  }
  
  char *seq_name[cterminals];

  unsigned char *Seqs = init_seq(cterminals, len_seq[0], len_name);
  init_seq_name(seq_name, cterminals, len_name);
  
  
  int status = raw2DNAbin(Seqs, seq_name, fasta_file, len_seq[0]);
  if(status != 0 )
  {
      printf("exit!\n");
      exit(1);
  }
  // for( int i = 0; i < cterminals; i++)
  // {
  //   printf("seq_name[%d]:%s", i, seq_name[i]);
  // }
  
  correct_tip_number( t, seq_name);
  transpose(Seqs, cterminals, len_seq[0]);
  d = dist_dna_api( Seqs, model, 0, 0, 0, 0, &cterminals, &num_distances,
                    len_seq[0], konw_bases);

  extern int cinnodes;
  // printf("cinnodes:%d\n", cinnodes);

  double sta_d[cinnodes];
  int parent[cinnodes + cterminals];

  // for( int i = 0; i < 2*cterminals -2; i++)
  // {
  //   printf("edge1:%d, edge2:%d\n", t->edge1[i], t->edge2[i]);
  // }
  
  int *stat;
  int bins;
  double stat_max, stat_min;
  double gap = get_gap( t, parent, sta_d, &stat, &bins, &stat_max, &stat_min);

  /* write distance matrix */
  int index_d;
  char output_dmat[ FILENAME_MAX] = {0};
  FILE *f_dmat;
  
  if( output_dir_len == 0)
  {
    strcat( output_dmat, fasta_file);
    f_dmat = fopen( strcat( output_dmat, ".dmat.csv"), "w");
  }
  else
  {
    strcat( output_dmat, output_dir);
    strcat( output_dmat, fas_file_name);
    f_dmat = fopen( strcat( output_dmat, "_dmat.csv"), "w");
  }
  if( f_dmat == NULL)
  {
    printf("Can't open the file:%s\n", output_dmat);
    return EXIT_FAILURE;
  }
    /* write first line */
  fprintf( f_dmat, "\t,");
  for( int i0 = 0; i0 < cterminals - 1; i0++)
  {
    fprintf( f_dmat, "%s,", seq_name[i0]);
  }
  fprintf( f_dmat, "%s\n", seq_name[cterminals - 1]);

    /* write the all left line */
  for( int i1 = 0; i1 < cterminals; i1++)
  {
    fprintf( f_dmat, "%s,", seq_name[i1]);
    for( int i2 = 0; i2 < cterminals; i2++)
    {
      if( i1 == i2)
      {
        fprintf( f_dmat, "%lf,", 0.000000);
      }
      else
      {
        get_d_index( i1+1, i2+1, cterminals, &index_d);
        fprintf( f_dmat, "%lf,", d[index_d]);
      }
    }
    fseek( f_dmat, -1L, 1);
    fprintf( f_dmat, "\n");
  }

  fclose( f_dmat);
  free(d);

  /* fixing tree */
  fixing_tree( t, gap, sta_d, parent);

  /* Write the species tree */
  char output_species_tree[ FILENAME_MAX] = {0};
  char output_species_table_tree[ FILENAME_MAX] = {0};
  if( output_dir_len == 0)
  {
    strcat( output_species_tree, fasta_file);
    strcat( output_species_table_tree, fasta_file);
    write_nwk( t->edge1, t->edge2, seq_name, strcat( output_species_tree, "_species.nwk"));
    write_species_table(t->edge1, t->edge2, seq_name, strcat( output_species_table_tree, "_table_species.csv"));
  }
  else
  {
    strcat( output_species_tree, output_dir);
    strcat( output_species_table_tree, output_dir);
    write_nwk( t->edge1, t->edge2, seq_name, strcat(strcat( output_species_tree, fas_file_name), "_species.nwk"));
    write_species_table(t->edge1, t->edge2, seq_name, strcat(strcat( output_species_table_tree, fas_file_name), "_table_species.csv"));
  }
  
  // for( int i = 0; i < 2*cterminals -2; i++)
  // {
  //   printf("edge1:%d, edge2:%d\n", t->edge1[i], t->edge2[i]);
  // }

  /* write stat */
  char output_stat[ FILENAME_MAX] = {0};
  if( output_dir_len == 0)
  {
    strcat( output_stat, fasta_file);
    FILE *f_stat = fopen( strcat( output_stat, "_stat.csv"),"w");
    if( f_stat == NULL)
    {
      printf("Can't open the file:%s\n", fasta_file);
      return EXIT_FAILURE;
    }

    for( j = 0; j < bins - 1; j++)
    {
      fprintf( f_stat, "%d", stat[j]);
      fputc(',', f_stat);
    }
    fprintf( f_stat, "%d\nstat_min:%lf,stat_max:%lf,gap_value:%lf", stat[bins], stat_min, stat_max, gap);
  }
  else
  {
    strcat( output_stat, output_dir);
    FILE *f_stat = fopen(strcat(strcat( output_stat, fas_file_name), "_stat.csv"),"w");
    
    if( f_stat == NULL)
    {
      printf("Can't open the file:%s, maybe some directories not existed in the path.\n", output_dir);
      return EXIT_FAILURE;
    }

    for( j = 0; j < bins - 1; j++)
    {
      fprintf( f_stat, "%d", stat[j]);
      fputc(',', f_stat);
    }
    fprintf( f_stat, "%d\nstat_min:%lf,stat_max:%lf,gap_value:%lf", stat[bins], stat_min, stat_max, gap);
  }

  /* write the internal nodes distance */
  if( output_dir_len == 0)
  {
    FILE *f_sta_d = fopen( strcat( fasta_file, "_innode_dist.csv"),"w");
    if( f_sta_d == NULL)
    {
      printf("Can't open the file:%s\n", fasta_file);
      return EXIT_FAILURE;
    }

    for( j = 0; j < cinnodes - 1; j++)
    {
      fprintf( f_sta_d, "%lf", sta_d[j]);
      fputc(',', f_sta_d);
    }
    fprintf( f_sta_d, "%lf", sta_d[cinnodes]);
  }
  else
  {
    FILE *f_sta_d = fopen(strcat(strcat( output_dir, fas_file_name), "_innode_dist.csv"),"w");
    
    if( f_sta_d == NULL)
    {
      printf("Can't open the file:%s, maybe some directories not existed in the path.\n", output_dir);
      return EXIT_FAILURE;
    }

    for( j = 0; j < cinnodes - 1; j++)
    {
      fprintf( f_sta_d, "%lf", sta_d[j]);
      fputc(',', f_sta_d);
    }
    fprintf( f_sta_d, "%lf", sta_d[cinnodes]);
  }
  
  /* free tree*/
  for( k = 0; k < 2*cterminals - 2; k++)
  {
    free( t->branch[k]);
    free( t->label[k]);
  }
  free( t->edge1);
  free( t->edge2);
  free( t);

  return EXIT_SUCCESS;
}
