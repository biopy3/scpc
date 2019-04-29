#include "correct_tip_number.h"

void correct_tip_number( Tree *tree, char *seq_name[])
{
  int i, j;
  if( tree == NULL)
  {
    printf("The pointer of tree is NULL.\n");
    exit(-1);
  }
  for( i = 0; i < tree->cterminals; i++)
  {
    if( tree->edge2[i] <= tree->cterminals)
    {
      for( j = 0; j < tree->cterminals; j++)
      {
        if( strcmp( tree->label[i], seq_name[j]) == 0)
        {
          tree->edge2[i] = j + 1;
          break;
        }
        if( j == tree->cterminals - 1)
        {
          printf("The given sequence name in fasta file is not same label of tips in tree\n.");
          exit(-1);
        }
      }

    }
  }
}