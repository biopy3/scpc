#include <ctype.h>
#include "correct_tip_number.h"

extern int cinnodes;

char *trimwhitespace(char *str)
{
  char *end;

  // Trim leading space
  while(isspace((unsigned char)*str))
  {
    str++;
  }
  if(*str == 0)  // All spaces?
  {
    return str;
  }
  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace((unsigned char)*end))
  {
    end--;
  }
  // Write new null terminator character
  end[1] = '\0';

  return str;
}

void correct_tip_number( Tree *tree, char *seq_name[])
{
  int i, j;
  if( tree == NULL)
  {
    printf("The pointer of tree is NULL.\n");
    exit(-1);
  }
  for( i = 0; i < tree->cterminals + cinnodes - 1; i++)
  {
    if( tree->edge2[i] <= tree->cterminals)
    {
      for( j = 0; j < tree->cterminals; j++)
      {
        tree->label[i] = trimwhitespace(tree->label[i]);
        seq_name[j] = trimwhitespace(seq_name[j]);
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