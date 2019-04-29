#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "newick.h"
#include "stack.h"

static long len_byte;
static char * byte_string;
static long pos = 0;
int cterminals = 1; /* the number of terminals equel number of comma plus one */
int cinnodes = 0;

int table[ASCII_SIZE] = {
/*      */   IGNORED,   IGNORED, IGNORED,   IGNORED,
/*      */   IGNORED,   IGNORED, IGNORED,   IGNORED,
/*      */   IGNORED,WHITESPACE, WHITESPACE,   IGNORED,
/*      */   IGNORED,WHITESPACE, IGNORED,   IGNORED,
/*      */   IGNORED, IGNORED, IGNORED,   IGNORED,
/*      */   IGNORED, IGNORED, IGNORED,   IGNORED,
/*      */   IGNORED, IGNORED, IGNORED,   IGNORED,
/*      */   IGNORED, IGNORED, IGNORED,   IGNORED,
/*  !"# */WHITESPACE,   LABEL,   LABEL,     LABEL,
/* $%&' */     LABEL,   LABEL,   LABEL,     LABEL,
/* ()*+ */    OPAREN,  CPAREN,   LABEL,      PLUS,
/* ,-./ */     COMMA,    DASH,     DOT,     LABEL,
/* 0123 */     DIGIT,   DIGIT,   DIGIT,     DIGIT,
/* 4567 */     DIGIT,   DIGIT,   DIGIT,     DIGIT,
/* 89:; */     DIGIT,   DIGIT,   COLON, SEMICOLON,
/* <=>? */     LABEL,   LABEL,   LABEL,     LABEL,
/* @ABC */     LABEL,    LABEL,    LABEL,      LABEL,
/* DEFG */      LABEL,    LABEL,    LABEL,      LABEL,
/* HIJK */      LABEL,    LABEL,    LABEL,      LABEL,
/* LMNO */      LABEL,    LABEL,    LABEL,      LABEL,
/* PQRS */      LABEL,    LABEL,    LABEL,      LABEL,
/* TUVW */      LABEL,    LABEL,    LABEL,      LABEL,
/* XYZ[ */      LABEL,    LABEL,    LABEL,   LSQUARE,
/* \]^_ */     LABEL, RSQUARE,   LABEL,     LABEL,
/* `abc */     LABEL,    LABEL,    LABEL,      LABEL,
/* defg */      LABEL,    LABEL,    LABEL,      LABEL,
/* hijk */      LABEL,    LABEL,    LABEL,      LABEL,
/* lmno */      LABEL,    LABEL,    LABEL,      LABEL,
/* pqrs */      LABEL,    LABEL,    LABEL,      LABEL,
/* tuvw */      LABEL,    LABEL,    LABEL,      LABEL,
/* xyz{ */      LABEL,    LABEL,    LABEL,     LABEL,
/* |}~  */     LABEL,   LABEL,   LABEL,    IGNORED
 };

char *read_file( const char *file_name, long *file_size)
{
  FILE *f;
  char *rawdata;
  int temp_char;

  printf( "Reading newick file %s ...\n", file_name);
  f = fopen( file_name, "rb");
  if( f == NULL)
  {
    printf( "Can't open file!\n");
    return(NULL);
  }
  
  /* get file size */
  if( fseek( f, 0, SEEK_END) == -1)
  {
    printf( "fseek the file :%s was wrong!\n", file_name);
    fclose(f);
    return(NULL);
  }
  
  *file_size = ftell(f);
  
  if( *file_size == -1)
  {
    printf( "Wrong, the file size is -1, %s\n", file_name);
    fclose(f);
    return(NULL);
  }
  
  rewind(f);
  while( ( temp_char = fgetc(f)) != EOF)
  {
    if( temp_char == ',')
    {
      cterminals++;
    }
    if( temp_char == '(')
    {
      cinnodes++;
    }
  }

  rawdata = ( char *)malloc( ( (*file_size)+1) * sizeof(char));
  if( rawdata != NULL)
  {
    rewind(f);
    if( fread( rawdata, sizeof(char), *file_size, f) != (size_t) *file_size)
    {
      printf( "The file size is not same!\n");
      free(rawdata);
      rawdata = NULL;
    }
    else
    {
      rawdata[*file_size] = 0;
    }
    
  }else
  {
    printf("Allocate the memory fail for rawdate!\n");
  }

  fclose(f);
  printf("Reading success!\n");
  return(rawdata);
}

int get_btye( void)
{
  if( pos == len_byte)
  {
    return(EOF);
  }
  
  return( table[ byte_string[pos]]);
}

TreeItem get_tree_item( TreeItem *prev_tree_item)
{
  long start_pos = pos;
  TreeItem tree_item;
  int temp_char;

  temp_char = get_btye();

  if( temp_char == EOF)
  {
    tree_item.type = EOF;
  }
  else if( temp_char == WHITESPACE)
  {
    tree_item.type = IGNORED;
    pos++;
  }
  else if( temp_char == LSQUARE)
  {
    tree_item.type = IGNORED;
    pos++;
    while( ( temp_char = get_btye()) != RSQUARE)
    {
      if( pos >= ASCII_SIZE)
      {
        printf("The comments start '[' and must end ']', can't find ']'.\n");
        exit(-1);
      }
      pos++;
    }
    pos++;
  }
  else if( temp_char == IGNORED)
  {
    tree_item.type = IGNORED;
    pos++;
  }
  else if( temp_char == OPAREN)
  {
    tree_item.type = OPAREN;
    pos++;
  }
  else if( temp_char == CPAREN)
  {
    tree_item.type = CPAREN;
    pos++;
  }
  else if( temp_char == COLON)
  {
    tree_item.type = COLON;
    pos++;
  }
  else if( temp_char == COMMA)
  {
    tree_item.type = COMMA;
    pos++;
  }
  else if( temp_char == SEMICOLON)
  {
    tree_item.type = SEMICOLON;
    pos++;
  }
  else if( temp_char == LABEL || temp_char == DIGIT ||
           temp_char == PLUS || temp_char == DASH ||
           temp_char == DOT)
  {
    if( (*prev_tree_item).type == COLON)
    {
      if( temp_char == DASH || temp_char == PLUS || temp_char == DIGIT)
      {
        pos++;
        while( (temp_char = get_btye()) == DIGIT)
        {
          pos++;
        }

        if( temp_char == DOT)
        {
          pos++;
          while( (temp_char = get_btye()) == DIGIT)
          {
            pos++;
          }  
          tree_item.type = BRANCH_LENGTH;
          tree_item.lexeme = byte_string + start_pos;
          tree_item.len = pos - start_pos - 1;
        }
        else if( byte_string[pos] == 'E' || byte_string[pos] == 'e')
        {
          pos++;
          temp_char = get_btye();
          if( temp_char == PLUS || temp_char == DASH || temp_char == DIGIT)
          {
            pos++;
            while( (temp_char = get_btye()) == DIGIT)
            {
              pos++;
            }
            tree_item.type = BRANCH_LENGTH;
            tree_item.len = pos - start_pos - 1;
            tree_item.lexeme = byte_string + start_pos;
          }
          else
          {
            printf("The newick file format is invalid, because the after of 'E' or 'e' in branch_length \n\
  must be digits, even if nothing, it's always required.\n");
            exit(-1);
          }

        }
        else 
        {
          tree_item.type = BRANCH_LENGTH;
          tree_item.len = pos - start_pos - 1;
          tree_item.lexeme = byte_string + start_pos;
        }
      }
      else
      {
        printf("The newick file format is invalid, because the first char of branch_length \n\
  must be '-', '+', digits or whitespace, and the after of colon character must \n\
  have branch length, even if 0.0, it's always required.\n");
        exit(-1);
      }
    }
    else
    {
      pos++;
      while( (temp_char = get_btye()) == LABEL || temp_char == DIGIT ||
              temp_char == DOT || temp_char == PLUS || temp_char == DASH)
      {
        pos++;
      }

      tree_item.type = LABEL;
      tree_item.len = pos - start_pos - 1;
      tree_item.lexeme = byte_string + start_pos;
    }
  }
  else
  {
    printf("The newick file format is invalid.\n");
    exit(-1);
  }

  return tree_item;
}

Tree *parse_newick_string()
{

  TreeItem tree_item;
  TreeItem prev_tree_item;
  Tree *tree;
  
  prev_tree_item.type = FIRST_TREE_ITEM;
  tree = ( Tree *)calloc( 1, sizeof(Tree));
  
  tree->cterminals = cterminals;
  tree->label = ( char **)calloc( ( 2*cterminals - 2), sizeof(char *));
  tree->edge1 = ( int *)calloc( ( 2*cterminals - 2), sizeof(int));
  tree->edge2 = ( int *)calloc( ( 2*cterminals - 2), sizeof(int));
  tree->branch = ( char **)calloc( ( 2*cterminals - 2), sizeof(char *));

  int nodes;
  int nop = 0; /* number of parethese*/
  int count_innode = cterminals;
  int count_tip = 0;
  int index_innode = 0;

  Stack *stack1 = NULL;
  Stack *stack2 = NULL;

  while( ( tree_item = get_tree_item( &prev_tree_item)).type != EOF)
  {
    if( tree_item.type == IGNORED)
    {
      continue;
    }
    
    if( tree_item.type == OPAREN)
    {
      if( prev_tree_item.type != COMMA &&
          prev_tree_item.type != OPAREN &&
          prev_tree_item.type != FIRST_TREE_ITEM)
      {
        printf("The newick file format is invalid, because the ahade of the\n\
  open parenthese must be open parenthese, comma or nothing.\n");
        exit(-1);
      }
      nop++;
      count_innode++;
      StackPush( &stack1, count_innode);
      StackPush( &stack2, count_innode);
    }
    else if( tree_item.type == LABEL)
    {
      if( prev_tree_item.type != OPAREN &&
          prev_tree_item.type != COMMA &&
          prev_tree_item.type != CPAREN)
      {
        printf("The newick file format is invalid, because the ahade of the\n\
  label must be open parenthese, comma or close parenthese.\n");
        exit(-1);
      }

      if( nop <= 0)
      {
        break;
      }

      if( prev_tree_item.type == CPAREN)
      {
        (tree->edge1)[index_innode] = stack1->chain;
        (tree->edge2)[index_innode] = stack2->chain;
        StackPop( &stack2);
      }
      else
      {
        count_tip++;
        (tree->edge1)[index_innode] = stack1->chain;
        (tree->edge2)[index_innode] = count_tip;
      }

      (tree->label)[index_innode] = ( char *)calloc( ( tree_item.len + 2), sizeof( char));
      if( (tree->label)[index_innode] == NULL)
      {
        printf("Can't allocate memory for tree label.\n");
        exit(-1);
      }
      memcpy( (tree->label)[index_innode], tree_item.lexeme, ( tree_item.len + 1)*sizeof(char));
      index_innode++;
    }
    else if( tree_item.type == COMMA)
    {
      if( prev_tree_item.type != LABEL &&
          prev_tree_item.type != COMMA &&
          prev_tree_item.type != CPAREN &&
          prev_tree_item.type != BRANCH_LENGTH &&
          prev_tree_item.type != OPAREN)
      {
        printf("The newick file format is invalid, because the ahade of the\n\
  comma must be open parenthese, comma, branch length, label or\n\
  close parenthese.\n");
        exit(-1);
      }
      if( prev_tree_item.type == CPAREN)
      {
        (tree->edge1)[index_innode] = stack1->chain;
        (tree->edge2)[index_innode] = stack2->chain;
        StackPop( &stack2);
        index_innode++;
      }
      else if( prev_tree_item.type == OPAREN || 
               prev_tree_item.type == COMMA)
      {
        count_tip++;
        (tree->edge1)[index_innode] = stack1->chain;
        (tree->edge2)[index_innode] = count_tip;
        index_innode++;
      }
      else
      {
        ;
      }
    }
    else if( tree_item.type == COLON)
    {
      if( prev_tree_item.type != LABEL &&
          prev_tree_item.type != COMMA &&
          prev_tree_item.type != CPAREN &&
          prev_tree_item.type != OPAREN)
      {
        printf("The newick file format is invalid, because the ahade of the\n\
  colon must be open parenthese, comma, label or  close parenthese.\n");
        exit(-1);
      }

      if( nop <= 0)
      {
        break;
      }
      if( prev_tree_item.type == CPAREN)
      {
        (tree->edge1)[index_innode] = stack1->chain;
        (tree->edge2)[index_innode] = stack2->chain;
        StackPop( &stack2);
        index_innode++;
      }
      else if( prev_tree_item.type != LABEL)
      {
        count_tip++;
        (tree->edge1)[index_innode] = stack1->chain;
        (tree->edge2)[index_innode] = count_tip;
        index_innode++;
      }
      else
      {
        ;
      }
    }
    else if( tree_item.type == BRANCH_LENGTH)
    {
      if( prev_tree_item.type != COLON)
      {
        printf("The newick file format is invalid, because the ahade of the\n\
                branch length must be colon.\n");
        exit(-1);
      }

      (tree->branch)[index_innode -1] = ( char *)calloc( ( tree_item.len + 2), sizeof(char));
      if( (tree->branch)[index_innode - 1] == NULL)
      {
        printf("Can't allocate memory for tree branch.\n");
        exit(-1);
      }
      memcpy( (tree->branch)[index_innode - 1], tree_item.lexeme, ( tree_item.len + 1)*sizeof(char));
    }
    else if( tree_item.type == CPAREN)
    {
      if( prev_tree_item.type != LABEL &&
          prev_tree_item.type != COMMA &&
          prev_tree_item.type != CPAREN &&
          prev_tree_item.type != BRANCH_LENGTH)
      {
        printf("The newick file format is invalid, because the ahade of the\n\
                close parenthese must be comma, branch length, label or\n\
                close parenthese.\n");
        exit(-1);
      }
      if( nop <= 0)
      {
        printf("The newick file format is invalid, the number of '(' not same ')'\n");
        exit(-1);
      }
      
      if( prev_tree_item.type == CPAREN)
      {
        (tree->edge1)[index_innode] = stack1->chain;
        (tree->edge2)[index_innode] = stack2->chain;
        StackPop( &stack2);
        index_innode++;
      }
      else if( prev_tree_item.type == COMMA)
      {
        count_tip++;
        (tree->edge1)[index_innode] = stack1->chain;
        (tree->edge2)[index_innode] = count_tip;
        index_innode++;
      }
      else
      {
        ;
      }
      ;
      nop--;
      StackPop( &stack1);
      
    }
    else if( tree_item.type == SEMICOLON)
    {
      if( prev_tree_item.type != LABEL &&
          prev_tree_item.type != CPAREN &&
          prev_tree_item.type != BRANCH_LENGTH)
      {
        printf("The newick file format is invalid, because the ahade of the\n\
                semicolon must be laber, branch length or close parenthese.\n");
        exit(-1);
      }
    }
    else
    {
      ;
    }
    
    memcpy( &prev_tree_item, &tree_item, sizeof(TreeItem));
  }

  return(tree);
}

Tree * parse_newick_file( const char * file_name)
{
  long file_size;
  char *raw_string;
  Tree *tree;

  raw_string = read_file( file_name, &file_size);
  if( raw_string == NULL)
  {
    fprintf (stderr, "Error while opening/reading file %s\n", file_name);
    return(NULL);
  }
  else
  {
    len_byte = strlen(raw_string);
    byte_string = raw_string;
    tree = parse_newick_string();
    free(raw_string);
    byte_string = NULL;
  }

  return(tree);
}

