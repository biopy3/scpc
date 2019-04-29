#ifndef __NEWICK_H__
#define __NEWICK_H__

#define ASCII_SIZE 128

#define PLUS '+'
#define DASH '-'
#define DIGIT 128
#define DOT '.'
#define LSQUARE '['
#define RSQUARE ']'
#define LABEL 129
#define WHITESPACE 130

#define FIRST_TREE_ITEM 131
#define IGNORED 132
#define OPAREN '('
#define CPAREN ')'
#define COLON ':'
#define COMMA ','
#define SEMICOLON ';'
#define BRANCH_LENGTH 133

typedef struct
{
  int type;
  const char *lexeme;
  int len;
} TreeItem;


typedef struct
{
  char **label;                                           
  int cterminals;
  char **branch;
  int *edge1;
  int *edge2;
} Tree;

/* API */
Tree * parse_newick_file( const char * file_name);
#endif