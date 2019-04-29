#include <stdio.h>
#include <stdlib.h>
#include "stack.h"

int StackSize( Stack **stack)
{
  Stack * top;
  int size = 0;
  top = *stack;

  while( top)
  {
    size++;
    top = top->next;
  }

  return( size);
}

void StackPush( Stack **head, int chain)
{
  Stack *new;

  new = ( Stack *)malloc( sizeof( Stack));
  if( new == NULL)
  {
    printf("Allocate memory is fail for new stack's chain!\n");
    exit(-1);
  }

  new->chain = chain;
  new->next = *head;
  *head = new;

}

int StackPop( Stack **head)
{
  int chain;
  Stack *temp;
  if( *head == NULL)
  {
    return(0);
  }

  temp = (*head);
  chain = (*head)->chain;
  (*head) = (*head)->next;
  free(temp);

  return(chain);
}

void StackClear( Stack **stack)
{
  while( *stack)
  {
    StackPop( stack);
  }
}