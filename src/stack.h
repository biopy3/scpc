#ifndef __STACK_H__
#define __STACK_H__

struct Stack
{
  int chain;
  struct Stack *next;
};

typedef struct Stack Stack;

void StackClear( Stack **stack);
int StackPop( Stack **head);
void StackPush( Stack **head, int chain);
int StackSize( Stack **stack);

#endif