/*
 *	Example of a C program with an environment
 *
 *       gcc -g -o maine maine.c
 *     
 */


#include <stdio.h>		


int main(int argc,char *argv[], char *env[])
{
  int i;

  for (i=0; env[i] != 0; i++)
    printf("%d: %s\n",i,env[i]);

  return 0;
}
