/*
 *	Example of a C program, only using NEMO I/O
 *      but not the initparam() because for example
 *      your package has its own argument parsing
 *
 *       On *Unix* to be compiled and linked as:
 *       
 *       cc -g -o mainc2 mainc2.c $NEMOLIB/libnemo.a -lm
 *     
 *       The flags '-g' are defensive programming flags
 *       and are optional (debugging in this case)
 */

#include <stdinc.h>
#include <filestruct.h>

#define NMAX  10

int main(int argc,char *argv[])
{
  int  nval = 15;
  stream ostr;

  initparam0("testing mainc2");


  ostr = stropen("mainc2.out","w");

  fwrite(&nval, sizeof(int), 1, ostr);

  strclose(ostr);
  return 0;
}
