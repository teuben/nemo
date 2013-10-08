/* 
 * File:  myline.c  
 * taken from the man page for the Testfile 
 */

#include <stdinc.h>

real func_line(real *x, real *p, int np)
{
  return p[0] + p[1]*x[0];
}

void derv_line(real *x, real *p, real *e, int np)
{
  e[0] = 1.0;
  e[1] = x[0];
}
