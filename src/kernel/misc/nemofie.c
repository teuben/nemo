/*
 * NEMOFIE:  NEMO driver for GIPSY fie routines
 *
 *	18-dec-88  created - Peter Teuben
 *	19-jun-89  adapted to F2C interface definitions - see nemofie(3) PJT
 *	18-feb-92  adapted to new 'real' usage in fie()			 PJT
 *       7-mar-92  happy gcc2.0, inifien is int, not void.               PJT
 *      16-dec-95  proper externs defined				 PJT
 */

#include <stdinc.h>


extern int inifie(string);
extern void dofie(real *, int *, real *, real *);
extern void dmpfie(void);


int inifien (string expr)
{
    return inifie(expr);
}

void dofien (real *pars, int n, real *result, real errval)
{
    dofie(pars, &n, result, &errval);
}

void dmpfien()
{ 
    dmpfie();
}
