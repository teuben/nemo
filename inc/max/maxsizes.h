/*
 * maxsizes.h:    your optional systemwide MAX* macros 
 *    this file should reside in $NEMOLIB, though the template 
 *    lives in $NEMOINC/max and will be copied during the 
 *    installation process
 *
 *    Programmers are encourages to assemble all MAX* macros in
 *    this single file and document their usage
 */
   

/* stropen.c */

#if !defined(NEMO_MAXFD)
#define NEMO_MAXFD  64
#endif


#ifndef MAX_COL
#define MAX_COL 256
#endif

#ifndef MAX_LINELEN
#define MAX_LINELEN  16384
#endif


/* getfunc.c */

#ifndef MAXN
#define MAXN  16384
#endif

/* math/nllsqfit.c:#define MAXPAR  32  */

/* math/odeint.c:#define MAXSTP  10000		*/
/* math/odeint.c:#define NMAX	10		*/

/* misc/fie.c:#define MAXPAR  32 */

/* misc/layout.c:#define MAXWORD 256 */

/* tab/gettab.c:#define MAXCOL   64 */

/* tab/tabdms.c:#define MAXCOL          256  */

/* tab/tablines.c:#define MAXLINELEN  2048 */

/* tab/tabplot.c:#define MAXYCOL 256 */
/* tab/tabplot.c:#define MAXCOORD 16 */


