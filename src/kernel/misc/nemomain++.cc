/*
 *	NEMOMAIN:  standard way a nemoprogram starts and finishes
 *	by calling your main program 'nemo_main' instead of
 *	'main' you get more of the user interface done
 *	transparently
 *              (1) calls 'initparam' to set up program
 *              (2) calls program 'nemo_main'
 *              (3) calls 'finiparam' to finish program
 *              (4) return normal status to shell
 *
 *      20-nov-93  added options MAIN_ symbol to fool the fortran compiler
 *	21-nov-94  correct C++ version for C++ main
 *	20-jun-01  gcc3, must return an int value
 */

#include <stdinc.h>
#include <getparam.h>

extern string defv[];		/* defv MUST be defined in user program ! */
extern string usage;		/* One line description of the program    */

extern void nemo_main(void);	/* this must be be a C++ nemo_main too !! */

int main(int argc,char *argv[])
{
    if(argv[argc] != NULL) 
       warning("Your C compiler seems to have an old-style (short) argv[]");
    initparam(argv,defv);		/* start                         */
    nemo_main();			/* call program                  */
    finiparam();			/* end                           */
    exit(0);                            /* return normal status to shell */
    return 0;                           /* also need return value        */
}

#if defined(MAIN)
int MAIN_()
{
  error("MAIN_ called; some fortran inconsistency");
  return 0;
}
#endif
