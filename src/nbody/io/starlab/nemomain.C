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
 */

#include <stdinc.h>
#include <getparam.h>

// typedef char *nemo_string;        /* should be in: stdinc.h */

extern nemo_string defv[];	  /* defv MUST be defined in user program ! */
extern void nemo_main(void);	  /* this is the programmer's 'main' */

int main(int argc,char *argv[])
{
  initparam(argv,defv);		/* start  */
  nemo_main();			/* call his/her main program */
  finiparam();			/* end */
  return 0;                     /* return normal status to shell */
}
