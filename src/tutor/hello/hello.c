/*
 *  HELLO:  Example of ``hello world'' using NEMO in C/C++
 *
 *      ??         V1.0    Created                 Peter Teuben
 *	23-nov-94  V1.3    C and C++ code now same 	    PJT
 *                    (but they need different file extensions)
 *      17-jun-01  updated source code documentation for NEMO V3
 */

#include <stdinc.h>              /* standard (NEMO) definitions */
#include <getparam.h>                         /* user interface */

string defv[] = {       /* standard keywords and default values */
    "verbose=true\n            Verbosity level (t|f)",  /* key1 */
    "VERSION=1.3\n             23-nov-94 PJT",          /* key2 */
    NULL,               /* standard terminator of defv[] vector */
};

string usage = "Example NEMO C/C++ program";      /* usage text */

void nemo_main()          /* standard start of any NEMO program */
{
    bool verbose;                  /* declaration of local var. */

#if defined(__cplusplus)
    printf("C++ main: ");
#else
    printf("C main: ");
#endif
    verbose = getbparam("verbose");         /* get that keyword */
    printf("Hello NEMO!\n");                /* do some work ... */
    if (verbose)                            /* and perhaps more */
        printf("Bye then.\n");
}
