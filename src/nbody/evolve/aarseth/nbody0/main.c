/*
 *  MAIN.C : NEMO driver for Aarseth's fortran nbody0 program
 *
 *	 3-jul-89	V1.1  
 *	23-jan-90	V1.2  ???  		PJT
 *	17-may-91	V1.2a bit more doc	PJT
 *	 6-mar-92       Happy GCC2.0		PJT
 *       6-aug-97       V1.3 options=           PJT
 *	11-feb-98	V1.3a fixed that bad index bug in nbody0.f/c	PJT
 *	 6-jan-00       V1.3b changed to F77_FUNC macros	PJT
 *	21-jan-00	V1.4 added reset=t|f	PJT
 *      21-feb-04       V1.5 protection for f2dot=0   PJT
 *      24-feb-04       V2.0 added f3dot=   ** only for nbody00 now ** PJT
 */

#include <stdinc.h>

string defv[] = {
    "in=???\n       Input (snapshot) file",
    "out=???\n      Output (snapshot) file",
    "eta=0.02\n     Accuracy parameter - determines timesteps",
    "deltat=0.25\n  When to dump major output",
    "tcrit=2\n      When to stop integrating",
    "eps=0.05\n     Softening length",
    "reset=t\n      Reset timestep after datadump (debug) (t|f)",
    "f3dot=f\n      Use more advanced timestep determination criterion?",
    "options=\n     Optional output of 'step' into AUX",
    "VERSION=2.0a\n 13-mar-04 PJT",
    NULL,
};

string usage = "NEMO driver for nbody0";

extern void pars_2_aarseth(void); 
extern void call_2_aarseth(void);

/*
 *  nemo_main is the main() type entry for the program. The functions
 *  called here live in nbodyio.c, which serves as a C-to-Fortran 
 *  layer.
 */

void nemo_main(void)
{
    pars_2_aarseth();       /* get parameters */
    call_2_aarseth();       /* run the code */
}

