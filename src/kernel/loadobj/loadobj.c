/* loadobj.c:  
 *	Merely includes the right loadobjXXX.c for a particular machine.
 *	The Makefile defaults compilation with -D$(MACH) 
 *	The first few #ifdef's below are ''upper'' case, and can be used
 *	on ''any'' machine, typically using a -D... compile switch,
 *      If those fail in the test, it will try one of the
 *	machine supplied triggers, the order of e.g. sun after sparc is
 *	in this respect important. Then, if those fail too, the
 *	STUB loader is used. This one can also be loaded
 *	specifically with -DSTUB.
 *
 *	At the end also includes the testloadobj.c code for TESTBED
 *	testing (-DTESTBED)
 *
 *	28-oct-90	finalized for new Nemo				PJT
 *	18-nov-90	added GNU dld					PJT
 *	25-may-91 	added the STUB option - more doc		PJT 
 *       5-nov-91       added NEXT option                               PJT
 *      24-nov-91       added AIX for R6000 (IBM)			PJT
 *	28-nov-91	formally added MIPS (note NEXT/AIX/MIPS all bad)PJT
 *	20-jan-94	NULL is now STUB				pjt
 *	22-oct-96	DL is default for linux				pjt
 *      12-jul-03       DL is default for darwin (though needs dlcompat) PJT
 *      04-Jul-07       DL for OpenBSD                             Woodchuck
 */

#include <stdinc.h>		/* standard Nemo include file */


/* patch up for some predefined defaults we understand: */

#if !defined(AUTO)
# if defined(sun) && defined(sparc) && defined(SYSV)
#  define DL 1
# elif defined(sgi) || defined(alpha) || defined(__alpha) 
#  define DL 1
# elif defined(linux) && defined(SYSV) || defined(darwin)
#  define DL 1
# elif defined(__OpenBSD__)
#  define DL 1
# endif
#endif

#if defined(DL)
char *loadobj_version="loadobjDL";
#include "loadobjDL.c"		/* new standard dlopen(3x) functions */

#elif defined(DLD)
char *loadobj_version="loadobjDLD";
#include "loadobjDLD.c"		/* GNU dld library */

#elif defined(XDL)
char *loadobj_version="loadobjXDL";
#include "loadobjXDL.c"		/* xdl (sun/dec) library */

#elif defined(STUB) || defined(sgi)
char *loadobj_version="loadobjSTUB";
#include "loadobjSTUB.c"	/* stubbed; gets you going on an odd machine */

#elif defined(sparc) && !defined(SYSV)
char *loadobj_version="loadobjSPARC";
#include "loadobjSPARC.c"	/* sun4 SPARC (works 99% of the time) */

#elif defined(sun) && !defined(SYSV)
char *loadobj_version="loadobjSUN";
#include "loadobjSUN.c"		/* general bsd (sun3), multiflow (?) */

#elif defined(_trace_)
char *loadobj_version="loadobj_mf";
#include "loadobj_mf.c"		/* MF hacking: does not work*/

#elif defined(COFF)
char *loadobj_version="loadobjCOFF";
#include "loadobjCOFF.c"	/* general sys5, not always works?? */

#elif defined(alliant)
char *loadobj_version="loadobjFX";
#include "loadobjFX.c"		/* alliant: does not work */

#elif defined(mc68k)
char *loadobj_version="loadobj3B1";
#include "loadobj3B1.c"		/* 3b1: doesn't always work - like COFF */

#elif defined(NeXT) || defined(NEXT)
char *loadobj_version="loadobjNEXT";
#include "loadobjNEXT.c"	/* NeXT: MACH O - doesn't work yet either */

#elif defined(aix)
char *loadobj_version="loadobjAIX";
#include "loadobjAIX.c"		/* IBM Risc - testing nov 91 */

#elif defined(MIPS) || defined(mips)
char *loadobj_version="loadobjMIPS";
#include "loadobjMIPS.c"	/* DEC5000 - testing nov 91 */

#elif defined(linux)
char *loadobj_version="loadobjLINUX";
#include "loadobjLINUX.c"		/* Linux - testing */

#elif defined(__CYGWIN__)
char *loadobj_version="loadobjCYGWIN";
#include "loadobjCYGWIN.c"		/* Cygwin  - testing by JM*/

#else
				/* if all else failed....*/
char *loadobj_version="loadobjSTUB";
#include "loadobjSTUB.c"	/* stubbed; gets you going on an odd machine */

#endif

#if defined(TESTBED)
# include "testloadobj.c"
#endif
