/*
 * YAPP: Yet Another Plotting Package.
 *
 *	This file includes one of the available YAPP's to go into
 *	the NEMO library. All other ones should be compiled
 *	with the installers discretion, and put into a separate .o
 *	files in $NEMOLIB. When special versions need to be
 *	linked, the underscrore appendix name should be used,
 *	e.g.
 *		yapp_ps.o 
 *	gets linked in as:
 *
 *	cc -o prg_ps prg.c $YAPP_PS -lnemo $YAPPLIB_PS 
 *
 * 	whereas the standard compiled version would be
 *
 *	cc -o prg prg.c $NEMOLIB/libnemo.a $YAPPLIB [-lm -ldl]
 *
 */
#include <stdinc.h>
#include <yapp.h>

#if defined(YAPP_CG)
#include "yapp_cg.c"            /* Sun Core: Color Graphics */

#elif defined(YAPP_CGI)
#include "yapp_cgi.c"           /* CGI (something obsolete?) */

#elif defined(YAPP_CORE)
#include "yapp_core.c"          /* Sun Core (becoming obsolete - see gks) */

#elif defined(YAPP_MONGO)
#include "yapp_mongo.c"         /* Tonry's Mongo (needs f-c interface) */

#elif defined(YAPP_CGM)
#include "yapp_cgm.c"           /* Computer Graphics Metafile (CGM) */

#elif defined(YAPP_GKS)
#include "yapp_gks.c"           /* GKS, for what it's worth */

#elif defined(YAPP_GSS)
#include "yapp_gss.c"           /* GSS, the UNIXPC version of GKS */

#elif defined(YAPP_PGPLOT)
#include "yapp_pgplot.c"        /* Caltech's pgplot (needs f-c interface) */

#elif defined(YAPP_PLPLOT)
#include "yapp_plplot.c"        /* U Texas plplot (a pgplot lookalike in C) */

#elif defined(YAPP_PS)
#include "yapp_ps.c"            /* PostScript */

#elif defined(YAPP_PSG)
#include "yapp_psg.c"            /* PostScript w/ simple greyscale */

#elif defined(YAPP_SM)
#include "yapp_sm.c"            /* Lupton & Monger's Super Mongo */

#elif defined(YAPP_SV)
#include "yapp_sv.c"            /* silly little sunview for movies mainly */

#elif defined(YAPP_SUNTOOLS)
#include "yapp_suntools.c"      /* Sun Core - with various extras */

#elif defined(YAPP_SUNVIEW)
#include "yapp_sunview.c"       /* proper sunview interface, uses mouse */

#elif defined(YAPP_X11)
#include "yapp_x11.c"           /* MIT's X11, R4 */

#elif defined(YAPP_GL)
#include "yapp_gl.c"            /* GL library - using VOGL */

#elif defined(YAPP_ZZZ)
#include "yapp_zzz.c"           /* free slot for Makefile */

#else
#include "yapp_null.c"          /* gets you going on any machine */
#endif

/**************************************************************************/
#if defined(TESTBED)
# include "testyapp.c"
#endif
/**************************************************************************/

