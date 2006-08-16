/*
 * PATHLABEL.C: label particles according to position along parametric path.
 *
 *	xx-xxx-xx  0.0 -- JEB
 *	20-feb-92  0.1 removed dex() - which is defined in the NEMO kernel PJT
 *	 7-mar-92  0.2 fixed pathw string bug  & happy gcc2.0	           pjt
 *      15-aug-06  0.3 
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <filefn.h>
#include <loadobj.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

#include <sys/types.h>
#include <unistd.h>

string defv[] = {
    "in=???\n			Input file",
    "out=???\n			Output file",
    "pathx=sin(TWO_PI*q)\n	X(q) coordinate of Path",
    "pathy=cos(TWO_PI*q)\n	Y(q) coordinate of Path",
    "pathw=0.1\n		--never used--",
    "npnts=128\n		Number of points",
    "VERSION=0.3\n		15-aug-06 PJT",
    NULL,
};

string usage = "label particles according to position along parametric path";

string cvsid="$Id$";



local Body *btab = NULL;
local int nbody;
local real tsnap = 0.0;

local rproc pathx, pathy;
local real pathw;
local int npnts;


local void loadfuns(void);	/* dummy local math loader */
local void pathlabel(void);
local rproc compile_func(string expr, string var);

void nemo_main(void)
{
    stream instr, outstr;
    rproc compile_func();
    int snapbits;

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    pathx = compile_func(getparam("pathx"), "q");
    pathy = compile_func(getparam("pathy"), "q");
    pathw = getdparam("pathw");
    npnts = getiparam("npnts");
    get_snap(instr, &btab, &nbody, &tsnap, &snapbits);
    if (snapbits & (PhaseSpaceBit == 0)) {
      error("no phasespace info");
      loadfuns();
    }
    pathlabel();
    snapbits |= AuxBit;
    put_snap(outstr, &btab, &nbody, &tsnap, &snapbits);

}

/*
 * COMPILE_FUNC: workhorse routine which compiles a real function of
 * a real variable, loads the object module, and returns a pointer.
 */


bool havesyms = FALSE;		/* TRUE if symbols have been loaded */
int funcmpld = 0;	       	/* count of functions compiled */

rproc compile_func(string expr, string var)
{
    char path[64], file[64], func[64], cmmd[128];
    stream cdstr;
    proc result;

    if (! havesyms)
	mysymbols(getargv0());
    havesyms = TRUE;
    fprintf(stderr, "[compile_func: invoking cc]\n");
    sprintf(path, "/tmp/cf_%d", getpid());
    sprintf(file, "%s.c", path);
    cdstr = fopen(file, "w");
    sprintf(func, "_cf_%d", ++funcmpld);
    fprintf(cdstr, "#include <stdinc.h>\n#include <math.h>\n");
    fprintf(cdstr, "real %s(%s)\n", &func[1], var);
    fprintf(cdstr, "real %s;\n", var);
    fprintf(cdstr, "{\n    return (%s);\n}\n", expr);
    fclose(cdstr);
    sprintf(cmmd, "cc -o %s.o -c %s.c", path, path);
    if (system(cmmd) != 0)
	error("compile_func in %s: cant %s\n", getargv0(), cmmd);
    sprintf(file, "%s.o", path);
    loadobj(file);
    sprintf(cmmd, "rm %s.c %s.o", path, path);
    if (system(cmmd) != 0)
	error("compile_func in %s: cant %s\n", getargv0(), cmmd);
    result = findfn(func);
    if (result == NULL)
	error("compile_func in %s: cant find %s\n", getargv0(), func);
    return ((rproc) result);
}

/*
 * LOADFUNS: this curious function is never called.  However, its mere
 * existence is sufficient to fool ld into loading the math functions.
 */

local void loadfuns(void)
{
#if 0
    real cbrt(), sqrt(), qbe(), sqr();
    real sin(), cos(), asin(), acos();
    real tan(), atan(), atan2();
    real exp(), dex(), log(), log10(), pow();
    real fabs(), floor(), ceil(), rint();
#endif

    (void) cbrt(1.0);
    (void) sqrt(1.0);
    (void) qbe(1.0);
    (void) sqr(1.0);
    (void) sin(1.0);
    (void) cos(1.0);
    (void) asin(1.0);
    (void) acos(1.0);
    (void) tan(1.0);
    (void) atan(1.0);
    (void) atan2(1.0, 1.0);
    (void) exp(1.0);
    (void) dex(1.0);
    (void) log(1.0);
    (void) log10(1.0);
    (void) pow(1.0, 1.0);
    (void) fabs(1.0);
    (void) floor(1.0);
    (void) ceil(1.0);
    (void) rint(1.0);
}

local void pathlabel(void)
{
    int j, i;
    real q, px, py, dx, dy, dij;

    for (j = 0; j <= npnts; j++) {
	q = ((real) j) / ((real) npnts);
	px = (pathx)(q);
	py = (pathy)(q);
	if (j%8==0) dprintf(1,"q, px, py = %f, %f, %f\n", q, px, py);
	for (i = 0; i < nbody; i++) {
	    dx = Pos(&btab[i])[0] - px;
	    dy = Pos(&btab[i])[1] - py;
	    dij = sqrt(dx*dx + dy*dy);
	    if (j == 0 || dij < Aux(&btab[i]))
		Aux(&btab[i]) = dij;
	}
    }
}
