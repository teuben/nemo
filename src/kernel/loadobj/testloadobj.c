/*
 *  TESTLOADOBJ:  TESTBED code for any of the loadobjXXX.c sets
 *
 *	This file is normally #included by loadobj.c, which also
 *	includes the correct loadobjXXX.c machine dependant file.
 *
 *	10-jan-95	ansi-fied code, also showed that one cannot
 *			the order in which function arguments are
 *			'executed' (on Linux the test 'failed',
 *			nmain was always lagging one iteration)
 *	15-dec-95 	use ldso
 *      27-mar-97       1.2 cleaned up old -DPIC code that we don't use anymore
 *	 1-apr-01       1.3 new .so based objects for NEMO V3
 */

#include <stdinc.h>
#include <getparam.h>
#include <loadobj.h>

string defv[] = {
    "func=(n==0?1:n*f(n-1))\n   Function to test, or to call",
    "low=1\n                    Low loop index for test run",
    "high=10\n                  High loop index for test run",
    "object=\n                  Name of object file to load & test",
    "type=\n                    Type of argument {r,i,v}",
    "cleanup=t\n                Cleanup local tmp files?",
    "VERSION=1.3\n		1-apr-01 pjt",
     NULL,
};

string usage="testbed for the loadobj package";

int nmain = 0;			/* !! a global symbol !! */

void nemo_main(void)
{
    string fdef, object, type;
    int    low, high, i, j, n;
    stream tmpfile;
    bool   Qcleanup = getbparam("cleanup");
    real   x;
    proc   fn;
    iproc  fni;
    rproc  fnr;
    extern char *loadobj_version;

    dprintf(0,"LoadobjTest: version/type = %s\n",loadobj_version);
    fdef = getparam("func");
    low = getiparam("low");
    high = getiparam("high");
    object = getparam("object");
    type = getparam("type");

    mysymbols(getparam("argv0"));

    if (*object) {
        dprintf(0,"Attempting to load %s\n", object);
        loadobj(object);         /* load it */
        if (*fdef) {             /* if a function name supplied ... */
            switch (*type) {     /* determine it's type */
              case 'r':
                    fnr = (rproc) findfn(fdef);  /* find name */
                    if (fnr==NULL) error("Error finding %s",fdef);
                    for (i=low; i<=high; i++) {
                        x = (*fnr)((real)i);          /* and call it */
                        dprintf(0,"%d ==> %f\n",i,x);
                    }
                    break;
              case 'i':
                    fni = (iproc) findfn(fdef);
                    if (fni==NULL) error("Error finding %s",fdef);
                    for (i=low; i<=high; i++) {
                        n = (*fni)(i);
                        dprintf(0,"%d ==> %d\n",i,n);
                    }
                    break;
              default:
                    fn = (proc) findfn(fdef);
                    if (fn==NULL) error("Error finding %s",fdef);
                    for (i=low; i<=high; i++)
                        (*fn)();
                    break;
            } /* switch */
        } /* if(func) */
        exit(0);
    } /* if(object) */

    if (*fdef != '0') {
        tmpfile = stropen("ld-tmp.c", "w!");             /****** TMPFILE ******/
        fprintf(tmpfile, "#include <stdio.h>\n");
        fprintf(tmpfile, "double sin(), cos(), tryext();\n");
        fprintf(tmpfile, "extern int nmain;\n");
        fprintf(tmpfile, "static int nlocal=0;\n");
        fprintf(tmpfile, "int f(n)\n");
        fprintf(tmpfile, "int n;\n");
        fprintf(tmpfile, "{\n");
        fprintf(tmpfile, "    dprintf(1,\" calling f(n=%%d)\\n\",n);\n");
        fprintf(tmpfile, "    nmain += n;\n");
        fprintf(tmpfile, "    tryext((double)n);\n");	
        fprintf(tmpfile, "    return (%s);\n", fdef);
        fprintf(tmpfile, "}\n");
        fclose(tmpfile);                                /**********************/
#if defined(LOADOBJ3)
        if (system("make -f $NEMOLIB/Makefile.lib ld-tmp.so") != 0)             
	    error("function %s does not parse\n", fdef);
#else
        if (system("cc -c ld-tmp.c; ldso ld-tmp") != 0)             
	    error("function %s does not parse\n", fdef);
#endif
    } else
	dprintf(0,"Attempt to load local file ld-tmp.o\n");
#if defined(LOADOBJ3)
    loadobj("./ld-tmp.so"); /* extra ./ is needed for Solaris/LD_LIBRARY_PATH */
#else
    loadobj("./ld-tmp.o"); /* extra ./ is needed for Solaris/LD_LIBRARY_PATH */
#endif
    if (Qcleanup)
        if (system("rm -f ld-tmp.*") != 0)
	    error("cannot rm ld-tmp.*\n");
    fni = (iproc) findfn("f");	   /* Try mapsys_SYSV */
    if (fni == NULL) {
	fni = (iproc) findfn("_f");   /* Try mapsys_BSD */
	if (fni == NULL)
		error("function not correctly defined");
        warning("mapsys should be defined with -Dmapsys_bsd");
    } else
        warning("mapsys is clean; no -Dmapsys_??? is needed (-Dmapsys_sysv)");
    for (i = low; i <= high; i++) {
	j = (*fni)(i);
	dprintf(0,"f(%d) = %d,  nmain=%d\n", i, j, nmain);
        x++;
    }
    dprintf(0,"###: The last line of the default output should be:\n");
    dprintf(0,"f(10) = 3628800,  nmain=220\n");
}

/* this routine will be called by the loaded object */

double tryext(double x) 
{
  printf("Tryext (%g) called\n",x);
  return (x+1.0);
}

/* ensures loading of some math functions - is never called however */

void math_loader(void)	
{
    (void)sin(1.0);
    (void)cos(1.0);
    (void)sqrt(1.0);
    (void)sqr(1.0);
    (void)sqr(1.0);
    (void)asin(1.0);
}
