/*
 *	Pattern: extract pattern speed as a user parameters
 *		from the command line:
 *
 *	Since we've been stupid enough to design the pattern
 *	into the first parameter of the POTPARS parameter,
 *	a simple getdparam("potpars") will extract it correctly,
 *	but not withouth a warning message from nemoinpd()
 *	This wrapper will prevent this.
 *	It's a kludge until better potential interface is designed.
 *
 *	 9-Jun-92	Created		PJT
 *	22-feb-94   documented that since oct-93 this routine is obsolete
 *		    however the get[dsub]param may be interesting enough
 *		    to save		PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>

extern string  *burststring(string, string);
extern void     freestrings string *);

real getdsubparam(string key, int n, real def)
{
    string *vp;
    int i;
    double val;

    n--;        /* input is 1 based, and local arrays are 0 based */
    vp = burststring(getparam(key)," ,");
    i = xstrlen(vp,sizeof(string))-1;
    if (n < i ) {
        if (nemoinpd(vp[n],&val,1) != 1) {
            warning("getdsubparam: error parsing #%d: \"%s\"",n+1,vp[n]);
            return def;
        } 
    } else {
       val = def;
    }
    freestrings(vp);
    return val;
}

#if defined(TESTBED)

string defv[] = {
    "pars=0,1,2,3,4,5,6,7:10\n      Some legal and illegal values",
    "VERSION=1.1\n      22-feb-94 PJT",
    NULL,
};

nemo_main()
{
    int i;

    if (!hasvalue("pars")) {
        putparam("pars","0,1,2,3,4,5,6,7:10");
        warning("This run should produce a warning");
    }
    for (i=1; i<10; i++)
        printf("%d -> %g\n",i,getdsubparam("pars",i,-1.0*i));
}
#endif
