/*
 * arbitrary function parser using loadobj/dlopen
 *	17-sep-90 created				pjt
 *	   apr-92 .. updated
 *	17-jan-95 fixed real=float problems		pjt
 *       13-dec-95 using new ldso command               pjt
 *	6-jun-96   fixed solaris search path problem	pjt   (D-day)
 *     16-sep-00   changed from TESTBED to TOOLBOX      pjt
 *      1-apr-01   new '.so' based for NEMO V3		pjt
 *     20-jun-01   fix header				pjt
 *     29-aug-06   use real_proc now                    pjt
 *
 *   TODO:     don't use K&R style definition in the produced code,  but ANSI
 */


#include <stdinc.h>
#include <loadobj.h>
#include <getparam.h>

extern string *burststring(string, string);

#define MAXNAMELEN  32

/*
 *  getrfunc:   
 *             defin=name(var1,var2,...)
 *             expr=f(var1,var2,...,par1,par2,..)   one line or multi if '{'
 *             pars=pas1,val1,par2,val2,...
 *
 *
 *  returns pointer to a safe string containing pointer to
 *  the real function. For safety the expewcted number of dependant
 *  variables in the argument list of the function is given
 *  too. If this is set to -1, no check is made however.
 */

real_proc getrfunc(string defin, string expr, string pars, int *nexpvar)
{
    stream fstr;
    char *cp, name[MAXNAMELEN+1];
    string *parts, np, vp;
    int i, npar=0, nvar=0, namelen = 0;
    double dvar;
    proc result;

    fstr=stropen("tmp.c","w!");
    
    fprintf(fstr,"/* ROUTINE CREATED BY GETRFUNC - DO NOT EDIT */\n");
    fprintf(fstr,"#include <stdinc.h>\n\n");
    
    fprintf(fstr,"\n\nreal %s\n",defin);


    cp=defin;
    while (*cp != 0 && *cp != '(')           /* make name */
        if (namelen < MAXNAMELEN)
            name[namelen++] = *cp++;
        else
            cp++;
    name[namelen] = '\0';
    if (namelen==MAXNAMELEN) 
        dprintf(0,"Warning: name too long in: %s\n",expr);
    if (*cp != '(') {
        dprintf(0,"setfunc_name: no arguments supplied (%s)\n",expr);
    } else
        cp++;

    parts = burststring(cp,"), ");          /* extract argument variables */
    for (nvar=0; parts[nvar] != NULL; i++) {
        fprintf(fstr,"real %s; /* var # %d */\n",parts[nvar],nvar+1);
        nvar++;
    }
    if (*nexpvar >= 0) {            /* check number of found variables */
        if (nvar != *nexpvar)
            error("getrfunc: Number of variables %d should have been %d",
                        nvar, *nexpvar);
    } else
        *nexpvar = nvar;
    fprintf(fstr,"{\n");

    /* free(parts) */

    parts = burststring(pars,", ");             /* Get fixed parameters */
    for (i=0; parts[i] != NULL;) {
        np = parts[i];
        vp = parts[i+1];
        if (vp==NULL) error("Parameter %s not given a value",np);
	if (nemoinpd(vp,&dvar,1)!=1) error("Error parsing %s",vp);
        fprintf(fstr,"real %s = %f; /* par # %d */\n",
                np,dvar,npar+1);
        i += 2;
        npar++;
    }
    /* free(parts) */

    cp = expr;
    if (*cp == '{')                 /* also allow more complicated ?? '}'*/
        fprintf(fstr,"%s\n",cp+1);
    else
        fprintf(fstr,"\n   return(%s);\n}\n",expr);
    strclose(fstr);

    dprintf(0,"[Compiling %s=%s]\n",defin,expr);	    /* Compile and load the function */
#if defined(LOADOBJ3)
    if (system("make -f $NEMOLIB/Makefile.lib tmp.so") != 0) 
        error("Cannot make with tmp.so");
    loadobj("./tmp.so");
#else
    if (system("cc -g -c tmp.c") != 0) 
        error("Cannot compile with cc");
    if (system("ldso tmp"))
        error("Cannot create sharable file");
    loadobj("./tmp.o");
#endif
    mapsys(name);
    result = findfn(name);

    if (result==NULL) error("getrfunc: Cannot find %s\n",name);

    return (real_proc) result;
}


#ifdef TOOLBOX

#include <getparam.h>

string defv[] = {
    "def=f(x)\n         Function name, and independant variables list",
    "expr=x+a\n         Function expression",
    "pars=a,1\n         List of fixed parameters to the function",
    "x=0:2:1\n          X-Values to test for...",
    "y=0:2:1\n         	Y-Values to test for...",
    "z=0:2:1\n         	Z-Values to test for...",
    "ndim=1\n           Dimensionality of output 1=x 2=x,y 3=x,y,z",
    "format=%g\n        Output format",
    "VERSION=3.1\n      29-aug-06 PJT",
    NULL,
};

string usage="Arbitrary Expression Parser tabulator";

#ifndef MAXN
#define MAXN  16384
#endif

nemo_main()
{
    int i, j, k, ndim, nvar, nx, ny, nz;
    real x[MAXN], y[MAXN], z[MAXN];
    char fmt[80];
    string def, fmt1;
    real_proc func;

    mysymbols(getargv0());

    ndim = getiparam("ndim");
    fmt1 = getparam("format");
    strcpy(fmt,fmt1);
    for (i=0; i<ndim; i++) {
        strcat(fmt," ");
        strcat(fmt,fmt1);
    }
    strcat(fmt,"\n");
    nvar = ndim;
    if (ndim<1 || ndim>3) error("Only ndim=1,2,3 are supported");
    def = getparam("def");
    func = getrfunc(getparam("def"),getparam("expr"),getparam("pars"),&nvar);

    nx = nemoinpr(getparam("x"),x,MAXN);
    ny = nemoinpr(getparam("y"),y,MAXN);
    nz = nemoinpr(getparam("z"),z,MAXN);
    dprintf(0,"Ndim=%d Nvar=%d Nxyz=(%d,%d,%d)\n",ndim,nvar,nx,ny,nz);

    switch(nvar) {
       case 1:
        for(i=0; i<nx; i++) {
           printf(fmt, x[i], (*func)(x[i]) );
	}
        break;
      case 2:
      case 3:
      default:
        error("Illegal case or un-implemented nvar = %d",nvar);
    }
}
#endif
