/*
 * Test of Nemo's i/o stuff
 *	3-oct-90    Created for new 'micro (nu)nemo'	PJT
 *     12-oct-90    Read or Write big files...          PJT
 *     13-oct-90    Random access                       PJT
 *     27-feb-94    ansi
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>


string defv[] = {
    "in=\n         Input file (if provided)",
    "out=\n        output file (if provided)",
    "real=1.234\n  Real number to write (if out=)",
    "count=1\n     Number of reals to write (if out=)",
    "random=1\n    Number of times to reverse (random) write count",
    "VERSION=1.3\n 31-may-05 PJT",
    NULL,
};

string usage = "testing NEMO's I/O";

string cvsid="$Id$";


extern void get_nanf(float *);
extern void get_nand(double *);

int get_number(string name, real *x)
{
    stream str;
    real *buf;
    int  n;

    str = stropen(name,"r");
    get_history(str);

    if (!get_tag_ok(str,"n"))
        error("%s: missing \"n\", Not a valid TESTIO input file",name);
    get_data(str,"n",IntType,&n,0);
    buf = (real *) allocate(n*sizeof(real));

    if (!get_tag_ok(str,"x"))
        error("%s: missing \"x\", Not a valid TESTIO input file",name);
    get_data(str,"x",RealType,buf,n,0);
    strclose(str);

    *x = *buf;
    dprintf(1,"Read number %f from file %s\n",*x,name);
    free( (char *) buf);
    return n;
}

void put_number(string name, real x, int n, int r)
{
    stream str;
    real *buf;
    int nout,i;

    str = stropen(name,"w");
    put_history(str);
    nout = n*r;
    //put_set(str,"Testing");
    put_data(str,"n",IntType,&nout,0);
    buf = (real *) allocate(n*sizeof(real));
    buf[0] = x;   /* only set first and third number */
    if (n>2) buf[2] = x;
    if (r==1) {
        put_data(str,"x",RealType,buf,0,n);
        put_data(str,"y",RealType,buf,0,n);
    } else {
        put_data_set(str,"x",RealType,nout,0);
        for (i=r;i>0;i--) {
            buf[i] = -1.0 * i;
            if (n>3) get_nand((double*)&buf[3]);
	    if (n>4) get_nanf((float *)&buf[4]);
            put_data_ran(str,"x",buf,(i-1)*n,n);
        }
        put_data_tes(str,"x");

        put_data_set(str,"y",RealType,nout,0);
        for (i=r;i>0;i--) {
            buf[i] = -1.0 * i;
            if (n>3) get_nand((double *)&buf[3]);
	    if (n>4) get_nanf((float  *)&buf[4]);
            put_data_ran(str,"y",buf,(i-1)*n,n);
        }
        put_data_tes(str,"y");

#if 1
    dprintf(0,"nout=%d\n",nout);
    for (i=0; i<nout; i++) buf[i]=10+i;
    put_data_set(str,"z1",RealType,nout,0);

    //    for (i=nout-1; i>=0; i--) put_data_ran(str,"z1",&buf[i],i,1);
    //    for (i=0; i<nout; i++) put_data_ran(str,"z1",&buf[i],i,1);
    put_data_ran(str,"z1",buf,0,1);
    put_data_tes(str,"z1");
#endif
    }
    //put_tes(str,"Testing");
    strclose(str);
    dprintf(1,"Wrote number %f to file %s\n",x,name);
    free( (char *) buf);
}

/* -------------------------------------------------------------------------- */

void nemo_main()
{
    real x,y = 0;
    int count, random;
    string name;
    int work = 0;

    x = getdparam("real");    
    count = getiparam("count");
    random = getiparam("random");
    
    name = getparam("in");
    if (*name) {
      work++;
      get_number(name,&y);
    } else
        dprintf(1,"No input file specified\n");
    
    name = getparam("out");
    if (*name) {
      work++;
      put_number(name,y+x,count,random);
    } else
        dprintf(1,"No output file specified\n");

    if (work==0) warning("No work done, use in=, out=, or both");
}
