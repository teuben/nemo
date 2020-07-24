/*
 *     Benchmark some common table I/O operations
 *
 *      3-jul-2020  V0.1    drafted
 */

//1    nbody=1000000      (for 10M there is a bug)
//2    mkplummer - $nbody | snapprint - format=%20.16f    > p1M.tab
//3    /usr/bin/time tabbench1 p1M.tab .
//4    /usr/bin/time tabtranspose p1M.tab p1Mt.tab $nbody
//5    /usr/bin/time tabbench1 p1Mt.tab .


#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "in=???\n	       input file",
    "out=???\n         output file",
    "nmax=10000\n      Default max allocation",
    "mode=1\n          Benchmark mode",
    "VERSION=0.1\n     3-jul-2020 PJT",
    NULL,
};

string usage="table I/O benchmark";

#ifndef MAX_LINELEN 
#define MAX_LINELEN  2048
#endif

#define MAXPAR 16

void nemo_main(void)
{
    stream istr, ostr;
    char line[MAX_LINELEN];
    int nmax,  *select = NULL;
    int nout, next = 0;
    int    i, j, npar, one = 1;
    string iname = getparam("in");
    real par[MAXPAR], retval, errval, sum;
    int mode = getiparam("mode");
    size_t linelen;
    //string *lineptr = line;
    char **lineptr = line;    

    npar = inifie("sqrt(%1*%1+%2*%2+%3*%3)");
    dprintf(0,"MAX_LINELEN=%d\n",MAX_LINELEN);

    istr = stropen(getparam("in"),"r");
    ostr = stropen(getparam("out"),"w");

    i = 0;
    sum = 0.0;
    linelen = 0;
    if (mode < 0) {
      while (fgets(line,MAX_LINELEN,istr) != NULL) {
	i++;
      }
    } else if (mode == 0) {
	while (fgets(line,MAX_LINELEN,istr) != NULL) {
	  i++;
	  npar = nemoinpr(line,par,MAXPAR);	
	}
    } else if (mode == 1) {
	while (fgets(line,MAX_LINELEN,istr) != NULL) {
	  i++;
	  npar = nemoinpr(line,par,MAXPAR);	
	  retval = sqrt(par[0]*par[0] + par[1]*par[1] + par[2]*par[2]);
	  sum += retval;	  
	}
    } else if (mode == 2) {
	while (fgets(line,MAX_LINELEN,istr) != NULL) {
	  i++;
	  npar = nemoinpr(line,par,MAXPAR);	
	  dofie(par,&one,&retval,&errval);
	  sum += retval;	  
	}
    }
	  
	
    //while (get_line(istr,line)) {
    //while (getline(&lineptr, &linelen, istr)) {
      //dprintf(1,"%ld \n",linelen);
	  //fprintf(ostr,"%s %g %g %g  %g\n",line,par[0],par[1],par[2],retval);
	  //fputs(line,ostr);

    dprintf(0,"sum=%g\n",retval);
    strclose(istr);
    strclose(ostr);
    dprintf(0,"Read %d lines\n",i);
}

// bench:
// tabgen tab1 1000000   3
// tabgen tab2 10000000  3
// tabgen tab3 100000000 3
//
// fgets method on tab3:
// -1               1.9  1.9  2.0
//  0 24.0 24.3    24.0 27.6 25.7 26.3 27.4 26.9
//  1 25.9 26.4    28.9 27.7 28.1 26.1 27.5 28.0
//  2 28.0         32.4 32.7 33.2 31.6 28.3 28.1
