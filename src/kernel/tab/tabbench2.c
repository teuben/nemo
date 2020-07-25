/*
 *     Benchmark some common table I/O operations
 *
 *     -  2" to read all lines using getline()  [100M benchmark see below]
 *     - 22" to break them using nemoinpr()
 *     -  6" to use dofie()
 *     - for same size file, longer lines (thus fewer to read) is more efficient
 *
 *     24-jul-2020  V0.1    drafted
 *     25-jul-2020  V0.2    cleaned up
 */


#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "in=???\n	       input file",
    "out=???\n         output file",
    "mode=1\n          Benchmark mode",
    "nmax=10000\n      Default max allocation",
    "VERSION=0.3\n     25-jul-2020 PJT",
    NULL,
};

string usage="table I/O benchmark";

#ifndef MAX_LINELEN 
#define MAX_LINELEN  2048
#endif

#define MAXPAR 16

extern int inifie(string);
extern void dofie(real *, int *, real *, real *);

void nemo_main(void)
{
    stream istr, ostr;
    int nmax,  *select = NULL;
    int nout, next = 0;
    int    npar, one = 1;
    size_t nlines;
    string iname = getparam("in");
    real par[MAXPAR], retval, errval, sum;
    int mode = getiparam("mode");
    size_t linelen = MAX_LINELEN;
    //char line[MAX_LINELEN];
    char *line = malloc((linelen) * sizeof(char));        

    npar = inifie("sqrt(%1*%1+%2*%2+%3*%3)");
    dprintf(0,"MAX_LINELEN=%d\n",MAX_LINELEN);

    istr = stropen(getparam("in"),"r");
    ostr = stropen(getparam("out"),"w");

    nlines = 0;
    sum = 0.0;
    linelen = 0;
    if (mode <= 0) {
      if (mode < 0)
	dprintf(0,"Just reading\n");
      else
	dprintf(0,"Just reading and writing\n");
      //while (fgets(line,MAX_LINELEN,istr) != NULL)
      while (getline(&line, &linelen, istr) != -1) {
	nlines++;
	if (mode == 0) fputs(line,ostr);
      }
    } else if (mode == 1) {
        dprintf(0,"nemoinp to split line\n");      
        while (getline(&line, &linelen, istr) != -1) {
	  nlines++;
	  npar = nemoinpr(line,par,MAXPAR);	
	}
    } else if (mode == 2) {
        dprintf(0,"nemoinp + sqrt()\n");            
        while (getline(&line, &linelen, istr) != -1) {	  
	  nlines++;
	  npar = nemoinpr(line,par,MAXPAR);	
	  retval = sqrt(par[0]*par[0] + par[1]*par[1] + par[2]*par[2]);
	  sum += retval;	  
	}
    } else if (mode == 3) {
        dprintf(0,"nemoinp + fie(sqrt())\n");                  
        while (getline(&line, &linelen, istr) != -1) {	        
	  nlines++;
	  npar = nemoinpr(line,par,MAXPAR);	
	  dofie(par,&one,&retval,&errval);
	  sum += retval;	  
	}
    }

    strclose(istr);
    strclose(ostr);
    
    dprintf(0,"sum=%g\n",sum);
    dprintf(0,"Read %ld lines\n",nlines);
}

// bench:
// tabgen tab1 1000000   3
// tabgen tab2 10000000  3
// tabgen tab3 100000000 3
// tabgen tab4 3 100000000
//
// fgets method on tab3:
// -1               1.9  1.9  2.0  
//  0               3.2  3.1  3.1
//  1 24.0 24.3    24.0 27.6 25.7 26.3 27.4 26.9 24.4 25.6 28.8 28.6
//  2 25.9 26.4    28.9 27.7 28.1 26.1 27.5 28.0 27.7 28.5 28.6 
//  3 28.0         32.4 32.7 33.2 31.6 28.3 28.1 33.4 33.6 34.4

// e.g. tabgen tab3    80 sec ->  34 MB/sec ( low level I/O :   435 MB/sec )

// T480  2.8 4.8 34.5 35.4 42.9
// X1Y4  1.9 3.1 28.6 28.6 34.4 (about 20% more for fie vs. compiled)
