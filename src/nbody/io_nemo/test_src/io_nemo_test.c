/* -------------------------------------------------------------- *\ 
|* $Id$
|*
|* NEMO program test using 'io_nemo()' function
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>

#include "io_nemo.h"

string defv[]=  
{ "in=???\n       input snapshot"      ,
  "out=???\n      output snapshot"     ,
  "select=all\n   selected particles"  ,
  "times=all\n    selected time"       ,
  "VERSION=1.2\n  19-Jul-02 jcl"       ,
  NULL
};

#ifndef SINGLEPREC
   #define FTYPE "double"
#else
   #define FTYPE "float"
#endif
string usage="NEMO program test using 'io_nemo()' function in "FTYPE" precision";
/* -------------------------------------------------------------- *\ 
|* Main program
\* -------------------------------------------------------------- */
nemo_main()
{
  real * pos=NULL, * vel=NULL, * mass=NULL;
  real * time=NULL;
  real * acc=NULL;
  real * pot=NULL;
  int  * keys=NULL;
  int   * nbody=NULL;
  string in,out,selt,selp;
  int i,j;

  in  =  getparam("in");
  out =  getparam("out");
  selt = getparam("times");
  selp = getparam("select");

  i = 1;
  while (i!=0) {
    /* read the NEMO snapshot */
    i=io_nemo(in,"n,t,x,v,a,p,m,k,read,info,st,sp,"FTYPE,
	      &nbody,&time,&pos,&vel,&acc,&pot,&mass,&keys,selt,selp);
    
    fprintf(stderr,"Nbody = %d\n",*nbody);
    
    /* save the NEMO snapshot */
    if (i != 0)
      j=io_nemo(out,"n,t,x,v,a,p,m,k,h,save,info,"FTYPE,
		&nbody,&time,&pos,&vel,&acc,&pot,&mass,&keys,in);
    
  }
  /* close the output NEMO snapshot */
  io_nemo(out,"close");
  
  return 1;
}

/* End of [io_nemo_test.c] */
