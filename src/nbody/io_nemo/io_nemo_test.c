/* -------------------------------------------------------------- *\ 
|* NEMO program test using 'io_nemo()' function
|*    16-April-97 jcl	original version
|*    16-jul-01   pjt   added maxio.h
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>

#include "maxio.h"
#include "io_nemo.h"

string defv[]=  
{ "in=???\n       input snapshot",
  "out=???\n      output snapshot",
  "select=all\n   selected particles",
  "times=all\n    selected time",
  "VERSION=1.0a\n 16-jul-01 PJT",
  NULL
};

string usage="NEMO program test using 'io_nemo()' function";
/* -------------------------------------------------------------- *\ 
|* Main program
\* -------------------------------------------------------------- */
nemo_main()
{
  double * pos=NULL, * vel=NULL, * mass=NULL;
  double * time=NULL;
  double * acc=NULL;
  double * pot=NULL;
  int   * nbody=NULL;
  string in,out,selt,selp;
  int i,j;

  in  =  getparam("in");
  out =  getparam("out");
  selt = getparam("times");
  selp = getparam("select");

  i = 1;
  while (i!=0)
    {   
      /* read the NEMO snapshot */
      i=io_nemo(in,"double,n,t,x,v,a,p,m,read,info,st,sp",
		&nbody,&time,&pos,&vel,&acc,&pot,&mass,selt,selp);

      fprintf(stderr,"Nbody = %d\n",*nbody);

      /* save the NEMO snapshot */
      if (i != 0)
	j=io_nemo(out,"double,n,t,x,v,a,p,m,h,save,info",
		  &nbody,&time,&pos,&vel,&acc,&pot,&mass,in);

    }
  /* close the output NEMO snapshot */
  io_nemo(out,"close");
}



