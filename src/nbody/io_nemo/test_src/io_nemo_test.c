/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2005                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* NEMO program test using 'io_nemo()' function                     
+----------------------------------------------------------------- */

#include <stdinc.h>
#include <getparam.h>
#include <snapshot/snapshot.h>
#include "io_nemo.h"

string defv[]=  
{ "in=???\n       input snapshot"      ,
  "out=???\n      output snapshot"     ,
  "select=all\n   selected particles"  ,
  "times=all\n    selected time"       ,
  "VERSION=2.0\n  08-Mar-05 jcl"       ,
  NULL
};

#ifndef SINGLEPREC
   #define FTYPE "double"
#else
   #define FTYPE "float"
#endif
string usage="NEMO program test using 'io_nemo()' function in "FTYPE" precision";
/* ----------------------------------------------------------------
|  Main program                                                    
+---------------------------------------------------------------- */ 
int nemo_main()
{
  real * pos=NULL, * vel=NULL, * mass=NULL;
  real * time=NULL;
  real * acc=NULL;
  real * pot=NULL;
  int  * keys=NULL;
  int   * nbody=NULL, *bits=NULL;
  string in,out,selt,selp;
  int i,j;

  in  =  getparam("in");
  out =  getparam("out");
  selt = getparam("times");
  selp = getparam("select");

  i = 1;
  while (i!=0) {
    /* ********************** */
    /* read the NEMO snapshot */
    /* ********************** */
    i=io_nemo(in,"n,t,x,v,a,p,m,k,b,read,info,st,sp,"FTYPE,
	      &nbody,&time,&pos,&vel,&acc,&pot,&mass,&keys,&bits,selt,selp);
    fprintf(stderr,"Nbody = %d\n",*nbody);


    /* HOW to use nemo bits control */
    if ( ! (*bits & TimeBit) ) {  /* ==> if true, means Time component does not exist */
      fprintf(stderr,"No <Time> component on snapshot.\n");
    } else {
      fprintf(stderr,"Time = %f\n",*time);
    }
    if ( ! (*bits & PosBit) ) {
      fprintf(stderr,"No <Positions> component on snapshot.\n");
    } 
    if ( ! (*bits & VelBit) ) {
      fprintf(stderr,"No <Velocities> component on snapshot.\n");
    } 
    if ( ! (*bits & MassBit) ) {
      fprintf(stderr,"No <Masses> component on snapshot.\n");
    } 
    if ( ! (*bits & AccelerationBit) ) {
      fprintf(stderr,"No <Accelerations> component on snapshot.\n");
    } 
    if ( ! (*bits & PotentialBit) ) {
      fprintf(stderr,"No <Potentials> component on snapshot.\n");
    } 
    if ( ! (*bits & KeyBit) ) {
      fprintf(stderr,"No <Keys> component on snapshot.\n");
    } 

    /* ********************** */
    /* save the NEMO snapshot */
    /* ********************** */
    if (i != 0) {
      /* to save the snapshot here, we still use the bit parameters 'b'.
       * It means that only components requested to save AND components 
       * which have their bit control positionned, during the previously
       * read operation (see above), will be saved.                     
       */
      j=io_nemo(out,"n,t,x,v,a,p,m,k,b,h,save,info,"FTYPE,
		&nbody,&time,&pos,&vel,&acc,&pot,&mass,&keys,&bits,in);
    }
  }
  /* close the output NEMO snapshot */
  io_nemo(out,"close");
  io_nemo(in,"close");
  
  return 1;
}
/* ----------------------------------------------------------------
|  End of [io_nemo_test.c]                                         
+---------------------------------------------------------------- */ 

