/*
 *  SNAPPAT: create a regular Particle-Attribute-Time cube from a snapshot
 *
 *       8-aug-2009    0.1 written                                    PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <image.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (snapshot)",
    "out=???\n                  Output file (image)",
    "options=x,y,z,vx,vy,vz\n	Things to output",
    "times=all\n		Times to select snapshot",
    "ntime=\n                   if used, pre-allocate this number of snapshots",
    "first=f\n                  only write first cube?",
    "VERSION=0.2\n		9-aug-09 PJT",
    NULL,
};

string usage="convert snapshot to a regular Particle-Attribute-Time image";

string cvsid="$Id$";

#define MAXOPT    64



void fixheader(imageptr iptr, string options, real t0, real dt);

void nemo_main()
{
    stream instr, outstr;
    real   tsnap, dr, aux, t0, dt;
    string times;
    Body *btab = NULL, *bp, *bq;
    int i, j, k, n, nbody, bits, nopt, ParticlesBit, ntime;
    char fmt[20],*pfmt;
    string *burststring(), *opt;
    rproc btrtrans(), fopt[MAXOPT];
    imageptr iptr = NULL; 
    bool Qfirst = getbparam("first");

    ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
            AuxBit | KeyBit | DensBit | EpsBit);
    instr = stropen(getparam("in"), "r");	
    outstr = stropen(getparam("out"), "w");

    opt = burststring(getparam("options"),", ");
    nopt = 0;					/* count options */
    while (opt[nopt]) {				/* scan through options */
        fopt[nopt] = btrtrans(opt[nopt]);
        nopt++;
        if (nopt==MAXOPT) {
            dprintf(0,"\n\nMaximum number of options = %d exhausted\n",MAXOPT);
            break;
        }
    }
    times = getparam("times");
    ntime = (hasvalue("ntime") ? getiparam("ntime") : 1);
    
    for(;;) {                /* repeating until first or all times are read */
	get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
            break;                                  /* done with work */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ( (bits & ParticlesBit) == 0)
            continue;                   /* skip work, only diagnostics here */
	dprintf(1,"Time=%g\n",tsnap);
	if (iptr == NULL) {
	  create_cube(&iptr,nbody,nopt,ntime);
	  k=0;
	  t0 = tsnap;
	  if (k==0) dt = 0.0;
	} 
	if (k==1) dt=tsnap-t0;    /* should be good for the whole snapshot */
	for (j=0; j<nopt; j++) {
	  for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
	    CubeValue(iptr,i,j,k) = fopt[j](bp,tsnap,i);
	  }
	}
	k++;

	if (k==ntime)  { /* cube is full */
	  fixheader(iptr,getparam("options"),t0,dt);
	  write_image(outstr,iptr);
	  free_image(iptr);
	  iptr = NULL;
	  k = 0;
	  if (Qfirst) break;
	}
    }
    if (!Qfirst && k) {
      warning("k=%d something not written yet, possible trailing garbage written",k);
      fixheader(iptr,getparam("options"),t0,dt);
      write_image(outstr,iptr);
    }
    strclose(instr);
    strclose(outstr);
}


void fixheader(imageptr iptr, string options, real t0, real dt)
{
  char *ap = allocate(sizeof(options) + sizeof("Attribute") + 10);
  sprintf(ap,"Attribute: %s",options);

  Namex(iptr) = strdup("Particle");
  Xmin(iptr)  = 0.0;
  Dx(iptr)    = 1.0;

  Namey(iptr) = strdup(ap);
  Ymin(iptr)  = 0.0;
  Dy(iptr)    = 1.0;

  Namez(iptr) = strdup("Time");
  Zmin(iptr)  = t0;
  Dz(iptr)    = dt;

  free(ap);
}
