/*
 *  SNAPIDL: tabulate a snapshot in binary mode for IDL
 *
 *      26-nov-2002    1.0 Test example for Walter Dehnen             Peter Teuben
 *                     1.1 add fortran= and output nvec also          PJT
 *       5-aug-2009    1.2 header 4 or 8 issue
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (snapshot)",
    "options=x,y,z,vx,vy,vz\n	Things to output",
    "times=all\n		Times to select snapshot",
    "out=-\n                    Output file - normally a pipe for IDL",
    "fortran=f\n                Add fortran unformatted I/O tags",
    "header=\n                  Header size (4, 8, or auto)",
    "VERSION=1.1\n		26-nov-02 PJT",
    NULL,
};

string usage="tabulate a snapshot for IDL in binary mode";

string cvsid="$Id$";

#define MAXOPT    50

void fortout(stream os, int nz);

int header= -1;

void nemo_main()
{
    stream instr, outstr;
    real   tsnap, dr, aux;
    string times;
    Body *btab = NULL, *bp, *bq;
    bool Qfort = getbparam("fortran");
    float  *vec, sca;
    int i, n, nbody, bits, nsep, isep, nopt, ParticlesBit, nfort;
    char fmt[20],*pfmt;
    string *burststring(), *opt;
    rproc btrtrans(), fopt[MAXOPT];

    ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
            AuxBit | KeyBit | DensBit | EpsBit);
    instr = stropen(getparam("in"), "r");	
    outstr = stropen(getparam("out"), "w");
    header = getiparam("header");
    setbuf(outstr,(char *)0);                   /* suggested by IDL, inc */

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

    for(;;) {                /* repeating until first or all times are read */
	get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
            break;                                  /* done with work */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ( (bits & ParticlesBit) == 0)
            continue;                   /* skip work, only diagnostics here */
	vec = (float *) allocate(nbody*sizeof(float));
	sca = tsnap;
	
	if (Qfort) fortout(outstr,2*sizeof(int)+sizeof(float));
	fwrite(&nbody, sizeof(int),   1, outstr);
	fwrite(&nopt,  sizeof(int),   1, outstr);
	fwrite(&sca,   sizeof(float), 1, outstr);
	if (Qfort) fortout(outstr,2*sizeof(int)+sizeof(float));
	for (n=0; n<nopt; n++) {
	  for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
	    vec[i] = (float) fopt[n](bp,tsnap,i);	  }
	  if (Qfort) fortout(outstr,nbody*sizeof(float));	  
	  fwrite(vec, sizeof(float), nbody, outstr);
	  if (Qfort) fortout(outstr,nbody*sizeof(float));	  
	}
	free(vec);
    }
    strclose(instr);
    strclose(outstr);
}

/*
 *  this is quite possibly quite unportable, but works on at least solaris
 *  and linux.
 *
 *  TODO: should really use unfio() routines for this !!!
 */

void fortout(stream os, int nz)
{
  int ivar4       = nz;
  long long ivar8 = nz;

  if (header < 0)
    fwrite(&nz,sizeof(int),1,os);
  else if (header == 0)
    dprintf(1,"Skip writing fortran block counter\n");
  else if (header == sizeof(ivar4))
    fwrite(&ivar4,sizeof(ivar4),1,os);
  else if (header == sizeof(ivar8))
    fwrite(&ivar8,sizeof(ivar8),1,os);
  else
    error("header=%d does not match a sizeof(type) for fortran I/O: %d or %d",
	  header,sizeof(ivar4),sizeof(ivar8));

}
