/*
 *  OTOS:  convert a binary orbit to a snapshot file
 *
 *	23-may-90	V1.0	created in great haste	PJT
 *	20-dec-90	V1.1	scan added - more verbose   PJT
 *	 7-mar-92	   a    gcc happy		   PJT
 *	25-may-92	V1.2    need <potential.h> now 	   PJT
 *	25-jun-92          a    allocate decl. error       PJT
 *	19-apr-96	   b    fixed stdout/sterr 	   pjt
 *       1-apr-01          c    compiler warning	   pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <potential.h>
#include <orbit.h>

#include <snapshot/body.h>
#include <snapshot/snapshot.h>
#include <snapshot/put_snap.c>


string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			   Orbit input file name",
    "out=???\n			   Snapshot output file name",
    "orbit=all\n		   Which orbits to take",
    "maxsteps=5000\n		   Orbit size allocation",
    "stride=1\n                    Orbit stride to snapshot",
    "nbody=10000\n                 Snapshot size allocation",
    "scan=f\n                      Scan to end of all orbits for counting?",
    "VERSION=1.2c\n		   1-apr-01 PJT",
    NULL,
};

string usage = "convert a binary orbit to a snapshot file";

#define MAXORBITS  512

stream instr,outstr;		/* file pointers */

	/* global snapshot stuff */
int    nbody, maxbody;
real   tsnap;			/* time of snapshot */
Body *btab=NULL;

	/* global orbit stuff */
int    maxsteps, stride;	/* allocation */
orbitptr optr;			/* pointer to orbit */

int oconvert[MAXORBITS];	/* ordinal of all orbits to convert (0,1,...) */

void nemo_main()
{
    int i, bits, norbit, nphase, nread, idx;
    Body *bp;
    bool oscan;
    string olist;

    maxsteps = getiparam("maxsteps");
    maxbody = getiparam("nbody");
    stride = getiparam("stride");    
    oscan = getbparam("scan");
    olist = getparam("orbit");
    if (*olist != 0) {
        if (streq(olist,"all"))
            nread = 0;
        else {
            nread = nemoinpi(olist,oconvert,MAXORBITS);
            if (nread<0) error("Parsing error orbit=%d",olist);
        }
    } else
        nread = 0;
    instr  = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");

    optr=NULL;
    allocate_orbit (&optr,NDIM,maxsteps);
    btab = (Body *) allocate(sizeof(Body)*maxbody);

    bp = btab;
    nbody=0;    /* count snapshot */
    nphase=0;   /* count how many we've from orbit */
    norbit=0;   /* count orbits in input */
    idx = 0;    /* index in oconvert[] array */
    while (read_orbit(instr,&optr)) {           /* go over all orbits */
        if (nread>0) {                      /* if to select orbits */
            if (oconvert[idx] == norbit) {  /* doesn't warn for ill order */
                idx++;          /* increase pointer */
                norbit++;       /* and orbit counter - and copy it */
            } else {
                norbit++;       /* increase orbit counter also */
                continue;       /* but do not copy - go to next orbit */
            }
        } else
            norbit++;           /* increase orbit counter */    
        for (i=0; i<Nsteps(optr); i+=stride) {
            nphase++;
            if (nbody < maxbody) {
                bp = &btab[nbody];
                Pos(bp)[0] = Xorb(optr,i);
	        Pos(bp)[1] = Yorb(optr,i);
	        Pos(bp)[2] = Zorb(optr,i);
	        Vel(bp)[0] = Uorb(optr,i);
	        Vel(bp)[1] = Vorb(optr,i);
	        Vel(bp)[2] = Worb(optr,i);
                bp++;       /* point to next particle */
                nbody++;
            } else {
                if (!oscan) {
                  warning("All %d particles in snapshot filled. Use scan=t",
                                maxbody);
                  nphase--;
                  break;
                }
            }
        }
        if (nbody >= maxbody && !oscan)
            break;  /* once more - just to avoid a goto */
	dprintf(1,"Read orbit, stored %d so far\n",nbody);
        if (nread>0 && idx>=nread && !oscan)     /* early bail out ? */
                break;                      /* if no scan, early bail out */ 
	Nsteps(optr) = maxsteps;
    }
    if (nbody<nphase)
        warning("Snapshot incomplete");
    dprintf(0,"Read total of %d phase points from %d orbits",nphase,norbit);
    if (nbody<nphase)
        dprintf(0,";\n %d are stored in snapshot.\n",nbody);
    else
        dprintf(0," into snapshot.\n");
    
    put_history(outstr);
    bits = (TimeBit | PhaseSpaceBit);
    tsnap = 0;
    put_snap(outstr,&btab,&nbody,&tsnap,&bits);

    strclose(instr);
    strclose(outstr);			/* close files */
}
