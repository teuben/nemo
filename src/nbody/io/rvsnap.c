/*
 *  RVSNAP: convert Carlberg's binary 'RV' format to snapshot format
 *	    Also known to be used by Dubinski
 *
 *	14-may-97	V1.0	Created				PJT
 *	
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <stdlib.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/put_snap.c>


string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (RV format)",
    "out=???\n			Output file (snapshot)",
    "VERSION=1.0\n		14-may-97",
    NULL,
};

string usage = "convert Carlberg's binary 'RV' format to snapshot format";

local int  rv_header(stream);
local void rv_data(stream, int, Body **, real *);

void nemo_main(void)
{
    stream instr, outstr;
    real   tsnap;
    string times;
    Body *btab = NULL, *bp;
    int j, nbody, bits;

            
    instr = stropen(getparam("in"), "r");	/* open input file (no pipe!) */
    outstr = stropen(getparam("out"), "w");	/* open output file */

    for(;;) {                /* repeating until first or all times are read */
        nbody = rv_header(instr);
        if (nbody <= 0) break;

	rv_data(instr, nbody, &btab, &tsnap);

	put_history(outstr);
        bits = (TimeBit | MassBit | PhaseSpaceBit);
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);

        break;
    }        

    strclose(instr);
    strclose(outstr);
}


/*
 * INITDATA: read initial conditions from source_file
 * Input is binary unformatted and assumed to be in the following form:
 * nobj,(nobj)masses,tnow,(nobj)(x,y,z,vx,vy,vz)
 * This file will be appended to and end up in the form:
 * nobj,nobj*masses,frames*(tnow,nobj*(x,y,z,vx,vy,vz))
 * Hence fixed data is in bytes [0,4*(1+nobj)-1]
 *   frame(0) starts at byte 4*(1+nobj)
 *	 frame(1) starts at byte 4*{(1+nobj) + (1+6*nobj)}
 *	 frame(2) starts at byte 4*{(1+nobj) + 2*(1+6*nobj)}
 *	 ...
 *	 frame(k) starts at byte 4*{(1+nobj) + k*(1+6*nobj)}
 * Hence a file containing N complete frames will have length
 *	4*{ 1+nobj + N*(1+6*nobj) } bytes.
 * Hence the number of frames in a file of size B bytes will be
 *	N = (B/4 - (1+nobj)) / (1+6*nobj)
 * This file thus should contain frames 0,1,2,...,N-1.
 * This thus allows for restarting from the last used frames.
 * NOTE: All real numbers are saved as type 'float' even if the simulation
 *       is run with double precision.  Also the source_file is assumed to
 *       consist of 'float' type numbers even if the simulation in run in
 *       double precision.
 */

local int rv_header(stream instr)
{
    int frame, nobj, i;
    float tnow_f, mass_f, pos_f[3], vel_f[3];

    fseek(instr,0,2);	        		/* go to eof of source_file */
    frame = ftell(instr);	        	/* get length of source_file */
    fseek(instr,0,0);		        	/* back to beginning */
    fread(&nobj,sizeof(nobj),1,instr);	        /* read number of bodies */
    dprintf(0,"Number of bodies found = %d\n",nobj);

    return nobj;
}


local void rv_data(stream instr, int nobj, Body **bodytab, real *tsnap)
{
    Body *p, *bp, *btab, *firstbody, *lastbody;
    real mtot;
    float tnow_f, mass_f, pos_f[3], vel_f[3];
    int i, frame;
    
    btab = (body *) allocate(nobj*sizeof(body));  /* allocate body table */
    
    firstbody = btab; 			/* Address of first body */	
    lastbody  = btab + nobj - 1;	/* Address of last body */

    mtot = 0.0;   /* Zero the total mass */
    for (p = btab; p < btab + nobj; p++) {	/* loop to get masses */
	fread(&mass_f,sizeof(float),1,instr);		/* read mass */
	Mass(p) = (real) mass_f;               /* Convert mass_f to real */
	mtot += Mass(p);		       /* Determine the system's mass */
    }

    frame = (frame/4 - nobj - 1)/(1 + 6*nobj); /* compute # of frames */

    fseek(instr,4*( 1+ nobj + (frame-1)*(1+6*nobj) ),0); 
							/* go to beginning of the last frame */
    fread(&tnow_f,sizeof(float),1,instr);		/* read current time */
    *tsnap= (real) tnow_f;        			/* Convert to real */
    dprintf(0,"Time: %f\n",*tsnap);
    for (p = btab; p < btab + nobj; p++) {
	fread(pos_f, sizeof(pos_f), 1, instr);	        /* x, y, z */
	fread(vel_f, sizeof(vel_f), 1, instr);  	/* vx, vy, vz */
	for(i=0; i<NDIM; i++) {
	    Pos(p)[i] = (real) pos_f[i];    /* Convert to real */
	    Vel(p)[i] = (real) vel_f[i];
	}
    }

    *bodytab = btab;
}
