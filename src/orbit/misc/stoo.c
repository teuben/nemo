/*
 *  STOO:  convert a binary snapshot file to a binary orbit
 *
 *	27-jul-87	V1.0 original version	Peter Teuben	
 *	28-jul-87	V2.0 new orbit(5)	PJT
 *	8-jun-88	V2.1 new filestruct	PJT
 *	30may-90	V2.2 fixed Masso()	PJT 
 *	21-nov-90	V2.3 new Nemo		PJT
 *	12-jun-91	V2.4 time->tsnap fixes bug, also process history PJT
 *	 7-mar-92	     make gcc happy
 *	25-may-92	V2.5 need <potential.h> now
 *	18-apr-95       V2.6 
 *      25-jul-2013     V3.0 allow multiple orbits output
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <potential.h>
#include <orbit.h>

string defv[] = {
  "in=???\n			snapshot input file name",
  "out=???\n			orbit output file name",
  "ibody=0\n			which particles to take (-1=all, 0=first)",
  "nsteps=10000\n		allocation default",
  "VERSION=3.0\n		25-jul-2013 PJT",
  NULL,
};

string usage="convert a snapshot file to an orbit";
string cvsid="$Id$";

string iname, oname;		/* file names */
stream instr,outstr;		/* file pointers */
string headline;

 	/* global snapshot stuff */

#define MOBJ 	8000
int    nobj, norb;
real   tsnap;			/* time of snapshot */
real   mass[MOBJ];
real   phase[MOBJ][2][NDIM];
int    key[MOBJ];
int    isel[MOBJ];            /* select for output orbit?*/
bool   Qtime, Qmass, Qkey;
	/* global orbit stuff */

int    ibody, nsteps;		/* allocation */
orbitptr optr;			/* pointer to orbit */

void setparams();
int read_snap();

void nemo_main()
{
  int i,j;
    
  setparams();
    
  instr  = stropen(iname,"r");		/* read from snapshot */
  outstr = stropen(oname,"w");		/* write to orbit */

  optr=NULL;
  allocate_orbit (&optr,NDIM,nsteps);
  
  for (j=0; j<norb; j++) {
    ibody = isel[j];
    i = 0;				/* counter of timesteps */
    if (j>0) rewind(instr);
    while (read_snap()) {		/* read until exausted */
      if (i==0) {			/* first time around */
	Masso(optr) = mass[ibody];
	I1(optr) = I2(optr) = I3(optr) = 0.0;
      }
      if (ibody>=nobj) {
	printf ("request to output ibody=%d, however nobj=%d\n",ibody,nobj);
	exit(1);
      }
      if (i>=nsteps) {
	printf ("Warning: too many timesteps requested, stopped at %d",i);
	break;
      }
      Torb(optr,i) = tsnap;
      Key(optr)    = key[ibody];
      Xorb(optr,i) = phase[ibody][0][0];	
      Yorb(optr,i) = phase[ibody][0][1];	
      Zorb(optr,i) = phase[ibody][0][2];	
      Uorb(optr,i) = phase[ibody][1][0];	
      Vorb(optr,i) = phase[ibody][1][1];	
      Worb(optr,i) = phase[ibody][1][2];
      i++;
    }
    Nsteps(optr) = i;			/* record actual number */
    if (j==0) put_history(outstr);
    write_orbit(outstr,optr);		/* write orbit to file */
  }
  strclose(outstr);			/* close files */
  strclose(instr);
}

void setparams()
{
  iname = getparam("in");
  oname = getparam("out");
  nsteps = getiparam("nsteps");
  norb = nemoinpi(getparam("ibody"),isel,MOBJ);
  if (norb < 0) error("ibody=%s bad",getparam("ibody"));
  
}

/*
 * READ_SNAP: read next snapshot from input stream
 */

int read_snap()
{				
  int i, cs;
    
  for(;;) {  /* loop until one snapshot found */
    get_history(instr);

    if (!get_tag_ok(instr,SnapShotTag)) {      /* we DO need a SnapShot */
	return 0;
    }
        
    get_set(instr, SnapShotTag);
      get_set(instr, ParametersTag);
        get_data(instr, NobjTag, IntType, &nobj, 0);
	if (nobj>MOBJ)
		error ("read_snap: not enough declared space to get data\n");
	if ((Qtime=get_tag_ok(instr,TimeTag)))
		get_data(instr,TimeTag,RealType,&tsnap,0);
	else
		tsnap=0.0;
      get_tes(instr,ParametersTag);
      
      if (!get_tag_ok(instr, ParticlesTag)) {		/* do it again, we need ParticlesTag */
						/*   normally happens because of DiagnosticsTag */
	 get_tes(instr,SnapShotTag);			/* close this SnapShotTag */
	 continue;					/* and loop again to read next one */
      }
      dprintf (1,"SnapShot: Time=%f ",tsnap);

      get_set(instr, ParticlesTag);
         get_data(instr, CoordSystemTag, IntType, &cs, 0);
         if (get_tag_ok(instr,MassTag)) {
           get_data(instr, MassTag, DoubleType, mass, nobj, 0);
           Qmass=TRUE;
         }
         else if (!Qmass) {	   
              dprintf (0,"Warning: no masses present: ASSUME 1/%d\n",nobj);
              for (i=0; i<nobj; i++)
                 mass[i] = 1.0/(double)nobj;
         }
         
         get_data(instr, PhaseSpaceTag, DoubleType, phase, nobj, 2, NDIM, 0);
	 /* need additional Qkey here ? */
	 if (get_tag_ok(instr,KeyTag))
	   get_data(instr, KeyTag, IntType, key, nobj, 0);
	 else {
	   dprintf(0,"added ordinal keys\n");
	   for (i=0; i<nobj; i++) key[i] = i;
	 }
      get_tes(instr,ParticlesTag);
    get_tes(instr,SnapShotTag);

    return 1;
  }
}

