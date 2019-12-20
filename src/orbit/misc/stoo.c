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
 *      25-jul-2013     V3.0 allow multiple orbits output;
 *      11-nov-2015     V3.1 now the efficient way
 *      19-dec-2019     V3.2 deal with ACC/PHI
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
  "nsteps=10000\n		max orbit length allocated",
  "VERSION=3.2\n		19-dec-2019 PJT",
  NULL,
};

string usage="convert a snapshot file to an orbit";
string cvsid="$Id$";

string iname, oname;		/* file names */
stream instr,outstr;		/* file pointers */
string headline;

 	/* global snapshot stuff */

#define MOBJ 	1000000
int    nobj;                    /* number of bodies per snapshot */
int    norb;                    /* number of orbits to extract */
real   tsnap;			/* time of snapshot */
real   mass[MOBJ];
real   phase[MOBJ][2][NDIM];
real   phi[MOBJ];
real   acc[MOBJ][NDIM];
int    key[MOBJ];
int    isel[MOBJ];              /* select for output orbit? */
bool   Qtime, Qmass, Qkey;
	/* global orbit stuff */

int    ibody, nsteps;		/* allocation */
orbitptr *optr;			/* pointer to orbit pointers */

void setparams();
int read_snap();

void nemo_main()
{
  int i,j;
    
  setparams();
    
  instr  = stropen(iname,"r");		/* read from snapshot */
  outstr = stropen(oname,"w");		/* write to orbit */

  optr=(orbitptr *) allocate(norb*sizeof(orbitptr));
  for (i=0; i<norb; i++) {
    optr[i] = NULL;
    allocate_orbit (&optr[i],NDIM,nsteps);
  }
  
      
  i = 0;				/* counter of timesteps */
  while (read_snap()) {		         /* read until exausted */
      if (i==0) {			/* first time around */
	for (j=0; j<norb; j++) {
	  ibody = isel[j];
	  if (ibody>=nobj)
	    error("request to output ibody=%d, however nobj=%d",ibody,nobj);
	  Masso(optr[j]) = mass[ibody];
	  I1(optr[j]) = I2(optr[j]) = I3(optr[j]) = 0.0;
	}
      } else {
	if (i>=nsteps) {
	  warning("too many timesteps requested, stopped at %d",i);
	  break;
	  }
      }
      for (j=0; j<norb; j++) {
	ibody = isel[j];
	Key(optr[0])    = key[ibody];
	Torb(optr[j],i) = tsnap;
	Xorb(optr[j],i) = phase[ibody][0][0];	
	Yorb(optr[j],i) = phase[ibody][0][1];	
	Zorb(optr[j],i) = phase[ibody][0][2];	
	Uorb(optr[j],i) = phase[ibody][1][0];	
	Vorb(optr[j],i) = phase[ibody][1][1];	
	Worb(optr[j],i) = phase[ibody][1][2];
#ifdef ORBIT_PHI
	Porb(optr[j],i) = phi[ibody];
	AXorb(optr[j],i) = acc[ibody][0];
	AYorb(optr[j],i) = acc[ibody][1];
	AZorb(optr[j],i) = acc[ibody][2];
#endif
      }
      i++;
      for (j=0; j<norb; j++)
	Nsteps(optr[j]) = i;			/* record actual number */      
  }
  put_history(outstr);
  for (j=0; j<norb; j++)
    write_orbit(outstr,optr[j]);		/* write orbit to file */

  
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
  if (norb == 0) warning("no orbits will be output");
  
}

/*
 * READ_SNAP: read next snapshot from input stream
 */

int read_snap()
{				
  int i, cs;
  static bool first = TRUE;
    
  for(;;) {  /* loop until one snapshot found */
    get_history(instr);

    if (!get_tag_ok(instr,SnapShotTag)) {      /* we DO need a SnapShot */
	return 0;
    }
        
    get_set(instr, SnapShotTag);
      get_set(instr, ParametersTag);
        get_data(instr, NobjTag, IntType, &nobj, 0);
	if (first) {
	  dprintf(0,"read_snap: Found first snapshot with %d bodies\n",nobj);
	  first = FALSE;
	} 
	dprintf(1,"."); fflush(stderr);
	if (nobj>MOBJ)
	  error ("read_snap: not enough declared space to get data");
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
           get_data(instr, MassTag, RealType, mass, nobj, 0);
           Qmass=TRUE;
         }
         else if (!Qmass) {	   
              dprintf (0,"Warning: no masses present: ASSUME 1/%d\n",nobj);
              for (i=0; i<nobj; i++)
                 mass[i] = 1.0/(double)nobj;
         }
         
         get_data(instr, PhaseSpaceTag, RealType, phase, nobj, 2, NDIM, 0);
	 /* need additional Qkey here ? */
	 if (get_tag_ok(instr,KeyTag))
	   get_data(instr, KeyTag, IntType, key, nobj, 0);
	 else {
	   dprintf(1,"added ordinal keys\n");
	   for (i=0; i<nobj; i++) key[i] = i;
	 }
#ifdef ORBIT_PHI
	 get_data(instr,PotentialTag,RealType,phi,nobj,0);
	 get_data(instr,AccTag,RealType,acc,nobj,NDIM,0);
#endif	 
      get_tes(instr,ParticlesTag);
    get_tes(instr,SnapShotTag);

    return 1;
  }
}

