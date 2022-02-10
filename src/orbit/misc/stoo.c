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
 *       6-feb-2022     V3.3 initialize I1, properly implement ibody=-1
 *       8-feb-2022     V3.4 stopped massive memory leak, add scan_snap() for files
 *                           use mdarray for dynamic snapshot
 *
 *  @todo    deal with no Potential/Accelleration
 */

#include <stdinc.h>
#include <unistd.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <mdarray.h>

#include <snapshot/snapshot.h>
#include <potential.h>
#include <orbit.h>

string defv[] = {
  "in=???\n			snapshot input file name",
  "out=???\n			orbit output file name",
  "ibody=-1\n			which particles to take (-1=all, 0=first)",
  "nsteps=0\n	                max orbit length allocated in case of pipe",
  "VERSION=3.4b\n		9-feb-2022 PJT",
  NULL,
};

string usage="convert a snapshot file to an orbit";
string cvsid="$Id$";

string iname, oname;		/* file names */
stream instr,outstr;		/* file pointers */
string headline;

int    nobj;                    /* number of bodies per snapshot */
int    norb;                    /* number of bodies to extract for orbits */
int    nsteps;                  /* number of snapshots to use for each orbit */
real   tsnap;			/* time of snapshot */

#define USE_MDARRAY

#if defined(USE_MDARRAY)
  int MOBJ = 0;
  real *mass;
  real ***phase;
  real *phi;
  real **acc;
  int  *key;
  int  *isel;
#else
  #define MOBJ 	1000000
  real   mass[MOBJ];
  real   phase[MOBJ][2][NDIM];
  real   phi[MOBJ];
  real   acc[MOBJ][NDIM];
  int    key[MOBJ];
  int    isel[MOBJ];              /* select for output orbit? */
#endif


bool   Qtime, Qmass, Qkey;

int    ibody;	
orbitptr *optr;			/* pointer to orbit pointers */

int read_snap();
int scan_snap();
void alloc_snap(int);

bool ispipe(stream instr)
{
  off_t try = lseek(fileno(instr), 0, SEEK_CUR);
  if (try < 0) return TRUE;
  return FALSE;
}

void nemo_main()
{
  int i,j;

  iname = getparam("in");
  oname = getparam("out");
  instr  = stropen(iname,"r");		/* read from snapshot */
  outstr = stropen(oname,"w");		/* write to orbit */
  nsteps = getiparam("nsteps");

  if (nsteps == 0)  {
    if (ispipe(instr))
      error("Need to set nsteps= for pipe files");

    nsteps = scan_snap();
    dprintf(0,"Found %d snapshots with nobj=%d\n",nsteps,nobj);
    alloc_snap(nobj);
    rewind(instr);
  } else {
    dprintf(0,"Setting orbits for %d snapshots\n",nsteps);
  }

  i = 0;				/* counter of snapshots/timesteps */
  while (read_snap()) {
      if (i==0) {			/* first time around */
	norb = nemoinpi(getparam("ibody"),isel,MOBJ);
	if (norb < 0)  error("%d: ibody=%s bad",norb,getparam("ibody"));
	if (norb == 0) error("no orbits will be output");
	
	if (isel[0] < 0) {              /* special case selecting all */
	  dprintf(1,"Selecting all %d bodies for conversion to an orbit\n", nobj);
	  norb = nobj;
	  for (j=0; j<norb; j++)
	    isel[j] = j;
	}
	
	optr=(orbitptr *) allocate(norb*sizeof(orbitptr));
	for (j=0; j<norb; j++) {
	  optr[j] = NULL;
	  allocate_orbit (&optr[j],NDIM,nsteps);
	  
	  ibody = isel[j];
	  if (ibody < 0 || ibody>=nobj)
	    error("illegal ibody=%d,  nobj=%d",ibody,nobj);
	  Masso(optr[j]) = mass[ibody];
	  I1(optr[j]) = I2(optr[j]) = I3(optr[j]) = 0.0;
	}
	dprintf(0,"Selecting %d bodies for orbits, %d..%d\n",norb,isel[0],isel[norb-1]);
  
      } else if (i>=nsteps) {
	warning("more snapshots found, stopped at %d",i);
	break;
      }
      for (j=0; j<norb; j++) {
	ibody = isel[j];
	Key(optr[j])    = key[ibody];
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
	if (i==0) 
	  I1(optr[j]) = phi[ibody] + 0.5*(sqr(phase[ibody][1][0]) + sqr(phase[ibody][1][1]) + sqr(phase[ibody][1][2]));
#endif
      }
      i++;
      for (j=0; j<norb; j++)
	Nsteps(optr[j]) = i;			/* record actual number */
      progress(1.0,"processed snapshot %d/%d", i,nsteps);      
  }
  strclose(instr);
  dprintf(0,"Writing %d orbits\n",norb);

  if (norb > 0) {
    put_history(outstr);
    for (j=0; j<norb; j++) {
      progress(1.0,"writing orbit %d", j);
      write_orbit(outstr,optr[j]);		/* write orbit to file */
      free_orbit(optr[j]);
    }
    free(optr);
  }
  dprintf(0,"Processed %d snapshots\n",i);
  
  strclose(outstr);		        	/* close files */
}



void alloc_snap(int nbody)
{
#if defined(USE_MDARRAY)
  dprintf(1,"alloc_snap(%d)\n",nbody);
  mass = (real *) allocate_mdarray1(nbody);
  phase = (real ***) allocate_mdarray3(nbody,2,NDIM);
  phi =  (real *) allocate_mdarray1(nbody);
  acc = (real **) allocate_mdarray2(nbody,NDIM);
  key =  (int *) allocate(nbody * sizeof(int));
  isel = (int *) allocate(nbody * sizeof(int));
  MOBJ = nbody;
#else
  warning("alloc_snap: already allocated statically");
#endif  
}
  

/*
 * SCAN_SNAP:   scan to find how many snapshots there are;
 *              also ensuring they are all the same size
 */

int scan_snap()
{
  int nsnap = 0;
  int maxobj = 0;

  dprintf(0,"Scanning snapshot\n");
  for(;;) {
    get_history(instr);
    if (!get_tag_ok(instr,SnapShotTag))
      return nsnap;

     get_set(instr,SnapShotTag);
     get_set(instr, ParametersTag);
     get_data(instr, NobjTag, IntType, &nobj, 0);
     if (maxobj == 0) 
       maxobj = nobj;
     if (maxobj != nobj)
       error("Cannot process snapshots with unequal number of bodies: %d != %d", nobj,maxobj);
     get_tes(instr, ParametersTag);
     if (!get_tag_ok(instr,ParticlesTag)) {
       get_tes(instr,SnapShotTag);
       continue;
     }
     get_tes(instr,SnapShotTag);
     nsnap++;
  }
}

/*
 * READ_SNAP: read next snapshot from input stream
 *
 */

int read_snap()
{				
  int i, cs;
  static bool first = TRUE;
  static int no_mass = 0;
    
  for(;;) {  /* loop until one snapshot found */
    get_history(instr);

    if (!get_tag_ok(instr,SnapShotTag)) {      /* we DO need a SnapShot */
        if (no_mass)
	    warning("There were %d snapshots with no masses. Total mass=1 was set" , no_mass);
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
	if (MOBJ == 0) {
	  alloc_snap(nobj);
	  // this will set MOBJ if allowed
	}
	if (nobj>MOBJ)
	  error ("read_snap: not enough memory; only space for %d bodies",MOBJ);
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
           get_data_coerced(instr, MassTag, RealType, mass, nobj, 0);
           Qmass=TRUE;
         }
         else if (!Qmass) {
     	      no_mass++;
              dprintf (1,"Warning: no masses present: ASSUME 1/%d\n",nobj);
              for (i=0; i<nobj; i++)
                 mass[i] = 1.0/(double)nobj;
         }
         if (get_tag_ok(instr,PhaseSpaceTag))
	   get_data(instr, PhaseSpaceTag, RealType, &phase[0][0][0], nobj, 2, NDIM, 0);
	 else {
	   get_data_coerced(instr, PosTag, RealType, &phase[0][0][0], nobj, NDIM, 0);
	   get_data_coerced(instr, VelTag, RealType, &phase[0][1][0], nobj, NDIM, 0);
	 }
	 /* need additional Qkey here ? */
	 if (get_tag_ok(instr,KeyTag))
	   get_data(instr, KeyTag, IntType, key, nobj, 0);
	 else {
	   dprintf(1,"added ordinal keys\n");
	   for (i=0; i<nobj; i++) key[i] = i;
	 }
#ifdef ORBIT_PHI
	 get_data_coerced(instr,PotentialTag,RealType,phi,nobj,0);
	 get_data_coerced(instr,AccTag,RealType,&acc[0][0],nobj,NDIM,0);
#endif	 
      get_tes(instr,ParticlesTag);
    get_tes(instr,SnapShotTag);

    return 1;
  }
}

