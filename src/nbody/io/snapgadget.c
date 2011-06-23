/*
 *  SNAPGADGET: convert a NEMO snapshot to GADGET format
 *
 * Jeremy Bailin September 2003, based on snaptipsy.c (by TRQ + PJT)
 * and JB's file_gadget.c.
 *
 * V0.2   fixed bug when no %t in output filename, re-arranged order of keys    PJT
 *
 * TODO:  should use mstropen() instead of the dangerous manual parsing of %t->%f
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <stdlib.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

struct io_header_1
{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;

#if 0
  /* newer versions may have this - need to recount the filler space now */
   int flag_multiphase;
   int flag_stellarage;
   int flag_sfrhistogram;
   int flag_metals;
   int flag_decouple;
   int flag_effmodel;
#endif

    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8]; /* fills to 256 bytes */
};


string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (NEMO snapshot)",
    "out=???\n                  Output file (GADGET format), %t for time",
    "N=???\n			Nhalo,Ndisk,Nbulge,Nstars (in that order)",
    "mass=\n                    Mh, Md, Mb, Ms (if no masses in input)",
    "times=all\n		Times to select snapshot",
    "swap=f\n                   Swap bytes on output?",
    "z=0\n                      ** Current Redshift",
    "box=128\n                  ** Box size (Mpc/h)",
    "omega0=0.3\n               ** Omega0",
    "lambda0=0.7\n              ** OmegaLambda0",
    "h=0.75\n                   ** HubbleParam",
    "VERSION=1.1\n		22-jun-2011 PJT",
    NULL,
};

string usage="convert snapshot to GADGET format";

string cvsid="$Id$";

#define TIMEFORMAT "%t"

extern void bswap(void *, int, int);
void write_gadget(stream outstr,real time,Body *bodies,int nhalo,int ndisk,
		  int nbulge,int nstars,bool swapp);


void nemo_main()
{
    stream instr, outstr;
    real   tsnap;
    string times, mode;
    char outfname[256],outstub[256];
    char *fnpt;
    bool Qswap;
    Body *btab = NULL, *bp;
    int i, ndim, nbody, bits, ParticlesBit, omode;
    int  nh,nd,nb,ns;
    real mh,md,mb,ms;
    string *burststring();
    rproc btrtrans();
    struct dark_particle *dark;

    ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
            AuxBit | KeyBit);
    instr = stropen(getparam("in"), "r");	/* open input file */

    times = getparam("times");
    Qswap = getbparam("swap");
    get_history(instr);                 /* read history */

    for(;;) {                /* repeating until first or all times are read */
	get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
            break;                                  /* done with work */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ( (bits & ParticlesBit) == 0)
            continue;                   /* skip work, only diagnostics here */

	/* replace all instances of '%t' with the snapshot time */
	strcpy(outstub, getparam("out"));
	strcpy(outfname,outstub);
	for(fnpt=strstr(outstub,TIMEFORMAT); fnpt!=NULL; ) {
	  *(fnpt+1)='f';
	  sprintf(outfname, outstub, tsnap);
	  strcpy(outstub, outfname);
	  fnpt=strstr(outstub,TIMEFORMAT);
	}
	dprintf(1,"Writing t=%f to \"%s\"\n",tsnap,outfname);
        outstr = stropen(outfname,"w");

	/* @todo:   should really use nemoinpi/nemoinpd here */
	sscanf(getparam("N"), "%d,%d,%d,%d", &nh, &nd, &nb, &ns);
	if(nbody != nh + nd  + nb + ns) {
	  error("nbody != nh + nd + nb + ns");
	}
	if (hasvalue("mass")) {
	  sscanf(getparam("mass"),"%lg,%lg,%lg,%lg", &mh,&md,&mb,&ms);
	  if (bits&MassBit) warning("Overwriting masses");
	  for (i=0; i<nh; i++)  Mass(btab+i) = mh;
    	  for (i=nh; i<nh+nd; i++) Mass(btab+i) = md;
    	  for (i=nh+nd; i<nh+nd+nb; i++) Mass(btab+i) = mb;
    	  for (i=nh+nd+nb; i<nh+nd+nb+ns; i++) Mass(btab+i) = ms;
	}
	write_gadget(outstr, tsnap, btab, nh, nd, nb, ns, Qswap);
	strclose(outstr);
    }
    strclose(instr);
}


/* write FORTRAN block length fields */
#define BLKLEN if(swapp) bswap(&blklen,sizeof(blklen),1); fwrite(&blklen,\
        sizeof(blklen), 1, outstr); if(swapp) \
        bswap(&blklen,sizeof(blklen),1)

void write_gadget(stream outstr,real time,Body *bodies,int nhalo,int ndisk,
		  int nbulge,int nstars,bool swapp)
{
  struct io_header_1 header;
  bool indivmass[4];
  double pmass;
  int i,k,blklen,np,ntotwithmass;
  float xyz[3];
  int ibuf;

  /* check if all the masses for a given particle type are the same,
     in which case we should use the mass field in the header rather
     than giving each particle an individual mass */

  /* 0: gas */
  header.mass[0] = 0.0;

  /* 1: halo */
  indivmass[0]=FALSE;
  if(nhalo>=1) {
    pmass = Mass(bodies+0);
    for(i=1 ; i<nhalo ; i++) {
      if( Mass(bodies+i) != pmass ) {
	/* individual masses */
	indivmass[0]=TRUE;
      }
    }
  } else {
    pmass=0.0;
  }
  if(indivmass[0]) {
    header.mass[1] = 0.0;
  } else {
    header.mass[1] = pmass;
  }

  /* 2: disk */
  indivmass[1]=FALSE;
  if(ndisk>=1) {
    pmass = Mass(bodies+nhalo);
    for(i=1 ; i<ndisk ; i++) {
      if( Mass(bodies + nhalo+ i) != pmass ) {
	/* individual masses */
	indivmass[1]=TRUE;
      }
    }
  } else {
    pmass=0.0;
  }
  if(indivmass[1]) {
    header.mass[2] = 0.0;
  } else {
    header.mass[2] = pmass;
  }

  /* 3: bulge */
  indivmass[2]=FALSE;
  if(nbulge>=1) {
    pmass = Mass(bodies + nhalo+ndisk);
    for(i=1 ; i<nbulge ; i++) {
      if( Mass(bodies + nhalo+ndisk+ i) != pmass ) {
	/* individual masses */
	indivmass[2]=TRUE;
      }
    }
  } else {
    pmass=0.0;
  }
  if(indivmass[2]) {
    header.mass[3] = 0.0;
  } else {
    header.mass[3] = pmass;
  }

  /* 4: stars */
  indivmass[3]=FALSE;
  if(nstars>=1) {
    pmass = Mass(bodies + nhalo+ndisk+nbulge);
    for(i=1 ; i<nstars ; i++) {
      if( Mass(bodies + nhalo+ndisk+nbulge+ i) != pmass ) {
	/* individual masses */
	indivmass[3]=TRUE;
      }
    }
  } else {
    pmass=0.0;
  }
  if(indivmass[3]) {
    header.mass[4] = 0.0;
  } else {
    header.mass[4] = pmass;
  }

  /* 5: boundary */
  header.mass[5] = 0.0;


  /* rest of the header */
  header.num_files    = 0;
  header.time         = time;              /* inherit from snapshot */
  header.redshift     = getrparam("z");
  header.flag_sfr     = 0;                 /* fake */
  header.flag_feedback= 0;                 /* fake */
  header.flag_cooling = 0;                 /* fake */
  header.BoxSize      = getdparam("box");
  header.Omega0       = getdparam("omega0");
  header.OmegaLambda  = getdparam("lambda0");
  header.HubbleParam  = getdparam("h");
  header.npart[0] = header.npartTotal[0] = 0;      /* no gas */
  header.npart[1] = header.npartTotal[1] = nhalo;
  header.npart[2] = header.npartTotal[2] = ndisk;
  header.npart[3] = header.npartTotal[3] = nbulge;
  header.npart[4] = header.npartTotal[4] = nstars;
  header.npart[5] = header.npartTotal[5] = 0;

  /* swap if necessary */
  if(swapp) {
    bswap(&header.time, sizeof(double), 1);
    bswap(&header.redshift, sizeof(double), 1);
    bswap(&header.flag_sfr, sizeof(int), 1);
    bswap(&header.flag_feedback, sizeof(int), 1);
    bswap(&header.flag_cooling, sizeof(int), 1);
    bswap(&header.num_files, sizeof(int), 1);
    bswap(&header.BoxSize, sizeof(double), 1);
    bswap(&header.Omega0, sizeof(double), 1);
    bswap(&header.OmegaLambda, sizeof(double), 1);
    bswap(&header.HubbleParam, sizeof(double), 1);
    bswap(&header.npart, sizeof(int), 6);
    bswap(&header.mass, sizeof(double), 6);
    bswap(&header.npartTotal, sizeof(int), 6);
  }

  /* write header */
  blklen=sizeof(header);
  BLKLEN;
  fwrite(&header, sizeof(header), 1, outstr);
  BLKLEN;

  /* positions */
  np = nhalo + ndisk + nbulge + nstars;
  blklen = 3*np*sizeof(float);
  BLKLEN;
  for(i=0 ; i<nhalo ; i++) {
    for(k=0 ; k<3 ; k++) {
      xyz[k] = Pos(bodies + i)[k];
    }
    if(swapp) {
      bswap(&xyz, sizeof(float), 3);
    }
    fwrite(xyz,sizeof(float),3,outstr);
  }
  for(i=0 ; i<ndisk ; i++) {
    for(k=0 ; k<3 ; k++) {
      xyz[k] = Pos(bodies + nhalo+ i)[k];
    }
    if(swapp) {
      bswap(&xyz, sizeof(float), 3);
    }
    fwrite(xyz,sizeof(float),3,outstr);
  }
  for(i=0 ; i<nbulge ; i++) {
    for(k=0 ; k<3 ; k++) {
      xyz[k] = Pos(bodies + nhalo+ndisk+ i)[k];
    }
    if(swapp) {
      bswap(&xyz, sizeof(float), 3);
    }
    fwrite(xyz,sizeof(float),3,outstr);
  }
  for(i=0 ; i<nstars ; i++) {
    for(k=0 ; k<3 ; k++) {
      xyz[k] = Pos(bodies + nhalo+ndisk+nbulge+ i)[k];
    }
    if(swapp) {
      bswap(&xyz, sizeof(float), 3);
    }
    fwrite(xyz,sizeof(float),3,outstr);
  }
  BLKLEN;

  /* velocities */
  BLKLEN;
  for(i=0 ; i<nhalo ; i++) {
    for(k=0 ; k<3 ; k++) {
      xyz[k] = Vel(bodies + i)[k];
    }
    if(swapp) {
      bswap(&xyz, sizeof(float), 3);
    }
    fwrite(xyz,sizeof(float),3,outstr);
  }
  for(i=0 ; i<ndisk ; i++) {
    for(k=0 ; k<3 ; k++) {
      xyz[k] = Vel(bodies + nhalo+ i)[k];
    }
    if(swapp) {
      bswap(&xyz, sizeof(float), 3);
    }
    fwrite(xyz,sizeof(float),3,outstr);
  }
  for(i=0 ; i<nbulge ; i++) {
    for(k=0 ; k<3 ; k++) {
      xyz[k] = Vel(bodies + nhalo+ndisk+ i)[k];
    }
    if(swapp) {
      bswap(&xyz, sizeof(float), 3);
    }
    fwrite(xyz,sizeof(float),3,outstr);
  }
  for(i=0 ; i<nstars ; i++) {
    for(k=0 ; k<3 ; k++) {
      xyz[k] = Vel(bodies + nhalo+ndisk+nbulge+ i)[k];
    }
    if(swapp) {
      bswap(&xyz, sizeof(float), 3);
    }
    fwrite(xyz,sizeof(float),3,outstr);
  }
  BLKLEN;

  /* id */
  blklen = np*sizeof(int);
  BLKLEN;
  for(i=0 ; i<nhalo ; i++) {
    ibuf = i;
    if(swapp) {
      bswap(&ibuf, sizeof(int), 1);
    }
    fwrite(&ibuf,sizeof(int),1,outstr);
  }
  for(i=0 ; i<ndisk ; i++) {
    ibuf = nhalo+ i;
    if(swapp) {
      bswap(&ibuf, sizeof(int), 1);
    }
    fwrite(&ibuf,sizeof(int),1,outstr);
  }
  for(i=0 ; i<nbulge ; i++) {
    ibuf = nhalo+ndisk+ i;
    if(swapp) {
      bswap(&ibuf, sizeof(int), 1);
    }
    fwrite(&ibuf,sizeof(int),1,outstr);
  }
  for(i=0 ; i<nstars ; i++) {
    ibuf = nhalo+ndisk+nbulge+ i;
    if(swapp) {
      bswap(&ibuf, sizeof(int), 1);
    }
    fwrite(&ibuf,sizeof(int),1,outstr);
  }
  BLKLEN;

  /* masses, where necessary */
  ntotwithmass=0;
  if(indivmass[0]) ntotwithmass+=nhalo;
  if(indivmass[1]) ntotwithmass+=ndisk;
  if(indivmass[2]) ntotwithmass+=nbulge;
  if(indivmass[3]) ntotwithmass+=nstars;

  blklen = ntotwithmass*sizeof(float);
  BLKLEN;
  if(indivmass[0]) {
    for(i=0 ; i<nhalo ; i++) {
      xyz[0] = Mass(bodies + i);
      if(swapp) {
	bswap(&xyz, sizeof(float), 1);
      }
      fwrite(xyz, sizeof(float), 1, outstr);
    }
  }
  if(indivmass[1]) {
    for(i=0 ; i<ndisk ; i++) {
      xyz[0] = Mass(bodies + nhalo+ i);
      if(swapp) {
	bswap(&xyz, sizeof(float), 1);
      }
      fwrite(xyz, sizeof(float), 1, outstr);
    }
  }
  if(indivmass[2]) {
    for(i=0 ; i<nbulge ; i++) {
      xyz[0] = Mass(bodies + nhalo+ndisk+ i);
      if(swapp) {
	bswap(&xyz, sizeof(float), 1);
      }
      fwrite(xyz, sizeof(float), 1, outstr);
    }
  }
  if(indivmass[3]) {
    for(i=0 ; i<nstars ; i++) {
      xyz[0] = Mass(bodies + nhalo+ndisk+nbulge+ i);
      if(swapp) {
	bswap(&xyz, sizeof(float), 1);
      }
      fwrite(xyz, sizeof(float), 1, outstr);
    }
  }
  BLKLEN;

}

