/*
 * gadgetsnap:  convert the Gadget nbody data format in nemo's 
 *              snapshot(5NEMO) format
 *
 *  2-jun-2003         based on their read_snapshot.c template example    - Peter Teuben
 * 17-mar-2006  V0.2   using header=
 *
 *
 */

#include <nemo.h>
#include <vectmath.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n       Input file (gadget format)",
    "out=???\n      Output file (snapshot format)",
    "swap=f\n       Force swaping bytes for non-native machines",
    "header=\n      Header size of unfio (4 or 8 or pre-configured)",
    "VERSION=0.2\n  17-mar-06 PJT",
    NULL,
};

string usage="convert gadget files to snapshot format";

#define HEADER_SIZE  256

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
  /* 
   *  the remainder is a 256byte space filler, and as the structure above
   *  here changes over time, will need to be adjusted 
   *  we're also making the usual assumption:
   *  int=4, double=8 and no padding between them
   */
  char     fill[HEADER_SIZE - 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8]; 
} header1;


int     NumPart, Ngas;
bool    Qswap;
int     header;

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *P;

int *Id;

double  Time, Redshift;

local int allocate_memory(void);
local int reordering(void);
local int do_what_you_want(stream);
local int unit_conversion(void);
local int load_snapshot(char *fname, int files);
                                                                                


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
void nemo_main()
{
  string input_fname;
  int  type, snapshot_number, files;
  stream outstr;

  if (hasvalue("header"))
    header = getiparam("header");
  else
    header = UNFIO_HDR_SIZE;

  if (sizeof(header1) != HEADER_SIZE)
      error("sizeof(io_header) = %d  != HEADER_SIZE=%d",
	    sizeof(header1),HEADER_SIZE);
  if (sizeof(int) != 4)    error("sizeof(int)    = %d (we assumed 4)",sizeof(int));
  if (sizeof(double) != 8) error("sizeof(double) = %d (we assumed 8)",sizeof(double));

  Qswap = getbparam("swap");

#if 0
  sprintf(path, "/home/vspringe/LCDM");
  sprintf(basename, "snapshot");
  snapshot_number= 6;                    /* number of snapshot */


  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
#endif
  input_fname = getparam("in");
  files=1;                               /* number of files per snapshot */
  load_snapshot(input_fname, files);


  reordering();  /* call this routine only if your ID's are set properly */

  unit_conversion();  /* optional stuff */


  outstr = stropen(getparam("out"),"w");
  put_history(outstr);
  do_what_you_want(outstr);
}





/* 
 *   write NEMO snapshot  (for now only all mass,pos,vel)
 */
local int do_what_you_want(stream outstr)
{
  Body *bp, *btab;
  real tsnap = Time;
  int i, bits = (TimeBit | MassBit | PhaseSpaceBit);


  btab = allocate(NumPart*sizeof(Body));

  for (i=1, bp=btab; i<=NumPart; i++, bp++) {
    SETV(Pos(bp), P[i].Pos);
    SETV(Vel(bp), P[i].Vel);
    Mass(bp) = P[i].Mass;
  }
  put_snap(outstr,&btab, &NumPart, &tsnap, &bits);

  free(btab);

}





/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
local int unit_conversion(void)
{
  double GRAVITY, BOLTZMANN, PROTONMASS;
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;  
  double G, Xh, HubbleParam;

  int i;
  double MeanWeight, u, gamma;

  /* physical constants in cgs units */
  GRAVITY   = 6.672e-8;
  BOLTZMANN = 1.3806e-16;
  PROTONMASS = 1.6726e-24;

  /* internal unit system of the code */
  UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h */
  UnitMass_in_g= 1.989e43;         /*  code mass unit in g/h */
  UnitVelocity_in_cm_per_s= 1.0e5;

  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
  UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
  UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

  G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);


  Xh= 0.76;  /* mass fraction of hydrogen */
  HubbleParam= 0.65;


  for(i=1; i<=NumPart; i++)
    {
      if(P[i].Type==0)  /* gas particle */
	{
	  MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

	  gamma= 5.0/3;
	 
	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	}
    }
}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
local int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200], dummy[8];
  int    i,j,k,l,ntot_withmasses;
  int    t,n,off,pc,pc_new,pc_sph;

#define SKIP fread(dummy, header, 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new) {
    if(files>1)
      sprintf(buf,"%s.%d",fname,i);
    else
      sprintf(buf,"%s",fname);
    
    if(!(fd=fopen(buf,"r")))
      error("can't open file `%s`",buf);

    dprintf(0,"reading `%s' ...\n",buf);

    /* this relies on the fact that fortran binary I/O
       dumps the record size before and after the actual
       data, the header in this case
       We could also use the unfio() routines, to
       simplify the use of the SKIP macro
    */
    SKIP;
    fread(&header1, sizeof(header1), 1, fd);
    SKIP;
    if (Qswap) {
      bswap(header1.npart,sizeof(int),6);
      bswap(header1.mass,sizeof(double),6);
      bswap(&header1.time,sizeof(double),1);
      bswap(&header1.redshift,sizeof(double),1);
      bswap(&header1.flag_sfr,sizeof(int),1);
      bswap(&header1.flag_feedback,sizeof(int),1);
      bswap(header1.npartTotal,sizeof(int),6);
      bswap(&header1.flag_cooling,sizeof(int),1);
      bswap(&header1.num_files,sizeof(int),1);
      bswap(&header1.BoxSize,sizeof(double),1);
      bswap(&header1.Omega0,sizeof(double),1);
      bswap(&header1.OmegaLambda,sizeof(double),1);
      bswap(&header1.HubbleParam,sizeof(double),1);
    }


    dprintf(0,"Header: \nnpart: ");
    for (l=0; l<6; l++) dprintf(0," %d",header1.npart[l]);
    dprintf(0,"\nmass: ");
    for (l=0; l<6; l++) dprintf(0," %g",header1.mass[l]);
    dprintf(0,"\ntime: %g",header1.time);
    dprintf(0,"\nredshift: %g",header1.redshift);
    dprintf(0,"\nflag_sfr: %d",header1.flag_sfr);
    dprintf(0,"\nflag_feedback: %d",header1.flag_feedback);
    dprintf(0,"\nnpartTotal: ");
    for (l=0; l<6; l++) dprintf(0," %d",header1.npartTotal[l]);
    dprintf(0,"\nflag_cooling: %d",header1.flag_cooling);
    dprintf(0,"\nnum_files: %d",header1.num_files);
    dprintf(0,"\nBoxSize: %g",header1.BoxSize);
    dprintf(0,"\nOmega0: %g",header1.Omega0);
    dprintf(0,"\nOmegaLambda: %g",header1.OmegaLambda);
    dprintf(0,"\nHubbleParam: %g",header1.HubbleParam);
    dprintf(0,"\n");

    if(files==1) {
      for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	NumPart+= header1.npart[k];
      Ngas= header1.npart[0];
    } else {
      for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	NumPart+= header1.npartTotal[k];
      Ngas= header1.npartTotal[0];
    }
    
    for(k=0, ntot_withmasses=0; k<5; k++) {
      if(header1.mass[k]==0)
	ntot_withmasses+= header1.npart[k];
    }
    
    if(i==0)
      allocate_memory();
    
    SKIP;
    for(k=0,pc_new=pc;k<6;k++) {
      for(n=0;n<header1.npart[k];n++) {
	fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	if (Qswap) bswap(&P[pc_new].Pos[0], sizeof(float), 3);
	pc_new++;
      }
    }
    SKIP;
    
    SKIP;
    for(k=0,pc_new=pc;k<6;k++) {
      for(n=0;n<header1.npart[k];n++) {
	fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	if (Qswap) bswap(&P[pc_new].Vel[0], sizeof(float), 3);
	pc_new++;
      }
    }
    SKIP;
    

    SKIP;
    for(k=0,pc_new=pc;k<6;k++) {
      for(n=0;n<header1.npart[k];n++) {
	fread(&Id[pc_new], sizeof(int), 1, fd);
	if (Qswap) bswap(&Id[pc_new], sizeof(int), 1);
	pc_new++;
      }
    }
    SKIP;
    

    if(ntot_withmasses>0)
      SKIP;
    for(k=0, pc_new=pc; k<6; k++) {
      for(n=0;n<header1.npart[k];n++) {
	P[pc_new].Type=k;
	
	if(header1.mass[k]==0) {
	  fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	  if (Qswap) bswap(&P[pc_new].Mass, sizeof(float), 1);
	} else
	  P[pc_new].Mass= header1.mass[k];
	pc_new++;
      }
    }
    if(ntot_withmasses>0)
      SKIP;
    
    
    if(header1.npart[0]>0) {
      SKIP;
      for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
	fread(&P[pc_sph].U, sizeof(float), 1, fd);
	pc_sph++;
      }
      SKIP;
      
      SKIP;
      for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
	fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	if (Qswap) bswap(&P[pc_sph].Rho, sizeof(float), 1);
	pc_sph++;
      }
      SKIP;
      
      if(header1.flag_cooling) {
	SKIP;
	for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
	  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
	  if (Qswap) bswap(&P[pc_sph].Ne, sizeof(float), 1);
	  pc_sph++;
	}
	SKIP;
      } else
	for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
	  P[pc_sph].Ne= 1.0;
	  pc_sph++;
	}
    }

    fclose(fd);
  } /* reading all file(s) */
  
  Time= header1.time;
  Redshift= header1.time;
}


/* 
 *   allocate memory for the particle data.
 */

local int allocate_memory(void)
{
  dprintf(1,"NumPart=%d : allocating memory...\n",NumPart);
  P=allocate(NumPart*sizeof(struct particle_data));
  P--;   /* start with offset 1 */
  Id=malloc(NumPart*sizeof(int));
  Id--;   /* start with offset 1 */
}



/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
local int reordering(void)
{
  int i,j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  dprintf(1,"reordering....\n");

  for(i=1; i<=NumPart; i++) {
    if(Id[i] != i) {
      psource= P[i];
      idsource=Id[i];
      dest=Id[i];
       
      do {
	psave= P[dest];
	idsave=Id[dest];

	P[dest]= psource;
	Id[dest]= idsource;
	
	if(dest == i) 
	  break;

	psource= psave;
	idsource=idsave;
	
	dest=idsource;
      }
      while(1);
    }
  }

  Id++;   
  free(Id);

  dprintf(1,"space for particle ID freed\n");
}






  











