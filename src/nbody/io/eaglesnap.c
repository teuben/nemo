/*
 *  EAGLESNAP: convert hdf5 files from the EAGLE simulation to NEMO snapshot
 *
 *
 *	2-nov-2019	V0.1	Created from read_eagle's example.c      PJT
 *	
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <stdlib.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

#include <read_eagle.h>


string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (RV format)",
    "out=???\n			Output file (snapshot)",
    "ptype=0\n                  Which particle type (0=gas, ...)",
    "region=0,10,0,10,0,10\n    Select region (xmin,xmax,ymin,ymax,zmin,zmax)",
    "VERSION=0.1\n		2-nov-2019",
    NULL,
};

string usage = "convert EAGLE hdf5 snapshots to NEMO snapshot format";

local int  rv_header(stream);
local void rv_data(stream, int, Body **, real *);

void nemo_main(void)
{
  stream instr, outstr;
  real   tsnap;
  string times;
  Body *btab = NULL, *p;
  int i, j, nbody, bits, ptype, nregion;
  real region[6];
  float *mass, *pos, *vel;
  long long *ids;
  char buf[200];
  int itype, iset, nset;

  ptype = getiparam("ptype");
  nregion = nemoinpr(getparam("region"),region,6);
  if (nregion == 1) {
    region[1] = region[3] = region[5] = region[0];
    region[0] = region[2] = region[4] = 0.0;
  } else if (nregion == 2) {
    region[2] = region[4] = region[0];
    region[3] = region[5] = region[1];
  } else if (nregion != 6)
    error("region=%s needs 6 values",getparam("region"));
  printf("Region: %g,%g %g,%g %g,%g\n",
	 region[0], region[1], region[2], region[3], region[4], region[5]);
  
  EagleSnapshot *snap = open_snapshot(getparam("in"));
  if (!snap) error("Cannot open EagleSnapshot %s", getparam("in"));
  
  for(itype=0; itype<6; itype++) {
    nset = get_dataset_count(snap, itype);
    for(iset=0; iset<nset; iset++) {
      get_dataset_name(snap, itype, iset, buf, 200);
      dprintf(1,"# Type %d has dataset: %s\n", itype, buf);
    }
  }

  /* Choose region to read */
  select_region(snap, region[0], region[1], region[2], region[3], region[4], region[5]);

  /* Find out how many particles we're going to get back */
  nbody = count_particles(snap, ptype);
  if (nbody > 0)
    printf("Found %d particles of type %d\n",nbody,ptype);
  else
    error("No valid particles of type %d\n",ptype);

  /* Read gas particle positions */
  mass= malloc(sizeof(float)*1*nbody);
  pos = malloc(sizeof(float)*3*nbody);
  vel = malloc(sizeof(float)*3*nbody);
  if(read_dataset_float(snap, ptype, "Coordinates", pos, nbody*3) < 0) {
      printf("read_dataset failed!\n");
      printf("Reason: %s\n", get_error());
      exit(1);
  }
  if(read_dataset_float(snap, ptype, "Velocity", vel, nbody*3) < 0) {
      printf("read_dataset failed!\n");
      printf("Reason: %s\n", get_error());
      exit(1);
  }

  /* Read particle IDs */
  ids = malloc(sizeof(long long)*nbody);
  read_dataset_long_long(snap, ptype, "ParticleIDs", ids, nbody);

  /* Write out the particles */
  for(i=0;i<nbody;i++)
    dprintf(2,"%i %14.6f %14.6f %14.6f\n", (int) ids[i], pos[3*i+0], pos[3*i+1], pos[3*i+2]);
  
  close_snapshot(snap);

            
  outstr = stropen(getparam("out"), "w");	/* open output file */

  btab = (body *) allocate(nbody*sizeof(body));  /* allocate body table */
  for (i=0, p = btab; i<nbody; i++, p++) {      
    Mass(p) = mass[i];
    for(j=0; j<NDIM; j++) {
      Pos(p)[j] = pos[3*i+j];
      Vel(p)[j] = vel[3*i+j];
    }
  }

  put_history(outstr);
  bits = (TimeBit | MassBit | PhaseSpaceBit);
  put_snap(outstr, &btab, &nbody, &tsnap, &bits);
  strclose(outstr);
}

