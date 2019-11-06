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
    "in=???\n			Input file (EAGLE HDF5 format)",
    "out=???\n			Output file (NEMO snapshot)",
    "ptype=0\n                  Which particle type (0=gas, ...)",
    "region=0,10,0,10,0,10\n    Select region (xmin,xmax,ymin,ymax,zmin,zmax)",
    "group=-1\n                 Select this group (-1 means all)",
    "subgroup=-1\n              Select this subgroup (-1 means all)",
    "units=\n                   Convert to Mpc, Msol, ....",
    "center=\n                  Recentering at this x0,y0,z0",
    "boxsize=\n                 If recentering, boxsize needed for periodic grid",
    "dm=1\n                     Cheat: give the total DM mass",
    "VERSION=0.6\n		5-nov-2019",
    NULL,
};

string usage = "convert EAGLE hdf5 snapshots to NEMO snapshot format";

local int  rv_header(stream);
local void rv_data(stream, int, Body **, real *);

void nemo_main(void)
{
  stream outstr;
  real   tsnap=0, rscale, vscale, mscale, center[NDIM], boxsize, boxsize2, dm;
  string times;
  Body *btab = NULL, *p;
  int i, j, k, nbody, bits, ptype, nregion, ncenter;
  real region[6];
  float *mass, *pos, *vel;
  int *grp, *sgrp, selgrp, selsgrp, grpmin, grpmax, grpnull;
  long long *ids;
  char buf[200];
  int itype, iset, nset;

  dm = getrparam("dm");
  ptype = getiparam("ptype");
  selgrp = getiparam("group");
  selsgrp = getiparam("subgroup");

  ncenter = nemoinpr(getparam("center"),center,NDIM);
  if (ncenter > 0) {
    if (ncenter != 3) error("Need 3 values for center=");
    if (!hasvalue("boxsize")) error("Need boxsize");
    boxsize = getrparam("boxsize");
    boxsize2 = 0.5 * boxsize;
  }
  
  nregion = nemoinpr(getparam("region"),region,6);
  if (nregion == 0) {
    error("No region is not supported");
  } else if (nregion == 1) {
    region[1] = region[3] = region[5] = region[0];
    region[0] = region[2] = region[4] = 0.0;
  } else if (nregion == 2) {
    region[2] = region[4] = region[0];
    region[3] = region[5] = region[1];
  } else if (nregion != 6)
    error("region=%s needs 6 values",getparam("region"));
  if (nregion > 0)
    dprintf(0,"Region: %g,%g %g,%g %g,%g\n",
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

  if (nregion > 0)
    select_region(snap, region[0], region[1], region[2], region[3], region[4], region[5]);

  /* Find out how many particles we're going to get back */
  nbody = count_particles(snap, ptype);
  if (nbody > 0)
    dprintf(0,"Found %d particles of type %d\n",nbody,ptype);
  else
    error("No valid particles of type %d\n",ptype);

  /* Read gas particle positions */
  mass= malloc(sizeof(float)*1*nbody);
  pos = malloc(sizeof(float)*3*nbody);
  vel = malloc(sizeof(float)*3*nbody);
  grp = malloc(sizeof(int)*1*nbody);
  sgrp= malloc(sizeof(int)*1*nbody);

  if(ptype != 1) {
    if(read_dataset_float(snap, ptype, "Mass", mass, nbody) < 0)
      error("failed reading Mass: %s", get_error());
  } else {
    // @todo grab this from the MassTab attribute
    warning("No masses for DM, setting total mass to 1.0");
    for (i=0; i<nbody; i++) mass[i] = dm/nbody;
  }
  if(read_dataset_float(snap, ptype, "Coordinates", pos, nbody*3) < 0)
    error("failed reading Coordinates: %s", get_error());
  if(read_dataset_float(snap, ptype, "Velocity", vel, nbody*3) < 0)
    error("failed reading Velocity: %s", get_error());    
  if(read_dataset_int(snap, ptype, "GroupNumber", grp, nbody) < 0)
    error("failed reading GroupNumber: %s", get_error());
  if(read_dataset_int(snap, ptype, "SubGroupNumber", sgrp, nbody) < 0)
    error("failed reading SubGroupNumber: %s", get_error());    

  ids = malloc(sizeof(long long)*nbody);
  read_dataset_long_long(snap, ptype, "ParticleIDs", ids, nbody);

  grpnull = 1073741824;   /* 2^30 */
  grpmin =  grpnull;
  grpmax = -grpnull;
  for(i=0;i<nbody;i++) {
    dprintf(2,"%i %d %d %14.6f %14.6f %14.6f\n", (int) ids[i], grp[i], sgrp[i], pos[3*i+0], pos[3*i+1], pos[3*i+2]);
    if (grp[i] != grpnull) {
      if (grp[i] > grpmax) grpmax = grp[i];
      if (grp[i] < grpmin) grpmin = grp[i];
    }
  }
  dprintf(0,"GroupNumber: %d - %d\n",grpmin,grpmax);
  close_snapshot(snap);
            
  outstr = stropen(getparam("out"), "w");	/* open output file */

  btab = (body *) allocate(nbody*sizeof(body));  /* allocate body table */
  grpmin =  grpnull;
  grpmax = -grpnull;
  for (i=0, k=0, p = btab; i<nbody; i++) {
    if (selgrp  >= 0 && selgrp  !=  grp[i]) continue;
    if (sgrp[i] > grpmax) grpmax = sgrp[i];
    if (sgrp[i] < grpmin) grpmin = sgrp[i];
    if (selsgrp >= 0 && selsgrp != sgrp[i]) continue;
    Mass(p) = mass[i];
    for(j=0; j<NDIM; j++) {
      Pos(p)[j] = pos[3*i+j];
      Vel(p)[j] = vel[3*i+j];
    }
    k++;
    p++;
  }
  nbody = k;
  dprintf(0,"SubGroupNumber: %d - %d\n",grpmin,grpmax);  
  dprintf(0,"Writing %d particles (group %d, subgroup %d)\n", nbody, selgrp, selsgrp);
  if (ncenter > 0) {
    dprintf(0,"Centering on %g %g %g with boxsize %g\n",center[0],center[1],center[2],boxsize);
    for (i=0, p = btab; i<nbody; i++, p++) {
      for(j=0; j<NDIM; j++) {
	// Pos(p)[j] = fmod(Pos(p)[j] - center[j] + boxsize2, boxsize) + center[j] - boxsize2;
	Pos(p)[j] = fmod(Pos(p)[j] - center[j] + boxsize2, boxsize) - boxsize2;	
      }
    }
  }

  put_history(outstr);
  bits = (TimeBit | MassBit | PhaseSpaceBit);
  put_snap(outstr, &btab, &k, &tsnap, &bits);
  strclose(outstr);
}

