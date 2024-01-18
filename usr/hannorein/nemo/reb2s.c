/*
 * Convert a Rebound SimulationArchive to NEMO Snapshot
 *
 */

#include <nemo.h>
#include <snapshot/snapshot.h>  
#include <snapshot/body.h>
#include <snapshot/put_snap.c>
#include "rebound.h"

string defv[] = {
    "in=???\n                     input SimulationArchive",
    "out=\n                       optional output snapshot; otherwise just scan",
    "VERSION=0.4\n                17-jan-2024 PJT",
    NULL,
};

string usage="convert rebound SimulationArchive to NEMO snapshot";

void nemo_main()
{
  string infile = getparam("in");
  string outfile = getparam("out");
  struct reb_simulationarchive *sa=reb_simulationarchive_create_from_file(infile);
  if (!sa) error("Error loading Simulationarchive from file `%s`.",infile);
  dprintf(0,"SimulationArchive `%s` with %ld snapshots\n",infile, sa->nblobs);
  Body    *btab, *bp;
  int nbody = -1;
  int i, bits;
  bool Qout = hasvalue("out");
  stream  outstr;
  if (Qout) {
    outstr = stropen(outfile,"w");
    put_history(outstr);
    bits = MassBit | PhaseSpaceBit | TimeBit;
  }

  // loop over all blobs (snapshots)
  for (int64_t j=0; j < sa->nblobs; j++) {
    struct reb_simulation* r = reb_simulation_create_from_simulationarchive(sa, j);
    dprintf(0,"Snapshot %ld at time %g with nbody %d\n", j, r->t, r->N );
    if (!r) error("Error loading Simulation from Simulationarchive.\n");

    if (Qout) {
      if (nbody < 0) {
	nbody = r->N;
	btab = (Body *) allocate(nbody * sizeof(Body));
      } else if (r->N > nbody)
	error("Cannot write snapshots largers than first yet");

      for (bp=btab, i = 0; i < nbody; bp++, i++) {
	struct reb_particle p = r->particles[i];
	if (i==0)
	  dprintf(1,"%ld %d %f %f %f %f\n", j, i, p.m, p.x, p.vx, p.ax);
	else
	  dprintf(2,"%ld %d %f %f %f %f\n", j, i, p.m, p.x, p.vx, p.ax);
	Mass(bp) = p.m;
	Pos(bp)[0] = p.x;
	Pos(bp)[1] = p.y;
	Pos(bp)[2] = p.z;
	Vel(bp)[0] = p.vx;
	Vel(bp)[1] = p.vy;
	Vel(bp)[2] = p.vz;
	// @todo   .ax .ay .az if requested (need new option)
	// note rebound doesn't carry the potential
      }
      put_snap(outstr, &btab, &nbody, &r->t, &bits);  
    }

    reb_simulation_free(r);
  }
  reb_simulationarchive_free(sa);
  if (Qout)
    strclose(outstr);
}
