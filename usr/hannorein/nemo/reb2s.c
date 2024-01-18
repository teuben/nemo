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
    "verbose=f\n                  show more info from the simulation",
    "VERSION=0.6\n                18-jan-2024 PJT",
    NULL,
};

string usage="convert rebound SimulationArchive to NEMO snapshot";

local void reb_verbose(struct reb_simulation* r);

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
  bool Qverbose = getbparam("verbose");
  stream  outstr;
  
  if (sizeof(real) != sizeof(double))
    warning("NEMO compiled in single precision");
  
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
    if (Qverbose) reb_verbose(r);

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

/*
 *   something like this should be in rebound
 */

local void reb_verbose(struct reb_simulation* r)
{
  printf("\n");
  printf("Time: %g\n", r->t);
  printf("N:    %d\n", r->N);
  printf("G:    %g\n", r->G);
  printf("Softening:          %g\n", r->softening);
  printf("dt_last_done:       %g\n", r->dt_last_done);
  printf("steps_done:         %ld\n", r->steps_done);
  printf("opening_angle2:     %g\n",r->opening_angle2);
  printf("output_timing_last: %g\n",r->output_timing_last);
  //double walltime;                // Cumulative walltime of entire integration.
  //double walltime_last_step;      // Wall time of last step.
  //double walltime_last_steps;     // Average wall time of last step (updated every 0.1s).
  //double walltime_last_steps_sum;
  //int walltime_last_steps_N;
    
  double e = reb_simulation_energy(r);
  printf("Total Energy:       %g\n",e);    

  struct reb_vec3d j = reb_simulation_angular_momentum(r);
  printf("Angular Momentum:   %g %g %g\n", j.x, j.y, j.z);

  struct reb_particle com = reb_simulation_com(r);
  printf("Center of Mass:     %g %g %g  %g %g %g\n", com.x, com.y, com.z, com.vx, com.vy, com.vz);

  
}
