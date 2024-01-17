/**
 * Convert Rebound to NEMO Snapshot
 *
 */
#include <nemo.h>
#include "rebound.h"

string defv[] = {
    "in=???\n                     input SimulationArchive",
    "out=\n                       optional output snapshot; otherwise just scan",
    "VERSION=0.1\n                16-jan-2024 PJT",
    NULL,
};

string usage="convert rebound SimulationArchive to NEMO snapshot";

void nemo_main()
{
  string infile = getparam("in");
  string outfile = getparam("out");
  struct reb_simulationarchive* sa = reb_simulationarchive_create_from_file(infile);
  if (!sa) error("Error loading Simulationarchive from file `%s`.",infile);
  dprintf(0,"Loaded Rebound SimulationArchive from file `%s` with %ld snapshots\n",infile, sa->nblobs);

  if (hasvalue("out"))
    warning("Not using out=%s yet",outfile);

  // loop over all blobs (snapshots)
  for (int64_t j=0; j < sa->nblobs; j++) {
    struct reb_simulation* r = reb_simulation_create_from_simulationarchive(sa, j);
    dprintf(0,"Snapshot %ld at time %g\n", j, r->t );
    if (!r) error("Error loading Simulation from Simulationarchive.\n");

    for (int i=0; i<r->N; i++){
      struct reb_particle p = r->particles[i];
      dprintf(1,"%ld %d %f %f %f %f\n", j, i, p.m, p.x, p.vx, p.ax);
    }

    reb_simulation_free(r);
  }
  reb_simulationarchive_free(sa);
}
