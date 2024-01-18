/*
 * Convert a NEMO Snapshot to a Rebound SimulationArchive
 *
 */

#include <nemo.h>
#include <filefn.h>
#include <snapshot/snapshot.h>  
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include "rebound.h"

string defv[] = {
  "in=???\n                     input snapshot",
  "out=???\n                    output SimulationArchive to be appended to",
  "VERSION=0.3\n                17-jan-2024 PJT",
  NULL,
};

string usage="convert a NEMO snapshot to a rebound SimulationArchive";

void nemo_main()
{
  string infile = getparam("in");
  string outfile = getparam("out");
  stream instr = stropen(infile, "r");
  Body    *btab = NULL, *bp;
  int i, nbody = -1, bits;
  double tsnap;
  struct reb_simulation* r = reb_simulation_create();

  if (sizeof(real) != sizeof(double))
    warning("NEMO compiled in single precision");

  if (fexist(outfile))
    warning("Existing file %s will be appended to", outfile);

  get_history(instr);
  for (;;) {                          /* loop through all snapshots */
    if (!get_tag_ok(instr, SnapShotTag))
      break;                           /* until done */
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
      dprintf (2,"Time= %f skipping, no particle data\n",tsnap);
      continue;       /* just skip - it may be diagnostics */
    }
    dprintf (1,"nbody=%d time=%f\n",nbody,tsnap);

    for (bp=btab, i = 0; i < nbody; bp++, i++) {
        struct reb_particle p = {0};
	p.m = Mass(bp);
	p.x = Pos(bp)[0];
	p.y = Pos(bp)[1];
	p.z = Pos(bp)[2];
	p.vx = Vel(bp)[0];
	p.vy = Vel(bp)[1];
	p.vz = Vel(bp)[2];
	reb_simulation_add(r, p);
	// @todo   .ax .ay .az if requested (need new option)
	// note rebound doesn't carry the potential
	if (i==0) 
	  dprintf(1,"%d %f %f %f %f\n", i, p.m, p.x, p.vx, p.ax);
	else
	  dprintf(2,"%d %f %f %f %f\n", i, p.m, p.x, p.vx, p.ax);
    }
    reb_simulation_save_to_file(r, outfile);
    
    free(btab);
    btab = NULL;
  }

}
