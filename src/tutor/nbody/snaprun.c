
/*
 *  SNAPRUN:     some C struct experiments
 *
   What's faster:
      particle[i].x       particle[i].y       particle[i].z
   or
      particle[i].pos[0]  particle[i].pos[1]  particle[i].pos[2]

   The (x,y,z) is a bit faster than pos[0],pos[1].pos[2]
   But looping over pos[] makes it 5x slower !!!
 
   9-jul-2024    inspired by how rebound structured particles
 */

#include <stdinc.h>
#include <getparam.h>
#include <timers.h>

string defv[] = {   /* Nemo_Keys: "key=val\n    help",       <--- required */
  "nbody=1000\n               number of bodies",
  "steps=1000\n               number of steps",
  "mode=0\n                   0: only init   1: x[3]  2:  x,y,z  3: both",
  "seed=0\n                   random number seed",
  "VERSION=0.1\n	      9-jul-2024 PJT",
  NULL,
};

string usage = "template for snapshot operations";

typedef struct _part1 {
  double x[3];
  double v[3];
} part1;

typedef struct _part2 {
  double  x,  y,  z;
  double vx, vy, vz;
} part2;
  

void nemo_main(void) /* this replaces main(argc,argv)        <--- required */
{
  int i,j,k;
  int nbody = getiparam("nbody");
  int mode = getiparam("mode");
  int nsteps = getiparam("steps");
  double dt = 1e-6;  // arbitrary
  part1 *p1 = (part1 *) allocate(nbody * sizeof(part1));
  part2 *p2 = (part2 *) allocate(nbody * sizeof(part2));

  if (mode < 0) return;

  init_timers(5);
  stamp_timers(0);
  
  // initialize when mode >=0 
  for (i=0; i<nbody; i++) {
    p1[i].x[0] = xrandom(0,1);
    p1[i].x[1] = xrandom(0,1);
    p1[i].x[2] = xrandom(0,1);
    p1[i].v[0] = xrandom(0,1);
    p1[i].v[1] = xrandom(0,1);
    p1[i].v[2] = xrandom(0,1);
    p2[i].x = xrandom(0,1);    
    p2[i].y = xrandom(0,1);    
    p2[i].z = xrandom(0,1);    
    p2[i].vx = xrandom(0,1);    
    p2[i].vy = xrandom(0,1);    
    p2[i].vz = xrandom(0,1);    
  }

  stamp_timers(1);  

  if (mode & 0x01) {
    dprintf(0,"running p1:  x[0],x[1],x[2]\n");
    for (k=0; k<nsteps; k++) {
      p1[i].x[0] += dt * p1[i].v[0];
      p1[i].x[1] += dt * p1[i].v[1];
      p1[i].x[2] += dt * p1[i].v[2];
    }
  }

  stamp_timers(2);
  
  if (mode & 0x02) {
    dprintf(0,"running p1:  x[j]\n");
    for (k=0; k<nsteps; k++) {
      for (j=0; j<3; j++) {
	p1[i].x[j] += dt * p1[i].v[j];
      }
    }
  }

  stamp_timers(3);  

  if (mode & 0x04) {
    dprintf(0,"running p2:  x,y,z\n");
    for (k=0; k<nsteps; k++) {
      p2[i].x += dt * p2[i].vx;
      p2[i].y += dt * p2[i].vy;
      p2[i].z += dt * p2[i].vz;
    }
  }

  stamp_timers(4);

  for (i=0; i<4; i++)
    printf("timer %Ld\n", diff_timers(i,i+1));
}  
