/* mkgalaxy2.c : a NEMOfied version of the 'galaxy' program used in
 *               in Xscreensaver 
 *
 * Originally done by Uli Siegmund <uli@wombat.okapi.sub.org> on Amiga
 *   for EGS in Cluster
 * Port from Cluster/EGS to C/Intuition by Harald Backert
 * Port to X11 and incorporation into xlockmore by Hubert Feyrer
 *   <hubert.feyrer@rz.uni-regensburg.de>
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose and without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and that
 * both that copyright notice and this permission notice appear in
 * supporting documentation.
 *
 * This file is provided AS IS with no warranties of any kind.  The author
 * shall have no liability with respect to the infringement of copyrights,
 * trade secrets or any patents by this file or any part thereof.  In no
 * event will the author be liable for any lost revenue or profits or
 * other special, indirect and consequential damages.
 *
 * Revision History:
 * 26-Aug-00: robert.nagtegaal@phil.uu.nl and roland@tschai.demon.nl:
 *            various improvements
 * 10-May-97: jwz@jwz.org: turned into a standalone program.
 * 18-Apr-97: Memory leak fixed by Tom Schmidt <tschmidt@micron.com>
 * 07-Apr-97: Modified by Dave Mitchell <davem@magnet.com>
 * 23-Oct-94: Modified by David Bagley <bagleyd@bigfoot.com>
 *  random star sizes
 *  colors change depending on velocity
 * 10-Oct-94: Add colors by Hubert Feyer
 * 30-Sep-94: Initial port by Hubert Feyer
 * 09-Mar-94: VMS can generate a random number 0.0 which results in a
 *            division by zero, corrected by Jouk Jansen
 *            <joukj@crys.chem.uva.nl>
 *
 * 13-may-02: (V2.0) NEMO-fied                   by Peter Teuben
 *  7-feb-2021:    renamed to mkgalaxy2 to avoid colliding w/ falcON's mkgalaxy script
 */

#include <nemo.h>
#include <mathfns.h>

#include <snapshot/snapshot.h>  
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string  defv[] = { 
  "out=???\n		      Output file name",
  "seed=0\n               Seed for the random number generator",
  "time=0.0\n             Time at which snapshot is taken",
  "ngal=3\n               Number of galaxies",
  "nstars=1000\n          Number of stars per galaxy",
  "ncolors=4\n            Number of colors",
  "galsize=0.15,0.25\n    Range in galaxy sizes",
  "size=2\n               Size of the universe",
  "niters=7500\n          Number of iterations",
  "zerocm=t\n             Centrate snapshot (t/f)?",
  "headline=\n	          Verbiage for output",
  "VERSION=3.0\n          7-feb-2021 PJT",
  NULL,
};

string usage="Create a bunch of spinning disk galaxies";


#define MINSIZE       1
#define MINGALAXIES    2
#define MAX_STARS    3000
#define MAX_IDELTAT    50
/* These come originally from the Cluster-version */
#define DEFAULT_GALAXIES  3
#define DEFAULT_STARS    1000
#define DEFAULT_HITITERATIONS  7500
#define DEFAULT_IDELTAT    200 /* 0.02 */
#define DELTAT (MAX_IDELTAT * 0.0001)

#define GALAXYRANGESIZE  0.1
#define GALAXYMINSIZE  0.15

#define NUMCOLORS 256

#define COLORBASE  16
  /* Colors for stars start here */
#define COLORSTEP  (NUMCOLORS/COLORBASE) /* NUMCOLORS / COLORBASE colors per galaxy */

typedef struct {                    /* massless disk particles in a galaxy */
  double      pos[3], vel[3];
} Star;

typedef struct {                    /* galaxy */
  double      mass;
  int         nstars;
  Star       *stars;
  double      pos[3], vel[3];
  int         galcol;
} Galaxy;

typedef struct {                    /* universe */
  double      mat[3][3]; 
  double      size;      
  Galaxy     *galaxies;  
  int         ngalaxies; 
} Universe;

static Universe universes[1];

void nemo_main()
{
  Universe   *u = &universes[0];
  int         i, j, seed;
  double      w1, w2, size;
  double      d, v, w, h;
  int ncolors = getiparam("ncolors");
  int nstars  = getiparam("nstars");
  int niters  = getiparam("niters");
  int ngal    = getiparam("ngal");
  int nbody = 0;
  Body *btab, *bp;

  warning("Program not completed/tested, use at own risk");

  seed = init_xrandom(getparam("seed"));
  size = getdparam("size");
  u->ngalaxies = ngal;
  u->galaxies = (Galaxy *) allocate(u->ngalaxies * sizeof (Galaxy));

  for (i = 0; i < u->ngalaxies; ++i) {
	Galaxy     *g = &u->galaxies[i];
	double      sinw1, sinw2, cosw1, cosw2;

	nbody++;             // add this galaxy
	g->galcol = ncolors;     // NRAND(COLORBASE - 2);
	if (g->galcol > 1)
	  g->galcol += 2; /* Mult 8; 16..31 no green stars */
	/* Galaxies still may have some green stars but are not all green. */

	g->nstars = nstars;    // (NRAND(MAX_STARS / 2)) + MAX_STARS / 2;
	g->stars = (Star *) allocate(g->nstars * sizeof (Star));

	w1 = xrandom(0.0, TWO_PI);
	w2 = xrandom(0.0, TWO_PI);
	sinw1 = sin(w1);
	sinw2 = sin(w2);
	cosw1 = cos(w1);
	cosw2 = cos(w2);
	
	u->mat[0][0] = cosw2;
	u->mat[0][1] = -sinw1 * sinw2;
	u->mat[0][2] = cosw1 * sinw2;
	u->mat[1][0] = 0.0;
	u->mat[1][1] = cosw1;
	u->mat[1][2] = sinw1;
	u->mat[2][0] = -sinw2;
	u->mat[2][1] = -sinw1 * cosw2;
	u->mat[2][2] = cosw1 * cosw2;
#if 0	
	g->vel[0] = xrandom(-1.0,1.0);
	g->vel[1] = xrandom(-1.0,1.0);
	g->vel[2] = xrandom(-1.0,1.0);
	g->pos[0] = -g->vel[0] * DELTAT * niters + xrandom(-0.5,0.5);
	g->pos[1] = -g->vel[1] * DELTAT * niters + xrandom(-0.5,0.5);
	g->pos[2] = -g->vel[2] * DELTAT * niters + xrandom(-0.5,0.5);
#else
	g->pos[0] = xrandom(-size,size);
	g->pos[1] = xrandom(-size,size);
	g->pos[2] = xrandom(-size,size);
	g->vel[0] = xrandom(-1.0,1.0);
	g->vel[1] = xrandom(-1.0,1.0);
	g->vel[2] = xrandom(-1.0,1.0);
#endif
	g->mass = xrandom(0.0,1.0);

	u->size = xrandom(GALAXYMINSIZE, GALAXYMINSIZE+GALAXYRANGESIZE);


	nemo_dprintf(0,"Galaxy: mass: %g   %g size %d stars\n",g->mass,u->size,g->nstars);
	nemo_dprintf(0,"Galaxy: pos: (%g,%g,%g)\n",g->pos[0], g->pos[1], g->pos[2]);
	nemo_dprintf(0,"Galaxy: vel: (%g,%g,%g)\n",g->vel[0], g->vel[1], g->vel[2]);

	nbody += g->nstars;            // add all the stars in this galaxy
	for (j = 0; j < g->nstars; j++) {
	  Star       *s = &g->stars[j];
	  double      sinw, cosw;

	  w = xrandom(0.0, TWO_PI);
	  sinw = sin(w);
	  cosw = cos(w);
	  d = xrandom(0.0, u->size);
	  h = xrandom(0.0,1.0) * exp(-2.0 * (d / u->size)) / 5.0 * u->size;
	  if (xrandom(0.0,1.0) < 0.5)
		h = -h;
	  s->pos[0] = u->mat[0][0] * d * cosw + u->mat[1][0] * d * sinw + u->mat[2][0] * h + g->pos[0];
	  s->pos[1] = u->mat[0][1] * d * cosw + u->mat[1][1] * d * sinw + u->mat[2][1] * h + g->pos[1];
	  s->pos[2] = u->mat[0][2] * d * cosw + u->mat[1][2] * d * sinw + u->mat[2][2] * h + g->pos[2];

	  v = sqrt(g->mass / sqrt(d * d + h * h));
	  s->vel[0] = -u->mat[0][0] * v * sinw + u->mat[1][0] * v * cosw + g->vel[0];
	  s->vel[1] = -u->mat[0][1] * v * sinw + u->mat[1][1] * v * cosw + g->vel[1];
	  s->vel[2] = -u->mat[0][2] * v * sinw + u->mat[1][2] * v * cosw + g->vel[2];

	  s->vel[0] *= DELTAT;
	  s->vel[1] *= DELTAT;
	  s->vel[2] *= DELTAT;
	}
  }
  if (hasvalue("out")) {
	stream outstr = stropen(getparam("out"), "w");
	real snap_time = getdparam("time");
	int bits = (MassBit | PhaseSpaceBit | TimeBit);

	put_history (outstr);     
	btab = (Body *) allocate (nbody * sizeof(Body));
	for (i = 0, bp=btab; i < u->ngalaxies; ++i) {
	  Galaxy     *g = &u->galaxies[i];

	  Mass(bp) = g->mass;
	  SETV(Pos(bp),g->pos);
	  SETV(Vel(bp),g->vel);
	  for (j = 0; j < g->nstars; ++j) {
		Star       *s = &g->stars[j];
		bp++;
		Mass(bp) = 0.0;
		SETV(Pos(bp),s->pos);
		SETV(Vel(bp),s->vel);
	  }
	  bp++;
	}
	put_snap (outstr, &btab, &nbody, &snap_time, &bits);
	strclose(outstr);
	free(btab);
  } 

}
