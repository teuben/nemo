/*
 * SNAPSHELL.C: compute various diagnostic properties in a set of shells
 *              
 * 
 *     13-nov-01     V1.0   derived from snapkinem                          PJT
 *     16-nov-01            implemente more....
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <moment.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {	
    "in=???\n			 Input file name (snapshot)",
    "radii=???\n                 (normalized) radii",
    "pvar=vr\n                   Variables to print statistics of",
    "svar=\n                     sorting variable if shells sorted by an expression",
    "weight=1\n			 weighting for particles",
    "axes=1,1,1\n                X,Y,Z axes for spatial spheroidal normalization",
    "stats=mean,disp,n\n         Statistics to print (mean,disp,skew,kurt,min,max,median,n)",
    "normalized=t\n              Use normalized radii is svar= is used?",
    "VERSION=1.0a\n		 16-nov-01 PJT",
    NULL,
};

string usage="compute various diagnostics on a set of shells";


Body *btab = NULL;		/* pointer to array of bodies		    */
int nbody;			/* number of bodies in array		    */
real tsnap = 0.0;		/* time associated with data		    */

rproc weight;			/* weighting function for bodies	    */
rproc pvar;
rproc svar;

real rcut = -1.0;
vector axes;                    /* normalization radii for shells           */
bool Qrad;

#define MAXRAD 10000
int nrad;
real radii[MAXRAD];


nemo_main()
{
    stream instr;
    rproc btrtrans();
    int bits, ndim;

    instr = stropen(getparam("in"), "r");
    nrad = nemoinpd(getparam("radii"),radii,MAXRAD);
    if (nrad<0) error("Parsing rad=");
    get_history(instr);
    weight = btrtrans(getparam("weight"));
    pvar = btrtrans(getparam("pvar"));
    if (hasvalue("svar")) {
      svar = btrtrans(getparam("svar"));
      Qrad = FALSE;
    } else {
      Qrad = TRUE;
    }
    ndim = nemoinpd(getparam("axes"),axes,3);
    if (ndim != NDIM) error("Not enough values for axes=");
    dprintf(1,"Axes: %g %g %g\n",axes[0],axes[1],axes[2]);
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if (bits & PhaseSpaceBit) {
      snapkinem();
    }
}

snapkinem()
{
    reshape();
    findmoment();
}

int n_tot;			/* number of bodies with positive weight    */

vector cm_pos;			/* rough center of mass position	    */

int n_sum;			/* number of bodies contributing	    */

real w_sum;			/* sum of body weights			    */

vector w_pos;			/* weighted center of mass position	    */
vector w_vel;			/* weighted center of mass velocity	    */
vector w_jvec;			/* specific angular momentum		    */

matrix w_qpole;			/* weighted quadrupole moment		    */
matrix w_keten;			/* weighted kinetic energy tensor	    */

reshape()
{
    int i;
    Body *b;
    vector tmpv;

    for (i = 0, b = btab; i < nbody; i++, b++) {
      MULVV(tmpv, Pos(b), axes);
      SETV(Pos(b),tmpv);
      /* printvec("pos: ",Pos(b)); */
    }
}


findmoment()
{
  int i, j, irad;
  Body *b;
  real rad;
  vector tmpv, pos_b, vel_b;
  matrix tmpm;
  Moment mr, mq;

  ini_moment(&mr,2);
  ini_moment(&mq,2);
  
  for (i = 0, b = btab, irad=0; i < nbody; ) {
    reset_moment(&mr);
    reset_moment(&mq);
    if (Qrad) {
      for (j=0; i < nbody ; b++, i++, j++) {
	rad = absv(Pos(b));
	dprintf(1,"%g checking %d[%g %g]\n",
		rad,irad,radii[irad],radii[irad+1]);
	if (rad >= radii[irad] && rad < radii[irad+1]) {
	  accum_moment(&mr, rad, 1.0);
	  accum_moment(&mq, (pvar)(b, tsnap, i), 1.0);
	} else if (rad < radii[irad]) {
	  dprintf(1,"Skipping %d\n",j);
	} else {
	  irad++;
	  break;
	}
      }
    } else {
      error("Not implemented yet");
    }
    if (n_moment(&mr))
      printf("%g %g %d : %g %g\n",
	   mean_moment(&mr), sigma_moment(&mr),
	   n_moment(&mr),
	   mean_moment(&mq), sigma_moment(&mq)
	   );
    if (irad >= nrad) break;
  }
}

printvec(name, vec)
string name;
vector vec;
{
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n",
	   name, absv(vec), vec[0], vec[1], vec[2]);
}

