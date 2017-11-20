/*
 * SNAPFOLD.C: fold radial lines to create a conical shape
 *
 *      20-nov-2017 V0.1 hack                                           PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {	
    "in=???\n	           Input snapshot file",
    "out=???\n             Output snapshot file",
    "theta=\n              Angle (degrees) to fold from XY into Z [-90..90]",
    "select=pos,vel,acc\n  Select the vectors for rotation",
    "VERSION=0.1\n         20-nov-2017 PJT",
    NULL,
};

string usage = "fold a snapshot";

string cvsid="$Id$";

#define DEG2RAD PI/180.0



void rotate(real angle,Body *btab,int nbody);

void nemo_main()
{
    stream instr, outstr;
    char *op, *order;
    Body *btab = NULL;
    int i, nop, nopt, nbody, bits;
    real tsnap, tscale, *ang, theta;
    bool   invert, need_hist = TRUE;
    string rotvects = getparam("select");
    bool Qpos, Qvel, Qacc;

    instr = stropen(getparam("in"), "r");           /* open input file */
    theta = getrparam("theta");
    outstr = stropen(getparam("out"), "w");

    get_history(instr);
    while (get_tag_ok(instr, SnapShotTag)) {
      get_snap(instr, &btab, &nbody, &tsnap, &bits);
      if (bits & PhaseSpaceBit) {
	dprintf(1,"Processing time %g bits=0x%x\n",tsnap,bits);
	rotate(theta, btab, nbody);
	if (need_hist) {
	  put_history(outstr);
	  need_hist = FALSE;
	}
	put_snap(outstr, &btab, &nbody, &tsnap, &bits);
      } /* if bits */
    } /* while */
    strclose(instr);
    strclose(outstr);
}

void rotate(real angle,Body *btab,int nbody)
{
    real x,y,vx,vy,r,vr;
    real sint, cost;
    int  i;

    /*
      1)  rotate the positions in the XY-Z plane
      2)  keep the VT
      3)  rotate the VR along the radius
     */
    cost = cos(angle * DEG2RAD);
    sint = sin(angle * DEG2RAD);
    
    for (i = 0; i < nbody; i++) {
      x = Pos(&btab[i])[0];
      y = Pos(&btab[i])[1];
      
      r = sqrt(x*x+y*y);
      if (r==0) continue;
      
      vx = Vel(&btab[i])[0];
      vy = Vel(&btab[i])[1];

      vr = (x*vx + y*vy) / r;

      Pos(&btab[i])[0] = x * cost;
      Pos(&btab[i])[1] = y * cost;
      Pos(&btab[i])[2] = r * sint;

      Vel(&btab[i])[0] += vr * cost * (x/r);
      Vel(&btab[i])[1] += vr * cost * (y/r);
      Vel(&btab[i])[2]  = vr * sint;
    }
}
