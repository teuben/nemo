/*
 * SNAPFOLD.C: fold radial lines to create a conical shape
 *
 *      20-nov-2017 V0.1 hack at ESO                                 PJT
 *      22-mar-2018 V0.3 once again at ESO, hack-2 filling the cone  PJT
 *      28-jul-2020 V0.4 add view=                                   PJT
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
    "theta=0\n              Angle (degrees) to fold from XY into Z [-90..90]",
    "select=pos,vel,acc\n  Select the vectors for rotation",
    "fill=f\n              Filling the cone uniformly?",
    "view=0\n              view all (0), back (-1) or front (1)",
    "VERSION=0.4\n         28-jul-2020 PJT",
    NULL,
};

string usage = "fold a snapshot along the Z axis into a cone";

string cvsid="$Id$";

#define DEG2RAD PI/180.0



void rotate(real angle, bool fill, int view, Body *btab,int nbody);

real my_fun(real x)
{
  return cos(x * DEG2RAD);
}

void nemo_main()
{
    stream instr, outstr;
    char *op, *order;
    Body *btab = NULL;
    int i, nop, nopt, nbody, bits;
    real tsnap, tscale, *ang, theta;
    bool   invert, fill, need_hist = TRUE;
    string rotvects = getparam("select");
    bool Qpos, Qvel, Qacc;
    int view = getiparam("view");

    instr = stropen(getparam("in"), "r");           /* open input file */
    theta = getrparam("theta");
    outstr = stropen(getparam("out"), "w");
    fill = getbparam("fill");

    get_history(instr);
    while (get_tag_ok(instr, SnapShotTag)) {
      get_snap(instr, &btab, &nbody, &tsnap, &bits);
      if (bits & PhaseSpaceBit) {
	dprintf(1,"Processing time %g bits=0x%x\n",tsnap,bits);
	rotate(theta, fill, view, btab, nbody);
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

void rotate(real angle, bool fill, int view, Body *btab,int nbody)
{
    real x,y,vx,vy,r,vr;
    real t, sint, cost, phi, cosp, sinp;
    int  i;

    /*
     * adjust the view ?
     */
    if (view != 0) {
      for (i = 0; i < nbody; i++) {
	x = Pos(&btab[i])[0];
	y = Pos(&btab[i])[1];
	vx = Vel(&btab[i])[0];
	vy = Vel(&btab[i])[1];
      
	r = sqrt(x*x+y*y);
	
	if (r==0) continue;

	if (view > 0) {
	  phi = atan2(x,-y) / (1+view);
	  cosp = cos(phi);
	  sinp = sin(phi);
	  Pos(&btab[i])[0] =   x * cosp +  y * sinp;
	  Pos(&btab[i])[1] =  -x * sinp +  y * cosp;
	  Vel(&btab[i])[0] =  vx * cosp + vy * sinp;
	  Vel(&btab[i])[1] = -vx * sinp + vy * cosp;
	  
	} else {
	  phi = atan2(x,y) / (1-view);
	  cosp = cos(phi);
	  sinp = sin(phi);	  
	  Pos(&btab[i])[0] =   x * cosp -  y * sinp;
	  Pos(&btab[i])[1] =   x * sinp +  y * cosp;
	  Vel(&btab[i])[0] =  vx * cosp - vy * sinp;
	  Vel(&btab[i])[1] =  vx * sinp + vy * cosp;
	  
	}

#if 0	
	Pos(&btab[i])[0] = x;
	Pos(&btab[i])[1] = y;
	

	vr = (x*vx + y*vy) / r;
	
	Vel(&btab[i])[0] += vr * cost * (x/r);
	Vel(&btab[i])[1] += vr * cost * (y/r);
	Vel(&btab[i])[2]  = vr * sint;
#endif	
      }
    }

    
    /*
      1)  rotate the positions in the XY-Z plane
      2)  keep the VT
      3)  rotate the VR along the radius
     */
    cost = cos(angle * DEG2RAD);
    sint = sin(angle * DEG2RAD);
    
    for (i = 0; i < nbody; i++) {
      if (fill) {
	// t = xrandom(angle, 90.0);
	t = frandom(angle, 90.0, my_fun);
	cost = cos(t * DEG2RAD);
	sint = sin(t * DEG2RAD);
      }
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
